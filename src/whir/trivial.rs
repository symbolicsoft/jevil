//! HVZK base case (Construction 7.2 of eprint 2026/391).
//!
//! Closes the WHIR recursion with a γ-batched mask reveal:
//!
//! 1. Prover commits to a fresh masking message `g` of the same length as
//!    the final folded message `f`.
//! 2. Prover sends `mask_sum = ⟨g, coefficients⟩`.
//! 3. Verifier sends `γ`.
//! 4. Prover sends `f* = g + γ · f` (this masks `f` because `g` is fresh).
//! 5. Both parties spot-check codeword consistency at random positions:
//!    `enc(f*)[pos] = enc(g)[pos] + γ · enc(f)[pos]`.
//! 6. Target check: `⟨f*, coefficients⟩ = mask_sum + γ · constraint.value`.

use spongefish::{ProverState, VerificationError, VerificationResult, VerifierState};

use super::code::{AdditiveCode, ReedSolomon};
use super::commitment::{CodeCommitmentHandle, CodeCommitmentProverHandle};
use super::linear_form::{FoldedFormHandle, LinearConstraint, LinearFormHandle};
use super::transcript_io::{
	read_opening, sample_positions_prover, sample_positions_verifier, write_opening,
};
use super::vc::{MerkleVc, VectorCommitment};
use crate::field::Goldilocks4;

/// Number of inner codeword symbols (and Merkle leaves) per outer query.
/// `g` is committed without interleaving, so each leaf is a single
/// Goldilocks4 wrapped as a `Vec<_>` of length 1 to fit the [`MerkleVc`]
/// alphabet.
const G_LEAF_WIDTH: usize = 1;

/// HVZK base case parameters: number of codeword spot-checks.
pub(crate) struct TrivialZk {
	pub(crate) queries: usize,
}

/// Prover-side intermediate state for the HVZK base case, held between
/// [`TrivialZk::commit_mask`] and [`TrivialZk::finish`].
pub(crate) struct TrivialZkProverState {
	g_msg: Vec<Goldilocks4>,
	g_vc: MerkleVc,
	g_vc_state: <MerkleVc as VectorCommitment>::CommitState,
}

impl TrivialZk {
	/// Commit to a fresh mask `g` matching the final folded message length,
	/// then write `mask_sum = ⟨g, coefficients⟩` to the transcript.
	pub(crate) fn commit_mask(
		&self,
		transcript: &mut ProverState,
		input_msg_len: usize,
		coefficients: &[Goldilocks4],
		mask_seed: &[u8; 32],
	) -> TrivialZkProverState {
		assert_eq!(coefficients.len(), input_msg_len);

		let g_msg = derive_g_msg(mask_seed, input_msg_len);

		let g_code = ReedSolomon::<Goldilocks4>::new(input_msg_len);
		let g_codeword = g_code.encode(&g_msg);
		let g_leaves: Vec<Vec<Goldilocks4>> = g_codeword.iter().map(|&x| vec![x]).collect();
		let g_vc = MerkleVc::new(g_codeword.len());
		let (g_root, g_vc_state) = g_vc.commit(&g_leaves);
		transcript.prover_message(&g_root);

		let mask_sum: Goldilocks4 = g_msg.iter().zip(coefficients).map(|(g, c)| *g * *c).sum();
		transcript.prover_message(&mask_sum);

		// `g_codeword` is consumed into `g_vc` via leaf hashes; no separate
		// caching is needed after commit.
		drop(g_codeword);
		TrivialZkProverState {
			g_msg,
			g_vc,
			g_vc_state,
		}
	}

	/// Steps 3–9 of Construction 7.2: read `γ`, send `f* = g + γ · f`, open
	/// both codewords at random positions, write openings to the transcript.
	pub(crate) fn finish<CCH>(
		&self,
		transcript: &mut ProverState,
		input: CCH,
		mask_state: TrivialZkProverState,
	) where
		CCH: CodeCommitmentProverHandle<
				VC = MerkleVc,
				Code: AdditiveCode<InputAlphabet = Goldilocks4>,
			>,
	{
		let gamma: Goldilocks4 = transcript.verifier_message();
		let f_msg = input.msg();
		assert_eq!(f_msg.len(), mask_state.g_msg.len());

		let f_star: Vec<Goldilocks4> = mask_state
			.g_msg
			.iter()
			.zip(f_msg)
			.map(|(g, f)| *g + gamma * *f)
			.collect();
		for el in &f_star {
			transcript.prover_message(el);
		}

		// Positions live in the inner codeword's index space (length
		// `g_codeword.len() == input.codeword_len()` after fold). The input
		// VC has `input.codeword_len()` leaves (outer codeword length); the
		// g VC has the same number of leaves but with 1-symbol alphabet.
		let positions = sample_positions_prover(transcript, self.queries, input.codeword_len());
		let input_opening = input.open(&positions);
		let g_opening = mask_state.g_vc.open(&mask_state.g_vc_state, &positions);

		write_opening(transcript, &input_opening);
		write_opening(transcript, &g_opening);
	}
}

/// Run the HVZK base case verifier inline against the closing folded
/// commitment. Returns `Ok(())` iff all consistency and target checks pass.
pub(crate) fn verify_trivial_zk<EC>(
	transcript: &mut VerifierState,
	queries: usize,
	folded_commitment: super::commitment::FoldedCodeCommitmentHandle<EC, MerkleVc>,
	constraint: LinearConstraint<FoldedFormHandle<Goldilocks4>>,
) -> VerificationResult<()>
where
	EC: super::code::LinearCode<Alphabet = Goldilocks4>,
	ReedSolomon<Goldilocks4>:
		AdditiveCode<InputAlphabet = Goldilocks4, OutputAlphabet = Goldilocks4>,
{
	let input_msg_len = folded_commitment.msg_len();

	let g_root: [u8; 32] = transcript.prover_message()?;
	let mask_sum: Goldilocks4 = transcript.prover_message()?;
	let gamma: Goldilocks4 = transcript.verifier_message();
	let f_star = transcript.prover_messages_vec::<Goldilocks4>(input_msg_len)?;

	// Target check: ⟨f*, coefficients⟩ = mask_sum + γ · constraint.value.
	let coefficients = constraint.linear_form_handle.folded_form(&[]);
	if coefficients.len() != f_star.len() {
		return Err(VerificationError);
	}
	let dot: Goldilocks4 = coefficients.iter().zip(&f_star).map(|(a, b)| *a * *b).sum();
	if dot != mask_sum + gamma * constraint.value {
		return Err(VerificationError);
	}

	let positions =
		sample_positions_verifier(transcript, queries, folded_commitment.codeword_len());

	// Read both openings. The input's openings carry `INTERLEAVING = 4`
	// symbols each (the interleaved outer alphabet); g's openings carry
	// exactly one symbol each (no interleaving).
	let path_len = folded_commitment
		.codeword_len()
		.next_power_of_two()
		.trailing_zeros() as usize;
	const INTERLEAVING: usize = 4;
	let input_opening = read_opening(transcript, queries, INTERLEAVING, queries * path_len)?;
	let g_opening = read_opening(transcript, queries, G_LEAF_WIDTH, queries * path_len)?;

	let input_folded_values = folded_commitment.verify_openings(&positions, &input_opening)?;
	let g_vc = MerkleVc::new(folded_commitment.codeword_len());
	if !g_vc.verify(&g_root, &positions, &g_opening) {
		return Err(VerificationError);
	}
	let g_values: Vec<Goldilocks4> = g_opening.openings.iter().map(|o| o[0]).collect();

	let inner_code = super::code::ReedSolomon::<Goldilocks4>::new(input_msg_len);
	let encoded_f_star = inner_code.encode(&f_star);

	for (i, &pos) in positions.iter().enumerate() {
		let expected = g_values[i] + gamma * input_folded_values[i];
		match encoded_f_star.get(pos) {
			Some(actual) if *actual == expected => {}
			_ => return Err(VerificationError),
		}
	}

	Ok(())
}

/// Derive `g_msg` deterministically from a 32-byte `mask_seed` via SHAKE256.
fn derive_g_msg(mask_seed: &[u8; 32], count: usize) -> Vec<Goldilocks4> {
	use crate::hash::{Family, JV_RZK, hash};
	let mut buffer_size = count * 32 * 2 + 32;
	let mut refill = 0u64;
	loop {
		let extra = refill.to_le_bytes();
		let stream = if refill == 0 {
			hash(Family::Xof, JV_RZK, &[mask_seed], buffer_size)
		} else {
			hash(Family::Xof, JV_RZK, &[mask_seed, &extra], buffer_size)
		};
		let mut out = Vec::with_capacity(count);
		let mut cursor = 0usize;
		while out.len() < count && cursor + 32 <= stream.len() {
			let chunk = &stream[cursor..cursor + 32];
			cursor += 32;
			if let Some(g) = Goldilocks4::from_bytes(chunk) {
				out.push(g);
			}
		}
		if out.len() == count {
			return out;
		}
		buffer_size *= 2;
		refill += 1;
	}
}
