//! Runtime-sized WHIR protocol specialised to Jevil's fixed parameter set.
//!
//! ## Fixed instantiation
//!
//! - Field: [`Goldilocks4`].
//! - Inner code: [`ReedSolomon`] (rate `1/4`) wrapped in
//!   [`InterleavedCode`] (factor 4 — i.e. 4 inner codewords per outer symbol).
//! - Vector commitment: [`MerkleVc`] (Poseidon2-Goldilocks Merkle tree).
//! - Zero evader: DEEP-FRI [`OodEvader`].
//! - Fold cap: stop folding when the inner message length reaches the
//!   `threshold` (Jevil uses 64).
//!
//! ## Entry points
//!
//! - [`ConcreteWhirProtocol::build`] / [`ConcreteWhirProtocol::prove_to_transcript`]
//!   — prover.
//! - [`ConcreteWhirVerifier::build`] / [`ConcreteWhirVerifier::verify_from_transcript`]
//!   — verifier.
//!
//! Both halves must be built with the *same* `(n_padded_base, queries,
//! threshold)`; mismatched parameters break interoperability silently. The
//! Jevil keygen/sign/verify pipeline drives these directly with the canonical
//! `(N, 32, 64)` triple.

use std::sync::Arc;

use spongefish::{ProverState, VerificationError, VerificationResult, VerifierState};

use super::code::{InterleavedCode, ReedSolomon};
use super::codeswitch::CodeswitchHandle;
use super::commitment::{
	CodeCommitment, CodeCommitmentHandle, CodeCommitmentProverHandle, ExplicitCodeCommitmentHandle,
	FoldedCodeCommitmentHandle, FoldedCodeCommitmentProverState,
};
use super::linear_form::{
	FoldedFormHandle, LinearCombinationForm, LinearConstraint, LinearForm, LinearFormHandle,
};
use super::sumcheck::{prove_sumcheck, verify_sumcheck};
use super::transcript_io::{
	read_opening, sample_positions_prover, sample_positions_verifier, write_opening,
};
use super::trivial::Trivial;
use super::vc::MerkleVc;
use super::zero_evader::{ETA, OodEvader, OodEvaderHandle};
use crate::field::Goldilocks4;

// ---------------------------------------------------------------------------
// Constants matching the reference WHIR parameter set (paper §2.3)
// ---------------------------------------------------------------------------

/// Interleaving factor for the outer code (groups of `INTERLEAVING` inner
/// codewords share one outer alphabet symbol).
const INTERLEAVING: usize = 4;

/// Reciprocal of the Reed–Solomon rate (`codeword_len = RATE_INV · msg_len`).
const RATE_INV: usize = 4;

// ---------------------------------------------------------------------------
// Concrete type aliases
// ---------------------------------------------------------------------------

type ProverFolded = FoldedCodeCommitmentProverState<ReedSolomon<Goldilocks4>, MerkleVc>;
type ProverCommit = CodeCommitment<InterleavedCode<ReedSolomon<Goldilocks4>>, MerkleVc>;

type VerifierFolded = FoldedCodeCommitmentHandle<ReedSolomon<Goldilocks4>, MerkleVc>;
type VerifierExplicit =
	ExplicitCodeCommitmentHandle<InterleavedCode<ReedSolomon<Goldilocks4>>, MerkleVc>;

// ===========================================================================
// PROVER
// ===========================================================================

/// One codeswitch round on the prover side.
pub(crate) struct ProverCodeswitch {
	ood_ze: OodEvader,
	queries: usize,
	output_commitment: ProverCommit,
}

impl ProverCodeswitch {
	/// Drive one codeswitch round. Writes its commitment, OOD answers, opening
	/// data, sumcheck transcripts inline into the prover transcript stream.
	fn prove(
		&self,
		transcript: &mut ProverState,
		input: ProverFolded,
		mut constraint: LinearForm<Goldilocks4>,
	) -> (LinearForm<Goldilocks4>, ProverFolded) {
		let output = self
			.output_commitment
			.commit(transcript, input.msg().to_vec());

		// Spec §2.3 fixes η = 2: draw two independent OOD seeds, then return
		// two evaluations and two corresponding linear-form constraints.
		let ood_seeds: Vec<Goldilocks4> = (0..ETA).map(|_| transcript.verifier_message()).collect();
		let ood_answers = self.ood_ze.apply(&output.msg, &ood_seeds);
		transcript.prover_message(&ood_answers);

		let positions = sample_positions_prover(transcript, self.queries, input.codeword_len());
		let openings = input.open(&positions);

		write_opening(transcript, &openings);

		let batch_rand = transcript.verifier_messages_vec(ood_answers.len() + self.queries);

		// Fold the OOD constraints into the accumulator.
		for (i, ze_constraint) in self
			.ood_ze
			.expanded_constraint(&ood_seeds)
			.into_iter()
			.enumerate()
		{
			constraint += LinearForm::new(ze_constraint) * batch_rand[i];
		}

		// Build the selector vector (one weight per opened position) and
		// pull its image under the transpose of the encoding map.
		//
		// IMPORTANT: use `+=` so duplicate positions accumulate both weights
		// (a plain `=` would overwrite the earlier weight and silently break
		// soundness, since the verifier's `LinearCombinationForm` adds both
		// contributions).
		use super::codeswitch::TransposeCode;
		let mut selector = vec![Goldilocks4::ZERO; input.codeword_len()];
		for (i, &pos) in positions.iter().enumerate() {
			selector[pos] += batch_rand[ood_answers.len() + i];
		}
		constraint += LinearForm::new(input.code().apply_transpose(&selector));

		let (folded_state, folded_constraint) = prove_sumcheck(transcript, output, constraint);
		(folded_constraint, folded_state)
	}
}

/// Runtime-sized WHIR prover specialised to (Goldilocks4 + RS + MerkleVc +
/// OOD).
pub(crate) struct ConcreteWhirProtocol {
	initial_commitment: ProverCommit,
	rounds: Vec<ProverCodeswitch>,
	final_trivial: Trivial,
}

impl ConcreteWhirProtocol {
	/// Build the protocol for a length-`n_padded_base` input vector.
	///
	/// Folds by `INTERLEAVING = 4` each round until the inner message length
	/// reaches `threshold` (Jevil uses 64). `queries` is the number of
	/// in-domain Merkle queries per round.
	pub(crate) fn build(n_padded_base: usize, queries: usize, threshold: usize) -> Self {
		let initial_inner_msg_len = n_padded_base / INTERLEAVING;
		let initial_commitment = make_commitment(initial_inner_msg_len);

		let mut rounds = Vec::new();
		let mut current_inner_msg_len = initial_inner_msg_len;
		while current_inner_msg_len > threshold {
			let next_inner_msg_len = current_inner_msg_len / INTERLEAVING;
			rounds.push(ProverCodeswitch {
				ood_ze: OodEvader::new(current_inner_msg_len),
				queries,
				output_commitment: make_commitment(next_inner_msg_len),
			});
			current_inner_msg_len = next_inner_msg_len;
		}

		Self {
			initial_commitment,
			rounds,
			final_trivial: Trivial { queries },
		}
	}

	/// Prove `⟨msg, α⟩ = v` against the Fiat–Shamir prover transcript.
	pub(crate) fn prove_to_transcript(
		&self,
		transcript: &mut ProverState,
		msg: Vec<Goldilocks4>,
		initial_constraint: LinearForm<Goldilocks4>,
	) {
		let initial_state = self.initial_commitment.commit(transcript, msg);
		let (mut folded_state, mut constraint) =
			prove_sumcheck(transcript, initial_state, initial_constraint);

		for round in &self.rounds {
			let (next_constraint, next_folded) = round.prove(transcript, folded_state, constraint);
			folded_state = next_folded;
			constraint = next_constraint;
		}

		let final_opening = self.final_trivial.prove(transcript, folded_state);
		write_opening(transcript, &final_opening);
	}
}

// ===========================================================================
// VERIFIER
// ===========================================================================

/// One codeswitch round on the verifier side.
pub(crate) struct VerifierCodeswitch {
	ood_ze: OodEvaderHandle,
	queries: usize,
	output_commitment: VerifierExplicit,
}

impl VerifierCodeswitch {
	/// Verify one codeswitch round. Reads the prover-written messages in
	/// exactly the order [`ProverCodeswitch::prove`] wrote them.
	fn verify(
		&self,
		transcript: &mut VerifierState,
		input: &VerifierFolded,
		constraint: LinearConstraint<FoldedFormHandle<Goldilocks4>>,
	) -> VerificationResult<(
		LinearConstraint<FoldedFormHandle<Goldilocks4>>,
		VerifierFolded,
	)> {
		let output = ExplicitCodeCommitmentHandle {
			code: self.output_commitment.code.clone(),
			vc: self.output_commitment.vc.clone(),
			commitment: transcript.prover_message()?,
		};

		let ood_seeds: Vec<Goldilocks4> = (0..ETA).map(|_| transcript.verifier_message()).collect();
		let ze_constraints = self.ood_ze.zero_evader_handles(&ood_seeds);
		let ood_answers = transcript.prover_messages_vec::<Goldilocks4>(ze_constraints.len())?;
		let positions = sample_positions_verifier(transcript, self.queries, input.codeword_len());

		// In a codeswitch round the queried commitment is the *previous* round's
		// folded interleaved code, so each opening is a length-`INTERLEAVING`
		// vector and the Merkle paths go up `log₂(codeword_len)` levels.
		let path_len_per_opening =
			input.codeword_len().next_power_of_two().trailing_zeros() as usize;
		let openings = read_opening(
			transcript,
			self.queries,
			INTERLEAVING,
			self.queries * path_len_per_opening,
		)?;
		let opened = input.verify_openings(&positions, &openings)?;

		let ood_len = ood_answers.len();
		let batch_rand: Vec<Goldilocks4> = (0..ood_len + self.queries)
			.map(|_| transcript.verifier_message())
			.collect();

		let mut value = constraint.value;
		let mut forms: Vec<Box<dyn LinearFormHandle<Alphabet = Goldilocks4>>> =
			vec![Box::new(constraint.linear_form_handle)];
		let mut coeffs = vec![<Goldilocks4 as effsc::field::SumcheckField>::ONE];

		for (i, (answer, ze_constraint)) in ood_answers.into_iter().zip(ze_constraints).enumerate()
		{
			value += answer * batch_rand[i];
			forms.push(Box::new(ze_constraint));
			coeffs.push(batch_rand[i]);
		}

		for (i, (&pos, opening)) in positions.iter().zip(opened).enumerate() {
			value += opening * batch_rand[ood_len + i];
			forms.push(Box::new(input.code().apply_transpose_handle(pos)));
			coeffs.push(batch_rand[ood_len + i]);
		}

		let batched_constraint = LinearConstraint {
			linear_form_handle: LinearCombinationForm {
				linear_form_handles: forms,
				combination_rand: coeffs,
			},
			value,
		};

		let (folded_output, folded_constraint) =
			verify_sumcheck(transcript, output, batched_constraint)?;
		Ok((folded_constraint, folded_output))
	}
}

/// Runtime-sized WHIR verifier specialised to (Goldilocks4 + RS + MerkleVc +
/// OOD).
pub(crate) struct ConcreteWhirVerifier {
	initial_commitment: VerifierExplicit,
	rounds: Vec<VerifierCodeswitch>,
	final_trivial: Trivial,
}

impl ConcreteWhirVerifier {
	/// Build the verifier. Parameters must match what
	/// [`ConcreteWhirProtocol::build`] was called with.
	pub(crate) fn build(n_padded_base: usize, queries: usize, threshold: usize) -> Self {
		let initial_inner_msg_len = n_padded_base / INTERLEAVING;
		let initial_commitment = make_explicit_handle(initial_inner_msg_len);

		let mut rounds = Vec::new();
		let mut current_inner_msg_len = initial_inner_msg_len;
		while current_inner_msg_len > threshold {
			let next_inner_msg_len = current_inner_msg_len / INTERLEAVING;
			rounds.push(VerifierCodeswitch {
				ood_ze: OodEvaderHandle::new(current_inner_msg_len),
				queries,
				output_commitment: make_explicit_handle(next_inner_msg_len),
			});
			current_inner_msg_len = next_inner_msg_len;
		}

		Self {
			initial_commitment,
			rounds,
			final_trivial: Trivial { queries },
		}
	}

	/// Verify `⟨msg, α⟩ = v` from the FS verifier transcript.
	pub(crate) fn verify_from_transcript<LFH>(
		&self,
		transcript: &mut VerifierState,
		initial_constraint: LinearConstraint<LFH>,
	) -> VerificationResult<()>
	where
		LFH: LinearFormHandle<Alphabet = Goldilocks4> + 'static,
	{
		let initial_commitment = ExplicitCodeCommitmentHandle {
			code: self.initial_commitment.code.clone(),
			vc: self.initial_commitment.vc.clone(),
			commitment: transcript.prover_message()?,
		};

		let (mut folded_commitment, mut constraint) =
			verify_sumcheck(transcript, initial_commitment, initial_constraint)?;

		for round in &self.rounds {
			let (next_constraint, next_folded) =
				round.verify(transcript, &folded_commitment, constraint)?;
			constraint = next_constraint;
			folded_commitment = next_folded;
		}

		// Final trivial step (inlined so the opening lands at the correct
		// transcript offset).
		let msg = transcript.prover_messages_vec::<Goldilocks4>(folded_commitment.msg_len())?;
		let encoded = folded_commitment.encode(&msg);
		let positions = sample_positions_verifier(
			transcript,
			self.final_trivial.queries,
			folded_commitment.codeword_len(),
		);

		let final_path_len_per_opening = folded_commitment
			.codeword_len()
			.next_power_of_two()
			.trailing_zeros() as usize;
		let final_openings = read_opening(
			transcript,
			self.final_trivial.queries,
			INTERLEAVING,
			self.final_trivial.queries * final_path_len_per_opening,
		)?;
		let opened = folded_commitment.verify_openings(&positions, &final_openings)?;

		for (&pos, opening) in positions.iter().zip(&opened) {
			if encoded.get(pos) != Some(opening) {
				return Err(VerificationError);
			}
		}

		let coefficients = constraint.linear_form_handle.folded_form(&[]);
		if coefficients.len() != msg.len() {
			return Err(VerificationError);
		}
		let dot_product: Goldilocks4 = coefficients.into_iter().zip(msg).map(|(a, b)| a * b).sum();
		if dot_product == constraint.value {
			Ok(())
		} else {
			Err(VerificationError)
		}
	}
}

// ---------------------------------------------------------------------------
// Constructor helpers
// ---------------------------------------------------------------------------

fn make_commitment(inner_msg_len: usize) -> ProverCommit {
	let rs = ReedSolomon::<Goldilocks4>::new(inner_msg_len);
	let code = Arc::new(InterleavedCode::new(rs, INTERLEAVING));
	let vc = Arc::new(MerkleVc::new(inner_msg_len * RATE_INV));
	CodeCommitment::new(code, vc)
}

fn make_explicit_handle(inner_msg_len: usize) -> VerifierExplicit {
	let rs = ReedSolomon::<Goldilocks4>::new(inner_msg_len);
	let code = Arc::new(InterleavedCode::new(rs, INTERLEAVING));
	let vc = Arc::new(MerkleVc::new(inner_msg_len * RATE_INV));
	ExplicitCodeCommitmentHandle::new(code, vc, <[u8; 32]>::default())
}
