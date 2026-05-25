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

use spongefish::{ProverState, VerificationResult, VerifierState};

use super::base_case::{BaseCase, verify_base_case};
use super::code::{InterleavedCode, ReedSolomon};
use super::codeswitch::CodeswitchHandle;
use super::commitment::{
	CodeCommitment, CodeCommitmentHandle, CodeCommitmentProverHandle, ExplicitCodeCommitmentHandle,
	FoldedCodeCommitmentHandle, FoldedCodeCommitmentProverState,
};
use super::evader::{ETA, OodEvader, OodEvaderHandle};
use super::linear_form::{
	FoldedFormHandle, LinearCombinationForm, LinearConstraint, LinearForm, LinearFormHandle,
};
use super::mask_stack::MaskStack;
use super::transcript_io::{
	read_opening, sample_positions_prover, sample_positions_verifier, write_opening,
};
use super::vc::{MerkleVc, VectorCommitment};
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
	/// Drive one codeswitch round + embedded HVZK sumcheck.
	///
	/// The `mask_seed`, `salt`, and `mask_stack` let `prove_sumcheck_zk`
	/// commit fresh sumcheck masks and push them onto the stack with
	/// their `(α, sl, target)`.
	fn prove(
		&self,
		transcript: &mut ProverState,
		input: ProverFolded,
		mut constraint: LinearForm<Goldilocks4>,
		mask_seed: &[u8; 32],
		salt: &[u8],
		mask_stack: &mut super::mask_stack::MaskStack,
	) -> (LinearForm<Goldilocks4>, ProverFolded) {
		let output = self
			.output_commitment
			.commit(transcript, input.msg().to_vec());

		// Construction 9.7 (privacy padding for OOD). Commit a fresh
		// C_zk-encoded padding mask whose first ETA entries (= r') will be
		// added to the ETA OOD answers to make them statistically
		// independent of the witness. The verifier reads the root,
		// reconstructs the OOD answers' sl_o contribution at the
		// joint-target check, and verifies codeword consistency at the
		// base case.
		let l_zk_inner = crate::params::Params::M_ZK - crate::params::Params::T_ZK;
		let t_zk = crate::params::Params::T_ZK;
		let zk_enc = super::encoding::ZkEncoding::new(l_zk_inner, t_zk);
		let padding_msg = super::base_case::derive_field_vec(
			mask_seed,
			&codeswitch_pad_msg_salt(salt),
			l_zk_inner,
		);
		let padding_r =
			super::base_case::derive_field_vec(mask_seed, &codeswitch_pad_r_salt(salt), t_zk);
		let padding_codeword = zk_enc.encode_with(&padding_msg, &padding_r);
		let padding_leaves: Vec<Vec<Goldilocks4>> =
			padding_codeword.iter().map(|&x| vec![x]).collect();
		let padding_vc = super::vc::MerkleVc::new(padding_codeword.len());
		let (padding_root, padding_vc_state) = padding_vc.commit(&padding_leaves);
		transcript.prover_message(&padding_root);

		// Spec §2.3 fixes η = 2: draw two independent OOD seeds, then return
		// two evaluations and two corresponding linear-form constraints.
		let ood_seeds: Vec<Goldilocks4> = (0..ETA).map(|_| transcript.verifier_message()).collect();
		// Privacy-padded OOD: y_i = ze(ρ_i)·folded_msg + padding_msg[i] for
		// i ∈ [ETA]. The first ETA entries of padding_msg are the fresh r'.
		let mut ood_answers = self.ood_ze.apply(&output.msg, &ood_seeds);
		for (i, y) in ood_answers.iter_mut().enumerate() {
			*y += padding_msg[i];
		}
		transcript.prover_message(&ood_answers);

		let positions = sample_positions_prover(transcript, self.queries, input.codeword_len());
		let openings = input.open(&positions);

		write_opening(transcript, &openings);

		let batch_rand = transcript.verifier_messages_vec(ood_answers.len() + self.queries);

		// Fold the OOD constraints (the ze part only) into the accumulator.
		// The privacy-padding part (r') is absorbed by the padding mask
		// oracle pushed below; this keeps main_constraint k-dimensional.
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

		// Build the padding mask's sl_o and target. sl_o has batch_rand[0..ETA]
		// in its first ETA slots and zeros elsewhere — this picks out the
		// r' = padding_msg[..ETA] entries when evaluated against the mask.
		let mut padding_sl_o = vec![Goldilocks4::ZERO; l_zk_inner];
		padding_sl_o[..ETA].copy_from_slice(&batch_rand[..ETA]);
		let padding_target: Goldilocks4 = padding_sl_o
			.iter()
			.zip(&padding_msg)
			.map(|(a, b)| *a * *b)
			.sum();
		transcript.prover_message(&padding_target);

		// Push the padding mask onto the stack with (α=1, sl_o, target=μ).
		let padding_handle = super::mask_stack::MaskOracleHandle::new_prover(
			padding_msg,
			padding_r,
			padding_vc,
			padding_vc_state,
		);
		mask_stack.push_padding_mask(padding_handle, Goldilocks4::ONE, padding_sl_o);
		// `push_padding_mask` leaves target=ZERO; patch it.
		if let Some(last_mc) = mask_stack.constraints.last_mut() {
			last_mc.target = padding_target;
		}

		// HVZK sumcheck. Carry-in handling: the prover's running sum for the
		// sumcheck is ⟨a, b⟩, which equals the verifier's claim MINUS the
		// carry-in mask contribution mask_carry_in_pre. The verifier
		// subtracts mask_carry_in_pre from its claim before its own HVZK
		// init step; both sides then run prove_sumcheck / verify_sumcheck
		// in sync. After the sumcheck, the verifier adds back ε · carry-in
		// to restore the joint claim. The prover doesn't track claim
		// explicitly; downstream IORs read coefficients/mask_stack.
		let (folded_state, folded_constraint, mask_handles, gammas, epsilon, mask_targets) =
			super::sumcheck::prove_sumcheck_zk(transcript, output, constraint, mask_seed, salt);
		mask_stack.scale_alphas(epsilon);
		mask_stack.push_sumcheck_masks(mask_handles, gammas, mask_targets);
		let folded_constraint = folded_constraint * epsilon;
		(folded_constraint, folded_state)
	}
}

/// Salt builders for the per-codeswitch padding mask derivations.
fn codeswitch_pad_msg_salt(salt: &[u8]) -> Vec<u8> {
	let mut s = salt.to_vec();
	s.extend_from_slice(b"::cs::pad_msg");
	s
}

fn codeswitch_pad_r_salt(salt: &[u8]) -> Vec<u8> {
	let mut s = salt.to_vec();
	s.extend_from_slice(b"::cs::pad_r");
	s
}

/// Runtime-sized WHIR prover specialised to (Goldilocks4 + RS + MerkleVc +
/// OOD).
pub(crate) struct ConcreteWhirProtocol {
	initial_commitment: ProverCommit,
	rounds: Vec<ProverCodeswitch>,
	final_base_case: BaseCase,
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
			final_base_case: BaseCase::new(queries, crate::params::Params::T_ZK),
		}
	}

	/// Prove `⟨msg, α⟩ = v` against the Fiat–Shamir prover transcript.
	///
	/// `mask_seed` is fresh per-signature randomness used to derive the
	/// HVZK masking polynomials (Constructions 6.3 and 7.2 of
	/// eprint 2026/391). Must be unique per signature for HVZK.
	pub(crate) fn prove_to_transcript(
		&self,
		transcript: &mut ProverState,
		msg: Vec<Goldilocks4>,
		initial_constraint: LinearForm<Goldilocks4>,
		mask_seed: &[u8; 32],
	) {
		// Each sumcheck commits k mask oracles, runs the HVZK round-poly
		// math, and pushes k masks onto a per-signature mask_stack. Each
		// codeswitch commits its Construction 9.7 padding mask and pushes
		// that onto the stack too. The base case consumes the full stack
		// via the joint Construction 7.2 target check.
		let initial_state = self.initial_commitment.commit(transcript, msg);
		let mut mask_stack = MaskStack::new();

		let (mut folded_state, mut constraint, mask_handles, gammas, epsilon, mask_targets) =
			super::sumcheck::prove_sumcheck_zk(
				transcript,
				initial_state,
				initial_constraint,
				mask_seed,
				b"sumcheck::initial",
			);
		// scale_alphas is a no-op on the empty stack (no carry-in masks).
		mask_stack.scale_alphas(epsilon);
		mask_stack.push_sumcheck_masks(mask_handles, gammas, mask_targets);
		constraint = constraint * epsilon;

		for (round_idx, round) in self.rounds.iter().enumerate() {
			let salt = format!("sumcheck::post-codeswitch-{round_idx}");
			let (next_constraint, next_folded) = round.prove(
				transcript,
				folded_state,
				constraint,
				mask_seed,
				salt.as_bytes(),
				&mut mask_stack,
			);
			folded_state = next_folded;
			constraint = next_constraint;
		}

		// HVZK base case (Construction 7.2). The mask stack contains every
		// sumcheck mask pushed across all sumchecks; the BaseCase consumes
		// them via the joint target check + per-mask local check + codeword
		// consistency.
		let coefficients: Vec<Goldilocks4> = constraint.coefficients().to_vec();
		self.final_base_case.prove(
			transcript,
			folded_state,
			&coefficients,
			&mask_stack,
			mask_seed,
		);
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
	/// Verify one codeswitch round + embedded HVZK sumcheck. Subtracts the
	/// carry-in mask contribution from constraint.value before calling
	/// verify_sumcheck (so the sumcheck's HVZK init scales just the main
	/// part by ε), then adds back ε · carry_in to restore the joint claim.
	fn verify(
		&self,
		transcript: &mut VerifierState,
		input: &VerifierFolded,
		constraint: LinearConstraint<FoldedFormHandle<Goldilocks4>>,
		mask_stack: &mut super::mask_stack::MaskStack,
	) -> VerificationResult<(
		LinearConstraint<FoldedFormHandle<Goldilocks4>>,
		VerifierFolded,
	)> {
		let output = ExplicitCodeCommitmentHandle {
			code: self.output_commitment.code.clone(),
			vc: self.output_commitment.vc.clone(),
			commitment: transcript.prover_message()?,
		};

		// Construction 9.7 padding mask root.
		let padding_root: [u8; 32] = transcript.prover_message()?;

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

		// Construction 9.7 padding mask: read its sent target and push onto
		// the mask stack with sl_o = (batch_rand[0..ETA], 0, ..., 0). This
		// captures the OOD answers' privacy-padding contribution into the
		// joint mask formulation so subsequent base case checks balance.
		let padding_target: Goldilocks4 = transcript.prover_message()?;
		let l_zk_inner = crate::params::Params::M_ZK - crate::params::Params::T_ZK;
		let mut padding_sl_o = vec![Goldilocks4::ZERO; l_zk_inner];
		padding_sl_o[..ETA].copy_from_slice(&batch_rand[..ETA]);
		let padding_mask_handle =
			super::mask_stack::MaskOracleHandle::verifier_root_only(padding_root);
		mask_stack.push_padding_mask(
			padding_mask_handle,
			<Goldilocks4 as effsc::field::SumcheckField>::ONE,
			padding_sl_o,
		);
		if let Some(last_mc) = mask_stack.constraints.last_mut() {
			last_mc.target = padding_target;
		}

		// Wrap in FoldedFormHandle so the type matches what verify_sumcheck_zk
		// expects (its scale field defaults to ONE here — no additional
		// scaling beyond the LinearCombinationForm batching).
		let mask_carry_in_pre = mask_stack.joint_mask_value();
		let batched_constraint = LinearConstraint {
			linear_form_handle: FoldedFormHandle {
				linear_form_handle: Box::new(LinearCombinationForm {
					linear_form_handles: forms,
					combination_rand: coeffs,
				}),
				rand: Vec::new(),
				scale: <Goldilocks4 as effsc::field::SumcheckField>::ONE,
			},
			value: value - mask_carry_in_pre,
		};

		// HVZK sumcheck (Construction 6.3).
		let (folded_output, mut folded_constraint, mask_handles, gammas, epsilon, mask_targets) =
			super::sumcheck::verify_sumcheck_zk(transcript, output, batched_constraint)?;
		// Restore joint claim: add back ε · carry_in_pre. The mask stack
		// alphas then scale by ε so Σ α_i · μ_i (post-scale) = ε · carry_in_pre.
		folded_constraint.value += epsilon * mask_carry_in_pre;
		mask_stack.scale_alphas(epsilon);
		mask_stack.push_sumcheck_masks(mask_handles, gammas, mask_targets);
		Ok((folded_constraint, folded_output))
	}
}

/// Runtime-sized WHIR verifier specialised to (Goldilocks4 + RS + MerkleVc +
/// OOD).
pub(crate) struct ConcreteWhirVerifier {
	initial_commitment: VerifierExplicit,
	rounds: Vec<VerifierCodeswitch>,
	final_base_case: BaseCase,
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
			final_base_case: BaseCase::new(queries, crate::params::Params::T_ZK),
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

		// Wrap the user-provided initial constraint in a FoldedFormHandle (with
		// scale=ONE and empty rand) so verify_sumcheck_zk's input type matches.
		// Initial mask stack is empty → mask_carry_in_pre = 0 → no subtraction.
		let wrapped_initial = LinearConstraint {
			linear_form_handle: FoldedFormHandle {
				linear_form_handle: Box::new(initial_constraint.linear_form_handle),
				rand: Vec::new(),
				scale: <Goldilocks4 as effsc::field::SumcheckField>::ONE,
			},
			value: initial_constraint.value,
		};

		let mut mask_stack = super::mask_stack::MaskStack::new();
		let (mut folded_commitment, mut constraint, mask_handles, gammas, epsilon, mask_targets) =
			super::sumcheck::verify_sumcheck_zk(transcript, initial_commitment, wrapped_initial)?;
		mask_stack.scale_alphas(epsilon); // no-op (empty)
		mask_stack.push_sumcheck_masks(mask_handles, gammas, mask_targets);

		for round in &self.rounds {
			let (next_constraint, next_folded) =
				round.verify(transcript, &folded_commitment, constraint, &mut mask_stack)?;
			constraint = next_constraint;
			folded_commitment = next_folded;
		}

		// HVZK base case (Construction 7.2 with non-empty mask stack).
		verify_base_case(
			transcript,
			self.final_base_case.queries,
			self.final_base_case.mask_queries,
			folded_commitment,
			constraint,
			&mask_stack.oracles,
			&mask_stack.constraints,
		)
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
