//! Degree-2 inner-product sumcheck, MSB half-split layout.
//!
//! Both prover and verifier live in this module so the wire format
//! (`[q(0), q(∞)]` per round) is defined exactly once.
//!
//! ## Round polynomial
//!
//! Given vectors `a` and `b`, the sumcheck claim is `Σ_x a(x)·b(x) = c`. The
//! prover sends, each round, the degree-2 polynomial
//!
//! ```text
//! q(X) = h0 + (h1 - h0 - h_inf) · X + h_inf · X²
//! ```
//!
//! transmitted as `[h0, h_inf]`; the verifier derives `h1 = claim − h0` and
//! updates the claim to `q(r)` for its sampled challenge `r`.

use spongefish::{ProverState, VerificationError, VerificationResult, VerifierState};

use super::code::{Field as WhirField, InterleavedCode, LinearCode};
use super::commitment::{
	CodeCommitmentProverState, ExplicitCodeCommitmentHandle, FoldedCodeCommitmentHandle,
};
use super::linear_form::{FoldedFormHandle, LinearConstraint, LinearForm, LinearFormHandle};
use super::vc::VectorCommitment;

// ---------------------------------------------------------------------------
// Shared MSB half-split fold
// ---------------------------------------------------------------------------

/// MSB half-split fold: `a[k] = a[k] + w · (a[k + half] − a[k])` for
/// `k ∈ [0, half)`. The tail entries (when `n` is not a power of two) are
/// treated as paired against an implicit `0`.
fn fold_in_place<F>(a: &mut Vec<F>, w: F)
where
	F: Copy + core::ops::Add<Output = F> + core::ops::Sub<Output = F> + core::ops::Mul<Output = F>,
{
	let n = a.len();
	if n <= 1 {
		return;
	}
	let half = n.next_power_of_two() >> 1;
	for k in 0..half.min(n - half) {
		let lo = a[k];
		let hi = a[k + half];
		a[k] = lo + w * (hi - lo);
	}
	a.truncate(half);
}

/// Compute `(q(0), q(∞))` for the inner-product sumcheck round on `(a, b)`.
fn round_poly<F: effsc::field::SumcheckField>(a: &[F], b: &[F]) -> (F, F) {
	let n = a.len();
	if n <= 1 {
		let v = if n == 1 { a[0] * b[0] } else { F::ZERO };
		return (v, F::ZERO);
	}
	let half = n.next_power_of_two() >> 1;
	let paired = half.min(n - half);

	let mut q0 = F::ZERO;
	let mut q_inf = F::ZERO;

	// Paired region (both lo and hi present).
	for k in 0..paired {
		let al = a[k];
		let ah = a[k + half];
		let bl = b[k];
		let bh = b[k + half];
		q0 += al * bl;
		q_inf += (ah - al) * (bh - bl);
	}

	// Tail region (lo present, hi implicit zero).
	for k in paired..half.min(n) {
		let al = a[k];
		let bl = b[k];
		let dot = al * bl;
		q0 += dot;
		q_inf += dot;
	}

	(q0, q_inf)
}

// ---------------------------------------------------------------------------
// Prover
// ---------------------------------------------------------------------------

/// Prove `⟨msg, α⟩ = claimed_value` (the sumcheck reduction step). Returns
/// the folded commitment state and folded linear form.
pub(crate) fn prove_sumcheck<EC, VC>(
	transcript: &mut ProverState,
	input: CodeCommitmentProverState<InterleavedCode<EC>, VC>,
	constraint: LinearForm<EC::InputAlphabet>,
) -> (
	super::commitment::FoldedCodeCommitmentProverState<EC, VC>,
	LinearForm<EC::InputAlphabet>,
)
where
	EC: LinearCode,
	EC::Alphabet: effsc::field::SumcheckField,
	VC: VectorCommitment<Alphabet = Vec<EC::OutputAlphabet>>,
{
	let interleaving = input.code.interleaving_factor();
	assert!(interleaving > 0 && interleaving.is_power_of_two());
	let num_rounds = interleaving.ilog2() as usize;
	assert_eq!(input.msg.len(), constraint.coefficients().len());

	let mut a = input.msg.clone();
	let mut b = constraint.into_coefficients();
	let mut prev_challenge: Option<EC::Alphabet> = None;

	for _ in 0..num_rounds {
		if let Some(w) = prev_challenge {
			fold_in_place(&mut a, w);
			fold_in_place(&mut b, w);
		}

		let (q0, q_inf) = round_poly(&a, &b);
		// Send element-by-element to match the verifier's per-element reads.
		transcript.prover_message(&q0);
		transcript.prover_message(&q_inf);

		let r: EC::Alphabet = transcript.verifier_message();
		prev_challenge = Some(r);
	}

	if let Some(w) = prev_challenge {
		fold_in_place(&mut a, w);
		fold_in_place(&mut b, w);
	}

	(
		super::commitment::FoldedCodeCommitmentProverState {
			inner: input,
			msg: a,
		},
		LinearForm::new(b),
	)
}

// ---------------------------------------------------------------------------
// Verifier
// ---------------------------------------------------------------------------

fn inline_sumcheck_verify<F: WhirField>(
	transcript: &mut VerifierState,
	claimed_sum: F,
	num_rounds: usize,
) -> VerificationResult<(F, Vec<F>)> {
	let mut claim = claimed_sum;
	let mut challenges = Vec::with_capacity(num_rounds);

	for _ in 0..num_rounds {
		let h0: F = transcript.prover_message()?;
		let h_inf: F = transcript.prover_message()?;

		let h1 = claim - h0;
		let r: F = transcript.verifier_message();
		challenges.push(r);

		// q(r) = h0 + (h1 - h0 - h_inf) · r + h_inf · r²
		let slope = h1 - h0 - h_inf;
		let q_r = h0 + slope * r;
		claim = q_r + h_inf * r * r;
	}

	Ok((claim, challenges))
}

/// Verify the sumcheck reduction. Returns the folded commitment handle and
/// the folded linear-form claim.
#[allow(clippy::type_complexity)] // the return is a generic tuple — aliasing it across both type
// parameters adds more noise than it saves.
pub(crate) fn verify_sumcheck<EC, VC, LFH>(
	transcript: &mut VerifierState,
	commitment: ExplicitCodeCommitmentHandle<InterleavedCode<EC>, VC>,
	constraint: LinearConstraint<LFH>,
) -> VerificationResult<(
	FoldedCodeCommitmentHandle<EC, VC>,
	LinearConstraint<FoldedFormHandle<EC::Alphabet>>,
)>
where
	EC: LinearCode,
	VC: VectorCommitment<Alphabet = Vec<EC::Alphabet>>,
	LFH: LinearFormHandle<Alphabet = EC::Alphabet> + 'static,
{
	use super::commitment::CodeCommitmentHandle;
	let n = commitment.code().interleaving_factor();
	if n == 0
		|| !n.is_power_of_two()
		|| constraint.linear_form_handle.form_size() != commitment.msg_len()
	{
		return Err(VerificationError);
	}

	let num_rounds = n.ilog2() as usize;
	let (final_claim, challenges) =
		inline_sumcheck_verify(transcript, constraint.value, num_rounds)?;

	Ok((
		FoldedCodeCommitmentHandle {
			inner: commitment,
			rand: challenges.clone(),
		},
		LinearConstraint {
			linear_form_handle: FoldedFormHandle {
				linear_form_handle: Box::new(constraint.linear_form_handle),
				rand: challenges,
			},
			value: final_claim,
		},
	))
}
