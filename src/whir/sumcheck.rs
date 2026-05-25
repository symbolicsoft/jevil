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

use effsc::field::SumcheckField;

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

/// HVZK sumcheck mask polynomial degree-bound `ℓ_zk` from Construction 6.3.
/// Set to 3 so the round polynomial degree is `max(2, ℓ_zk − 1) = 2` — same
/// wire format (`[h(0), h(∞)]` per round) as the non-ZK path; only `+1` field
/// element (`mask_sum`) prepended before the rounds.
pub(crate) const HVZK_MASK_LENGTH: usize = 3;

/// Prove `⟨msg, α⟩ = claimed_value` (the sumcheck reduction step). Returns
/// the folded commitment state and folded linear form.
///
/// `mask_coeffs`: if empty, vanilla sumcheck. If non-empty, must be exactly
/// `num_rounds · HVZK_MASK_LENGTH` field elements; runs the HVZK variant
/// (Construction 6.3 of eprint 2026/391).
pub(crate) fn prove_sumcheck<EC, VC>(
	transcript: &mut ProverState,
	input: CodeCommitmentProverState<InterleavedCode<EC>, VC>,
	constraint: LinearForm<EC::InputAlphabet>,
	mask_coeffs: &[EC::Alphabet],
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

	let masks_present = !mask_coeffs.is_empty();
	if masks_present {
		assert_eq!(mask_coeffs.len(), num_rounds * HVZK_MASK_LENGTH);
	}

	let mut a = input.msg.clone();
	let mut b = constraint.into_coefficients();

	// HVZK setup: send mask_sum, receive ε.
	// `sum_multiple_initial = 2^(num_rounds - 1)`; `eval_01(s_j) = s_j(0) +
	// s_j(1) = s_j[0] + (s_j[0]+s_j[1]+s_j[2]) = 2·s_j[0] + s_j[1] + s_j[2]`.
	let mut mask_sum = EC::Alphabet::ZERO;
	let mut mask_rlc = EC::Alphabet::ONE;
	let mut running_sum = if masks_present {
		// initial_sum = ⟨a, b⟩.
		let s: EC::Alphabet = a.iter().zip(b.iter()).map(|(x, y)| *x * *y).sum();
		Some(s)
	} else {
		None
	};
	if masks_present {
		let sum_multiple_initial = pow2::<EC::Alphabet>(num_rounds.saturating_sub(1));
		let total_eval01: EC::Alphabet = mask_coeffs
			.chunks_exact(HVZK_MASK_LENGTH)
			.map(eval_01)
			.fold(EC::Alphabet::ZERO, |acc, x| acc + x);
		mask_sum = total_eval01 * sum_multiple_initial;
		transcript.prover_message(&mask_sum);
		mask_rlc = transcript.verifier_message();
	}

	let mut prev_challenge: Option<EC::Alphabet> = None;
	let half_inv: EC::Alphabet = {
		// half = 2^{-1}. char(F) ≠ 2 is required; Goldilocks satisfies this.
		(EC::Alphabet::ONE + EC::Alphabet::ONE)
			.inverse()
			.expect("char ≠ 2")
	};

	for round_idx in 0..num_rounds {
		if let Some(w) = prev_challenge {
			fold_in_place(&mut a, w);
			fold_in_place(&mut b, w);
		}

		let (q0, q_inf) = round_poly(&a, &b);

		if !masks_present {
			transcript.prover_message(&q0);
			transcript.prover_message(&q_inf);
			let r: EC::Alphabet = transcript.verifier_message();
			prev_challenge = Some(r);
			continue;
		}

		// HVZK round: build the modified univariate and send (h0, h_inf).
		let mask = &mask_coeffs[round_idx * HVZK_MASK_LENGTH..(round_idx + 1) * HVZK_MASK_LENGTH];
		let sum_multiple = pow2::<EC::Alphabet>(num_rounds.saturating_sub(round_idx + 1));

		// q1 derived from sumcheck running claim:
		// q(0) + q(1) = current_sum ⇒ q1 = current_sum − 2·q0 − q_inf.
		let current_sum = running_sum.expect("running sum tracked");
		let q1 = current_sum - q0.double() - q_inf;

		// univariate = sum_multiple · mask + (mask_sum − sum_multiple · eval_01(mask))/2 · e_0
		//            + mask_rlc · (q0, q1, q_inf).
		let constant_adj = (mask_sum - sum_multiple * eval_01(mask)) * half_inv;
		let h0 = sum_multiple * mask[0] + constant_adj + mask_rlc * q0;
		let h1 = sum_multiple * mask[1] + mask_rlc * q1;
		let h_inf = sum_multiple * mask[2] + mask_rlc * q_inf;

		transcript.prover_message(&h0);
		transcript.prover_message(&h_inf);

		let r: EC::Alphabet = transcript.verifier_message();
		// Update sum_running (sumcheck part) and mask_sum (mask part).
		// new_sum = q(r) = q0 + q1·r + q_inf·r²
		// new_mask_sum = univariate(r) − mask_rlc · new_sum
		let new_sum = q0 + q1 * r + q_inf * r * r;
		let new_univariate_at_r = h0 + h1 * r + h_inf * r * r;
		mask_sum = new_univariate_at_r - mask_rlc * new_sum;
		running_sum = Some(new_sum);
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

/// `eval_01(p) = p(0) + p(1)` for a length-`HVZK_MASK_LENGTH` univariate
/// polynomial `p` given by coefficients `(p[0], p[1], p[2])`.
fn eval_01<F: effsc::field::SumcheckField>(coeffs: &[F]) -> F {
	if coeffs.is_empty() {
		return F::ZERO;
	}
	// p(0) + p(1) = p[0] + (p[0] + p[1] + p[2]) = 2·p[0] + p[1] + p[2]
	let mut sum = F::ZERO;
	for c in coeffs {
		sum += *c;
	}
	coeffs[0] + sum
}

/// Compute `2^n` as a field element by repeated doubling.
fn pow2<F: effsc::field::SumcheckField>(n: usize) -> F {
	let mut acc = F::ONE;
	for _ in 0..n {
		acc = acc.double();
	}
	acc
}

// ---------------------------------------------------------------------------
// Verifier
// ---------------------------------------------------------------------------

fn inline_sumcheck_verify<F: WhirField>(
	transcript: &mut VerifierState,
	claimed_sum: F,
	num_rounds: usize,
	hvzk: bool,
) -> VerificationResult<(F, Vec<F>)> {
	let mut claim = claimed_sum;

	// HVZK setup: read mask_sum, send ε; adjust running claim accordingly.
	if hvzk {
		let mask_sum: F = transcript.prover_message()?;
		let mask_rlc: F = transcript.verifier_message();
		// claim = mask_sum + mask_rlc · claim (the new running claim for the
		// combined polynomial Σ s_j + ε · G).
		claim = mask_sum + mask_rlc * claim;
	}

	let mut challenges = Vec::with_capacity(num_rounds);

	for _ in 0..num_rounds {
		let h0: F = transcript.prover_message()?;
		let h_inf: F = transcript.prover_message()?;

		// Sumcheck consistency: q(0) + q(1) = current_claim ⇒
		// h0 + (h0 + h1 + h_inf) = claim ⇒ h1 = claim − 2·h0 − h_inf.
		let h1 = claim - h0.double() - h_inf;
		let r: F = transcript.verifier_message();
		challenges.push(r);

		// q(r) = h0 + h1·r + h_inf·r²
		claim = h0 + h1 * r + h_inf * r * r;
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
	hvzk: bool,
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
		inline_sumcheck_verify(transcript, constraint.value, num_rounds, hvzk)?;

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
