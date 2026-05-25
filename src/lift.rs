//! The monomial-basis lift `u(x)` and its multilinear-extension folding.
//!
//! The lift is the length-`M` vector
//!
//! ```text
//! u(x) := (1, x, x², …, x^{M-1})  ∈ F^M,
//! ```
//!
//! whose inner product with the coefficient vector `c ∈ F^M` is the
//! polynomial evaluation `⟨c, u(x)⟩ = f(x)`. The MLE of `u(x)` on
//! `{0,1}^ν` (with `M = 2^ν`, LSB-first storage) factors as
//!
//! ```text
//! u(x)(s₁, …, s_ν) = ∏_{j=1}^{ν} F_j(s_j),   F_j(s) = 1 − s + s · x^{2^{j-1}}.
//! ```
//!
//! ## Embedding into the WHIR primitive's `F^N` namespace
//!
//! The underlying WHIR primitive (`crate::whir`) operates on length-`N`
//! coefficient vectors, where `N = 2^ν' ≥ M` is the Prop. 3.19 ZK-encoded
//! message length the primitive uses internally to attach its encoding
//! randomness. The length-`M` lift is embedded into the length-`N`
//! namespace by zero-padding the trailing `N − M` coordinates; this is
//! purely a wire-format step at the WHIR API boundary and is invisible to
//! the lift's user-facing semantics. The MLE on `{0,1}^{ν'}` picks up one
//! extra factor `∏_{j=ν+1}^{ν'} (1 − s_j)` from the zero-padding region.
//!
//! [`MonomialLift`] materialises either the length-`M` lift
//! ([`MonomialLift::materialize`]) or its `ν'`-variable WHIR-embedded
//! folded form ([`MonomialLift::folded`]) in `O(ν')` instead of `O(N)`.
//! The verifier uses the folded form via [`crate::alpha::BatchedAlpha`]
//! to avoid materialising any length-`N` vector.

use crate::field::Goldilocks4;

/// A single monomial-basis lift `u(x)` parameterised by `(ν, ν')` and the
/// evaluation point `x`.
pub(crate) struct MonomialLift {
	/// `x_powers[j] = x^{2^j}` for `j ∈ {0, …, ν − 1}`. Empty when `ν = 0`.
	x_powers: Vec<Goldilocks4>,
	nu: u32,
	nu_prime: u32,
}

impl MonomialLift {
	/// Build a lift for `x` with dimensions `M = 2^ν`, `N = 2^ν'`.
	pub(crate) fn new(x: Goldilocks4, nu: u32, nu_prime: u32) -> Self {
		assert!(nu <= nu_prime, "nu={nu} > nu_prime={nu_prime}");
		let mut x_powers = Vec::with_capacity(nu as usize);
		if nu > 0 {
			x_powers.push(x);
			for _ in 1..nu {
				let last = *x_powers.last().unwrap();
				x_powers.push(last * last);
			}
		}
		Self {
			x_powers,
			nu,
			nu_prime,
		}
	}

	/// Build the WHIR-embedded length-`2^ν'` form of the lift,
	/// `(1, x, x², …, x^{M-1}, 0, 0, …, 0)`. The leading `M = 2^ν` slots
	/// are the lift proper; the trailing `2^ν' − M` slots are the
	/// zero-padding by which the lift embeds into the WHIR primitive's
	/// length-`N` namespace (`lift.rs` module docs).
	///
	/// Used by unit tests as a reference implementation against which the
	/// symbolic [`MonomialLift::folded`] path is checked; production
	/// callers (signer and verifier) never materialise a single lift
	/// directly.
	#[allow(dead_code)]
	pub(crate) fn materialize(&self) -> Vec<Goldilocks4> {
		// Append-doubling: start with [1]; at step `j ∈ {1, …, ν'}`, multiply
		// each new "right half" entry by `a_j = x^{2^{j-1}}` (or zero in the
		// WHIR zero-pad region for j > ν).
		let mut v = vec![Goldilocks4::ONE];
		for j in 1..=self.nu_prime as usize {
			let a_j = if (j as u32) <= self.nu {
				self.x_powers[j - 1]
			} else {
				Goldilocks4::ZERO
			};
			let half = v.len();
			v.reserve(half);
			for k in 0..half {
				let prod = v[k] * a_j;
				v.push(prod);
			}
		}
		v
	}

	/// Partial-fold the lift's MLE with `rand` (MSB-first), returning a
	/// length-`2^(ν' − rand.len())` vector. Equivalent to
	/// `fold_evaluations(materialize(), rand)` but `O(ν')` instead of `O(N)`.
	pub(crate) fn folded(&self, rand: &[Goldilocks4]) -> Vec<Goldilocks4> {
		let r = rand.len() as u32;
		assert!(
			r <= self.nu_prime,
			"rand.len()={r} > nu_prime={}",
			self.nu_prime
		);

		// WHIR's MSB half-split binds rand[i] to spec-variable s_{ν' - i}.
		// The "bound" scalar is therefore ∏_{i<r} F_{ν' − i}(rand[i]).
		let mut scalar = Goldilocks4::ONE;
		for (i, &r_i) in rand.iter().enumerate() {
			let j = self.nu_prime - i as u32;
			let a_j = if j <= self.nu {
				self.x_powers[j as usize - 1]
			} else {
				Goldilocks4::ZERO
			};
			// F_j(s) = 1 + s · (a_j − 1).
			let f = Goldilocks4::ONE + r_i * (a_j - Goldilocks4::ONE);
			scalar *= f;
		}

		// Free tensor product over the remaining variables (s_1, …, s_{ν' − r}).
		let l = self.nu_prime - r;
		let mut v = vec![scalar];
		for j in 1..=l as usize {
			let a_j = if (j as u32) <= self.nu {
				self.x_powers[j - 1]
			} else {
				Goldilocks4::ZERO
			};
			let half = v.len();
			v.reserve(half);
			for k in 0..half {
				let prod = v[k] * a_j;
				v.push(prod);
			}
		}
		v
	}
}

#[cfg(test)]
mod tests {
	use super::*;
	use crate::field::Goldilocks;

	fn g(n: u64) -> Goldilocks4 {
		Goldilocks4::new([
			Goldilocks::new(n),
			Goldilocks::new(0),
			Goldilocks::new(0),
			Goldilocks::new(0),
		])
	}

	fn manual_fold(mut v: Vec<Goldilocks4>, rand: &[Goldilocks4]) -> Vec<Goldilocks4> {
		for &w in rand {
			let half = v.len() / 2;
			let mut new = Vec::with_capacity(half);
			for k in 0..half {
				new.push(v[k] + (v[k + half] - v[k]) * w);
			}
			v = new;
		}
		v
	}

	#[test]
	fn materialize_zero_on_whir_pad_region() {
		let lift = MonomialLift::new(g(7), 3, 5);
		let v = lift.materialize();
		assert_eq!(v.len(), 32);
		for (k, val) in v.iter().enumerate() {
			if k >= 8 {
				assert!(val.is_zero(), "WHIR-pad position {k} should be zero");
			}
		}
	}

	#[test]
	fn materialize_matches_horner_below_m() {
		let x = g(7);
		let lift = MonomialLift::new(x, 4, 6);
		let v = lift.materialize();
		let mut x_pow = Goldilocks4::ONE;
		for (k, v_k) in v.iter().take(16).enumerate() {
			assert_eq!(*v_k, x_pow, "mismatch at k={k}");
			x_pow *= x;
		}
	}

	#[test]
	fn folded_matches_materialize_plus_explicit_fold() {
		for (nu, nu_prime) in [(2u32, 3), (3, 5), (4, 5), (5, 7)] {
			let lift = MonomialLift::new(g(11), nu, nu_prime);
			let materialised = lift.materialize();
			for r in 0..=nu_prime {
				let rand: Vec<Goldilocks4> = (0..r).map(|i| g(13 + i as u64)).collect();
				let symbolic = lift.folded(&rand);
				let explicit = manual_fold(materialised.clone(), &rand);
				assert_eq!(symbolic, explicit, "nu={nu} nu'={nu_prime} R={r}");
			}
		}
	}

	#[test]
	fn folded_full_returns_single_value() {
		let lift = MonomialLift::new(g(5), 3, 4);
		let rand: Vec<Goldilocks4> = (0..4).map(|i| g(20 + i)).collect();
		let v = lift.folded(&rand);
		assert_eq!(v.len(), 1);
	}
}
