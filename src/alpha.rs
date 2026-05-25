//! The batched linear form `α = Σ_t β_t · u(x_t)` used by the verifier.
//!
//! Implements [`crate::whir::LinearFormHandle`] symbolically so the verifier
//! never materialises the length-`M` α vector. The fold is computed
//! per-`MonomialLift` and then linearly combined with the `β_t` weights.

use crate::field::Goldilocks4;
use crate::lift::MonomialLift;
use crate::whir::LinearFormHandle;

/// Symbolic representation of `α = Σ_t β_t · u(x_t)` as seen by the verifier.
pub(crate) struct BatchedAlpha {
	lifts: Vec<MonomialLift>,
	betas: Vec<Goldilocks4>,
	nu: u32,
}

impl BatchedAlpha {
	/// Build the symbolic batched α from `K` evaluation points `xs` and their
	/// Fiat–Shamir weights `betas`.
	pub(crate) fn new(xs: &[Goldilocks4], betas: Vec<Goldilocks4>, nu: u32) -> Self {
		assert_eq!(xs.len(), betas.len(), "xs / betas length mismatch");
		let lifts = xs.iter().map(|&x| MonomialLift::new(x, nu)).collect();
		Self { lifts, betas, nu }
	}
}

impl LinearFormHandle for BatchedAlpha {
	type Alphabet = Goldilocks4;

	fn form_size(&self) -> usize {
		1usize << self.nu
	}

	fn folded_form(&self, rand: &[Self::Alphabet]) -> Vec<Self::Alphabet> {
		let l = self.nu as usize - rand.len();
		let len = 1usize << l;
		let mut acc = vec![Goldilocks4::ZERO; len];
		for (lift, &beta) in self.lifts.iter().zip(self.betas.iter()) {
			let folded = lift.folded(rand);
			debug_assert_eq!(folded.len(), len);
			for (a, f) in acc.iter_mut().zip(folded) {
				*a += beta * f;
			}
		}
		acc
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

	/// The critical pin: symbolic α folding equals materialised α folding,
	/// across many (ν, K, R) configurations.
	#[test]
	fn symbolic_fold_matches_materialized_fold() {
		for nu in [2u32, 3, 4, 5, 7] {
			for k in [1usize, 2, 4, 8] {
				let xs: Vec<Goldilocks4> = (0..k).map(|i| g(7 + i as u64)).collect();
				let betas: Vec<Goldilocks4> = (0..k).map(|i| g(101 + i as u64)).collect();
				let alpha = BatchedAlpha::new(&xs, betas.clone(), nu);

				let n = 1usize << nu;
				let mut explicit_alpha = vec![Goldilocks4::ZERO; n];
				for (x, &beta) in xs.iter().zip(betas.iter()) {
					let lift = MonomialLift::new(*x, nu);
					let u = lift.materialize();
					for (a, &uk) in explicit_alpha.iter_mut().zip(u.iter()) {
						*a += beta * uk;
					}
				}

				for r in 0..=nu {
					let rand: Vec<Goldilocks4> = (0..r).map(|i| g(2000 + i as u64)).collect();
					let symbolic = alpha.folded_form(&rand);
					let explicit = manual_fold(explicit_alpha.clone(), &rand);
					assert_eq!(symbolic, explicit, "ν={nu} K={k} R={r}");
				}
			}
		}
	}
}
