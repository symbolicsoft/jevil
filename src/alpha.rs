//! The batched linear form `α = Σ_t β_t · u(x_t)` used by the verifier.
//!
//! Implements [`crate::whir::LinearFormHandle`] symbolically so the verifier
//! never materialises the (WHIR-embedded) length-`N` α vector. The fold is
//! computed per-`MonomialLift` and then linearly combined with the `β_t`
//! weights. `BatchedAlpha` is logically a length-`M` form (just `Σ_t β_t
//! u(x_t)`); the length-`N` form_size it reports is the embedding into the
//! WHIR primitive's wire format (see `lift.rs` module docs).

use crate::field::Goldilocks4;
use crate::lift::MonomialLift;
use crate::whir::LinearFormHandle;

/// Symbolic representation of `α = Σ_t β_t · u(x_t)` as seen by the verifier.
pub(crate) struct BatchedAlpha {
	lifts: Vec<MonomialLift>,
	betas: Vec<Goldilocks4>,
	nu_prime: u32,
}

impl BatchedAlpha {
	/// Build the symbolic batched α from `K` evaluation points `xs` and their
	/// Fiat–Shamir weights `betas`.
	pub(crate) fn new(xs: &[Goldilocks4], betas: Vec<Goldilocks4>, nu: u32, nu_prime: u32) -> Self {
		assert_eq!(xs.len(), betas.len(), "xs / betas length mismatch");
		let lifts = xs
			.iter()
			.map(|&x| MonomialLift::new(x, nu, nu_prime))
			.collect();
		Self {
			lifts,
			betas,
			nu_prime,
		}
	}
}

impl LinearFormHandle for BatchedAlpha {
	type Alphabet = Goldilocks4;

	fn form_size(&self) -> usize {
		1usize << self.nu_prime
	}

	fn folded_form(&self, rand: &[Self::Alphabet]) -> Vec<Self::Alphabet> {
		let l = self.nu_prime as usize - rand.len();
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

	#[test]
	fn symbolic_fold_matches_materialized_fold() {
		for (nu, nu_prime) in [(2u32, 4), (3, 5), (4, 6), (4, 9), (5, 7)] {
			for k in [1usize, 2, 4, 8] {
				let xs: Vec<Goldilocks4> = (0..k).map(|i| g(7 + i as u64)).collect();
				let betas: Vec<Goldilocks4> = (0..k).map(|i| g(101 + i as u64)).collect();
				let alpha = BatchedAlpha::new(&xs, betas.clone(), nu, nu_prime);

				let n = 1usize << nu_prime;
				let mut explicit_alpha = vec![Goldilocks4::ZERO; n];
				for (x, &beta) in xs.iter().zip(betas.iter()) {
					let lift = MonomialLift::new(*x, nu, nu_prime);
					let u = lift.materialize();
					for (a, &uk) in explicit_alpha.iter_mut().zip(u.iter()) {
						*a += beta * uk;
					}
				}

				for r in 0..=nu_prime {
					let rand: Vec<Goldilocks4> = (0..r).map(|i| g(2000 + i as u64)).collect();
					let symbolic = alpha.folded_form(&rand);
					let explicit = manual_fold(explicit_alpha.clone(), &rand);
					assert_eq!(symbolic, explicit, "ν={nu} ν'={nu_prime} K={k} R={r}");
				}
			}
		}
	}

	#[test]
	fn alpha_value_independent_of_whir_pad_region() {
		let nu = 3u32;
		let nu_prime = 5u32;
		let m = 1usize << nu;
		let n = 1usize << nu_prime;
		let xs = vec![g(7), g(11)];
		let betas = vec![g(101), g(103)];
		let alpha = BatchedAlpha::new(&xs, betas, nu, nu_prime);
		let alpha_vec = alpha.folded_form(&[]);
		assert_eq!(alpha_vec.len(), n);

		// Two WHIR-embedded vectors that agree on the first M entries (the
		// coefficient vector) but differ in the trailing WHIR-pad region
		// (the encoding-randomness slots) must produce the same inner
		// product.
		let mut m_a = vec![g(0); n];
		let mut m_b = vec![g(0); n];
		for k in 0..m {
			m_a[k] = g(31 + k as u64);
			m_b[k] = g(31 + k as u64);
		}
		for k in m..n {
			m_a[k] = g(1_000 + k as u64);
			m_b[k] = g(9_999 + k as u64);
		}
		let dot_a: Goldilocks4 = m_a.iter().zip(alpha_vec.iter()).map(|(c, a)| *c * *a).sum();
		let dot_b: Goldilocks4 = m_b.iter().zip(alpha_vec.iter()).map(|(c, a)| *c * *a).sum();
		assert_eq!(dot_a, dot_b);
	}
}
