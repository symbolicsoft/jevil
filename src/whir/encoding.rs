//! Reed–Solomon ZK encoding per Proposition 3.19 of eprint 2026/391.
//!
//! For Reed–Solomon `C := RS[F, L, ℓ + t]` the Prop 3.19 ZK encoding is
//!
//! ```text
//! Enc_C(msg, r) := NTT(msg ‖ r) padded to codeword_len = 4·(ℓ + t).
//! ```
//!
//! Any subset `S ⊆ [codeword_len]` with `|S| ≤ t` is statistically independent
//! of `msg` over uniformly random `r` — perfect ZK (`ζ_C = 0` for RS).
//!
//! The simulator [`ZkEncoding::simulate`] exploits this by returning `|S|`
//! uniformly random field elements, which is exactly the joint distribution of
//! `Enc_C(msg, r)[S]` over fresh `r` whenever `|S| ≤ t`.

use rand::{CryptoRng, RngCore};

use crate::field::Goldilocks4;

/// Reed–Solomon ZK encoding with Prop 3.19 randomness sizing.
pub(crate) struct ZkEncoding {
	/// Length of the honest message (call it `ℓ`).
	pub(crate) msg_len: usize,
	/// Length of the Prop 3.19 encoding randomness (call it `t`). Any subset
	/// of ≤ `t` codeword positions is statistically independent of the honest
	/// message over uniform `r`.
	pub(crate) rand_len: usize,
	/// Codeword length, `4 · (msg_len + rand_len)`. Must be a power of two.
	pub(crate) codeword_len: usize,
}

#[allow(dead_code)] // `encode`/`simulate`/`generator_row` are Stage B / simulator-test surfaces.
impl ZkEncoding {
	/// Construct a ZK encoding. Panics unless `msg_len + rand_len` is a power
	/// of two (required by the underlying NTT).
	pub(crate) fn new(msg_len: usize, rand_len: usize) -> Self {
		let total = msg_len + rand_len;
		assert!(
			total.is_power_of_two(),
			"ZkEncoding::new: msg_len + rand_len must be a power of two, got {total}"
		);
		Self {
			msg_len,
			rand_len,
			codeword_len: total * 4,
		}
	}

	/// The Prop 3.19 query budget `t`: any subset of ≤ `t` codeword positions
	/// is perfectly hidden over uniform `r`.
	pub(crate) fn query_budget(&self) -> usize {
		self.rand_len
	}

	/// Honest encode: `NTT(msg ‖ r)` with `r` freshly sampled from `rng`.
	pub(crate) fn encode<R: RngCore + CryptoRng>(
		&self,
		msg: &[Goldilocks4],
		rng: &mut R,
	) -> Vec<Goldilocks4> {
		let mut r = Vec::with_capacity(self.rand_len);
		while r.len() < self.rand_len {
			let mut bytes = [0u8; 32];
			rng.fill_bytes(&mut bytes);
			if let Some(g) = Goldilocks4::from_bytes(&bytes) {
				r.push(g);
			}
		}
		self.encode_with(msg, &r)
	}

	/// Encode with caller-supplied randomness. Used by [`crate::keygen`] which
	/// derives `r` deterministically from the secret seed via `JV-RZK`.
	pub(crate) fn encode_with(&self, msg: &[Goldilocks4], r: &[Goldilocks4]) -> Vec<Goldilocks4> {
		assert_eq!(msg.len(), self.msg_len);
		assert_eq!(r.len(), self.rand_len);
		let mut padded: Vec<Goldilocks4> = Vec::with_capacity(self.codeword_len);
		padded.extend_from_slice(msg);
		padded.extend_from_slice(r);
		padded.resize(self.codeword_len, Goldilocks4::default());
		Goldilocks4::ntt(padded)
	}

	/// `Sim_C(S)` per Prop 3.19: sample `|S|` uniformly random field elements
	/// from `rng`. For Reed–Solomon, this matches the joint distribution of
	/// `encode_with(msg, r)[S]` over uniform `r` whenever `|S| ≤ query_budget()`.
	/// Panics if `|S| > query_budget()`.
	pub(crate) fn simulate<R: RngCore + CryptoRng>(
		&self,
		query_set: &[usize],
		rng: &mut R,
	) -> Vec<Goldilocks4> {
		assert!(
			query_set.len() <= self.query_budget(),
			"ZkEncoding::simulate: |S|={} exceeds query budget t={}",
			query_set.len(),
			self.query_budget()
		);
		let mut out = Vec::with_capacity(query_set.len());
		while out.len() < query_set.len() {
			let mut bytes = [0u8; 32];
			rng.fill_bytes(&mut bytes);
			if let Some(g) = Goldilocks4::from_bytes(&bytes) {
				out.push(g);
			}
		}
		out
	}

	/// Row `G_C[pos]` of the generator matrix: the length-`(msg_len + rand_len)`
	/// vector such that `encode_with(msg, r)[pos] = ⟨G_C[pos], msg ‖ r⟩`.
	///
	/// For an NTT-based Reed–Solomon code, `G[pos, j] = ω^{pos · j}` where
	/// `ω` is the 2-adic generator of order `codeword_len`. The first
	/// `msg_len + rand_len` columns of the matrix are exactly the columns
	/// that participate in the inner product with `(msg, r)`.
	pub(crate) fn generator_row(&self, pos: usize) -> Vec<Goldilocks4> {
		assert!(pos < self.codeword_len);
		let log_n = self.codeword_len.trailing_zeros() as usize;
		let omega = Goldilocks4::two_adic_generator(log_n);
		let base = omega.pow(pos as u64);
		let total = self.msg_len + self.rand_len;
		let mut row = Vec::with_capacity(total);
		let mut p = Goldilocks4::ONE;
		for _ in 0..total {
			row.push(p);
			p *= base;
		}
		row
	}
}

#[cfg(test)]
mod tests {
	use super::*;
	use crate::field::Goldilocks;
	use rand::SeedableRng;
	use rand_chacha::ChaCha20Rng;

	fn g4(v: u64) -> Goldilocks4 {
		Goldilocks4::new([
			Goldilocks::new(v),
			Goldilocks::new(0),
			Goldilocks::new(0),
			Goldilocks::new(0),
		])
	}

	#[test]
	fn encode_decomposes_via_generator_row() {
		let enc = ZkEncoding::new(8, 8);
		let msg: Vec<Goldilocks4> = (0..8).map(|i| g4(i as u64)).collect();
		let r: Vec<Goldilocks4> = (0..8).map(|i| g4(100 + i as u64)).collect();
		let codeword = enc.encode_with(&msg, &r);
		let combined: Vec<Goldilocks4> = msg.iter().chain(r.iter()).copied().collect();
		for pos in [0usize, 3, 7, 15, 31, 47] {
			let row = enc.generator_row(pos);
			let dot: Goldilocks4 = row.iter().zip(&combined).map(|(a, b)| *a * *b).sum();
			assert_eq!(codeword[pos], dot, "generator_row at pos={pos}");
		}
	}

	#[test]
	fn simulator_returns_query_set_length() {
		// The strong simulator-vs-encoding distributional equivalence is the
		// flagship multi-opening test (`tests/multi_opening_hvzk.rs`). Here we
		// only sanity-check the shape: simulate(S) returns |S| field elements
		// drawn from the RNG.
		let enc = ZkEncoding::new(8, 8); // t = 8.
		let query_set = [0usize, 5, 12, 17, 23, 28]; // |S| = 6 ≤ t = 8.
		let mut rng = ChaCha20Rng::seed_from_u64(7);
		let sim = enc.simulate(&query_set, &mut rng);
		assert_eq!(sim.len(), query_set.len());
	}

	#[test]
	#[should_panic(expected = "exceeds query budget")]
	fn simulator_panics_outside_budget() {
		let enc = ZkEncoding::new(4, 4);
		let too_many: Vec<usize> = (0..5).collect();
		let mut rng = ChaCha20Rng::seed_from_u64(0);
		let _ = enc.simulate(&too_many, &mut rng);
	}

	#[test]
	fn encoding_is_linear() {
		// Enc(m1 + m2, r1 + r2) = Enc(m1, r1) + Enc(m2, r2). A trivial sanity
		// check that the NTT is linear, which is the foundation of the
		// generator-row decomposition.
		let enc = ZkEncoding::new(4, 4);
		let m1: Vec<_> = (0..4).map(|i| g4(i as u64)).collect();
		let m2: Vec<_> = (0..4).map(|i| g4(10 + i as u64)).collect();
		let r1: Vec<_> = (0..4).map(|i| g4(100 + i as u64)).collect();
		let r2: Vec<_> = (0..4).map(|i| g4(200 + i as u64)).collect();
		let c1 = enc.encode_with(&m1, &r1);
		let c2 = enc.encode_with(&m2, &r2);
		let m_sum: Vec<_> = m1.iter().zip(&m2).map(|(a, b)| *a + *b).collect();
		let r_sum: Vec<_> = r1.iter().zip(&r2).map(|(a, b)| *a + *b).collect();
		let c_sum = enc.encode_with(&m_sum, &r_sum);
		for i in 0..enc.codeword_len {
			assert_eq!(c1[i] + c2[i], c_sum[i], "linearity at i={i}");
		}
	}
}
