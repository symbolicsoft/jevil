//! Jevil parameters and their derivations.
//!
//! Jevil takes a single deployment-time input: the signing budget
//! [`Params::n_star`]. From it, every other parameter is derived by the rules
//! laid out in the paper. The positions-per-signature parameter `K` is a
//! **global constant** (16), not a per-deployment choice — see the
//! [`Params::K`] doc comment for the rationale.

/// Configuration for a Jevil signer/verifier.
///
/// The only field is the signing budget [`Params::n_star`]. All derived
/// quantities (`M`, `N`, `T`, `ν`, `ν'`, `n_cliff`) are computed on demand by
/// methods on `Params` from `n_star` and the global constant `K`.
///
/// **`n_star + 1` must be a power of two** — that is, `n_star ∈ {1, 3, 7, 15,
/// 31, 63, 127, 255, 511, 1023, …}`. This is the recommended regime of the
/// paper, where `(n_star + 1) · K` is itself a power of two, `M = (n_star +
/// 1) · K` exactly, and the cliff fires at signature `n_star + 1`. Outside
/// this regime `M` would round up to the next power of two, leaving a gap
/// between `n_star` and `n_cliff` that erodes the HORS coverage margin
/// (Theorem~7, §5.3 of the paper); [`Params::new`] panics rather than letting
///    a caller deploy into the bad regime by accident.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub struct Params {
	/// The signing budget: at most `n_star` honest signatures should be issued
	/// per `(PublicKey, SecretKey)` pair before the cliff is reached. Must
	/// satisfy `n_star + 1` is a power of two; see the type-level docs.
	pub n_star: u32,
}

impl Params {
	/// Positions revealed per signature. Fixed at **16** for every Jevil
	/// deployment.
	///
	/// `K = 16` is the sweet spot of the HORS / signature-size trade-off:
	/// smaller `K` forces the position-sampling space `T` infeasibly large,
	/// while larger `K` balloons signatures. Making `K` a global constant —
	/// rather than a deployment-time choice — removes a footgun: signers and
	/// verifiers with mismatched `K` would otherwise produce silently
	/// incompatible artifacts.
	pub const K: u32 = 16;

	/// Convenience constructor. Equivalent to `Params { n_star }` with a
	/// recommended-regime check.
	///
	/// Panics if `n_star == 0` or if `n_star + 1` is not a power of two. The
	/// latter restricts callers to the recommended regime `n_star ∈ {1, 3, 7,
	/// 15, 31, 63, 127, 255, 511, 1023, …}` so that `M = (n_star + 1) · K`
	/// exactly and the cliff fires at signature `n_star + 1`. Outside this
	/// regime `M` rounds up to the next power of two, leaving a gap between
	/// `n_star` and `n_cliff` that erodes the HORS coverage margin (§5.3 of
	/// the paper). A `const`-time assertion catches the mistake at
	/// construction rather than producing a silently-misconfigured deployment.
	pub const fn new(n_star: u32) -> Self {
		assert!(n_star >= 1, "Params::new: n_star must be ≥ 1 (per spec §3)");
		// (n_star + 1) must be a power of two ⇔ (n_star + 1) & n_star == 0.
		// `n_star ≥ 1` means `n_star + 1 ≥ 2` so this is never the special-cased
		// `0 & ?`. We pre-check the high bit to keep the addition in u32 range:
		// `u32::MAX` is rejected here because `u32::MAX + 1` overflows, and we
		// would otherwise miss the rejection.
		assert!(
			n_star < u32::MAX,
			"Params::new: n_star + 1 must be a power of two (n_star = u32::MAX overflows)"
		);
		assert!(
			((n_star + 1) & n_star) == 0,
			"Params::new: n_star + 1 must be a power of two (recommended regime: \
			 n_star ∈ {{1, 3, 7, 15, 31, 63, 127, 255, 511, 1023, …}})"
		);
		Self { n_star }
	}

	/// `ν = log₂((n* + 1) · K)` — the **cliff dimension** exponent. Since
	/// `Params::new` enforces that `n* + 1` is a power of two and `K = 16 =
	/// 2⁴`, `(n* + 1) · K` is exactly a power of two and the `⌈·⌉` of the
	/// paper's formula is a no-op.
	pub fn nu(&self) -> u32 {
		let prod = (self.n_star as u64 + 1) * Self::K as u64;
		ceil_log2_u64(prod)
	}

	/// `M = 2^ν`. The number of honest polynomial coefficients in `c^pad`.
	/// The secret polynomial `f` has degree `D = M − 1`.
	pub fn m(&self) -> usize {
		1usize << self.nu()
	}

	/// `D = M − 1`. The degree bound on the secret polynomial `f`.
	pub fn d(&self) -> usize {
		self.m() - 1
	}

	/// `ν' = ⌈log₂(M + 384)⌉` — the **commit dimension** exponent.
	///
	/// The extra 384 entries hold uniform random masks that absorb WHIR's
	/// codeword leakage. For `M ≥ 384` this gives `ν' = ν + 1`.
	pub fn nu_prime(&self) -> u32 {
		ceil_log2_u64(self.m() as u64 + 384)
	}

	/// `N = 2^ν'`. The length of the padded coefficient vector that is
	/// committed via WHIR.
	pub fn n(&self) -> usize {
		1usize << self.nu_prime()
	}

	/// `T = nextpow2(n* · K · 2^((λ + λ_H)/K))` — the size of the
	/// position-sampling space.
	///
	/// With `λ = λ_H = 128` (target classical security level and adversary
	/// random-oracle query budget) and `K = 16`, this reduces to
	/// `nextpow2(n* · 2²⁰) = n* · 1_048_576` rounded up to the next power of
	/// two. The `λ_H/K` term in the exponent absorbs the adversary's grinding
	/// multiplier `q_H ≤ 2^{λ_H}` on the per-attempt HORS coverage bound
	/// `(n* K / T)^K`, so the post-grinding bound stays below `2^{-λ}` at
	/// `n* = T / (K · 2^{(λ + λ_H)/K})`.
	pub fn t(&self) -> usize {
		// K = 16, λ + λ_H = 256 ⇒ 2^((λ + λ_H)/K) = 2^16 = 65_536 exactly.
		let base = (self.n_star as u64) * (Self::K as u64) * 65_536;
		(base.max(1) as usize).next_power_of_two()
	}

	/// `n_cliff = M / K` — the first signature index at which an outsider
	/// holds ≥ `M` distinct `(x, f(x))` pairs (worst case) and can therefore
	/// recover `f` by Lagrange interpolation.
	///
	/// Because `Params::new` restricts to the recommended regime,
	/// `n_cliff = n_star + 1` exactly.
	pub fn n_cliff(&self) -> usize {
		self.m().div_ceil(Self::K as usize)
	}

	/// 4-byte canonical encoding of `n_star` used in Fiat–Shamir transcripts.
	/// `K` is implicit since it is a global constant.
	pub(crate) fn canonical_bytes(&self) -> [u8; 4] {
		self.n_star.to_le_bytes()
	}
}

/// `⌈log₂(x)⌉` for positive `x`.
fn ceil_log2_u64(x: u64) -> u32 {
	assert!(x > 0, "ceil_log2_u64: input must be positive");
	if x == 1 {
		0
	} else {
		64 - (x - 1).leading_zeros()
	}
}

#[cfg(test)]
mod tests {
	use super::*;

	#[test]
	fn k_is_sixteen() {
		assert_eq!(Params::K, 16);
	}

	#[test]
	fn reference_n_star_1023() {
		let p = Params::new(1023);
		assert_eq!(p.nu(), 14);
		assert_eq!(p.m(), 1 << 14);
		assert_eq!(p.d(), (1 << 14) - 1);
		assert_eq!(p.nu_prime(), 15);
		assert_eq!(p.n(), 1 << 15);
		assert_eq!(p.n_cliff(), 1024);
	}

	#[test]
	fn n_cliff_equals_n_star_plus_one_in_recommended_regime() {
		for n_star in [1u32, 3, 7, 15, 31, 63, 127, 255, 511, 1023] {
			let p = Params::new(n_star);
			assert_eq!(p.n_cliff(), (n_star + 1) as usize, "n_star={n_star}");
		}
	}

	#[test]
	fn nu_prime_satisfies_zk_invariant_above_min() {
		for n_star in [31u32, 63, 127, 1023] {
			let p = Params::new(n_star);
			assert!(p.m() >= 384);
			assert!(
				p.n() - p.m() >= 384,
				"n_star={n_star}: N - M = {}",
				p.n() - p.m()
			);
		}
	}

	#[test]
	fn t_is_power_of_two() {
		for n_star in [1u32, 3, 15, 63, 1023] {
			let p = Params::new(n_star);
			assert!(p.t().is_power_of_two(), "n*={n_star} t={}", p.t());
		}
	}

	#[test]
	fn t_matches_paper_for_reference_configs() {
		// Paper §6.2 (sample sizes at K = 16, λ = λ_H = 128):
		// n*=127 ⇒ T = 2^27,  n*=1023 ⇒ T = 2^30,  n*=16,383 ⇒ T = 2^34.
		assert_eq!(Params::new(127).t(), 1 << 27);
		assert_eq!(Params::new(1023).t(), 1 << 30);
		assert_eq!(Params::new(16_383).t(), 1 << 34);
	}

	#[test]
	fn canonical_bytes_layout() {
		assert_eq!(Params::new(1023).canonical_bytes(), 0x3ff_u32.to_le_bytes());
	}

	#[test]
	#[should_panic(expected = "n_star must be ≥ 1")]
	fn new_rejects_zero() {
		let _ = Params::new(0);
	}

	#[test]
	fn new_rejects_non_recommended_regime() {
		// Every value where (n_star + 1) is not a power of two must panic.
		// We use std::panic::catch_unwind to sweep a handful of bad values
		// without writing one #[should_panic] test per case.
		for bad in [2u32, 4, 5, 6, 8, 9, 16, 64, 100, 128, 1000, 16_384] {
			let result = std::panic::catch_unwind(|| Params::new(bad));
			assert!(
				result.is_err(),
				"Params::new({bad}) should have panicked (n_star + 1 not a power of two)"
			);
		}
	}

	#[test]
	fn new_accepts_full_recommended_set() {
		// Sweep every n_star ∈ {1, 3, 7, …, 16_383} (the deployable range).
		for k in 1u32..=14 {
			let n_star = (1u32 << k) - 1;
			let p = Params::new(n_star);
			assert_eq!(
				p.n_cliff(),
				(n_star + 1) as usize,
				"n_cliff must equal n_star + 1 for every recommended n_star"
			);
		}
	}

	#[test]
	fn ceil_log2_smoke() {
		assert_eq!(ceil_log2_u64(1), 0);
		assert_eq!(ceil_log2_u64(2), 1);
		assert_eq!(ceil_log2_u64(3), 2);
		assert_eq!(ceil_log2_u64(4), 2);
		assert_eq!(ceil_log2_u64(1 << 14), 14);
		assert_eq!(ceil_log2_u64((1 << 14) + 1), 15);
	}
}
