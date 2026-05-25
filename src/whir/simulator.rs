//! HVZK simulator for `WhirPcs` — composed per Theorem 4.5 of ePrint 2026/391
//! and extended to multi-opening per Lemma 11 of the Jevil paper.
//!
//! # Status
//!
//! This module realises the **algebraic-helper layer** of the HVZK simulator:
//!
//! - [`PublicClaim`] / [`SimulatedTranscript`] — the public API surface.
//! - [`derive_joint_seed`] — domain-tagged hash combining `n*` per-opening
//!   challenge seeds into the single seed `Sim_C` is invoked with.
//! - [`joint_query_set_budget`] / [`lemma_11_hypothesis_holds`] — Lemma 11
//!   parameter-sizing checks. Both are guaranteed by `Params::nu_prime` for
//!   every deployable `n_star`; the runtime check is here so that downstream
//!   simulator code can panic loudly on any violation.
//! - [`simulate_joint_codeword`] — the multi-opening counterpart of
//!   Prop 3.19's `Sim_C`: samples `|S|` uniformly random field elements from
//!   a caller-provided CSPRNG, representing the public-key codeword's values
//!   at the joint query set. Callers typically seed the CSPRNG with
//!   [`derive_joint_seed`] applied to the per-opening challenge seeds.
//!
//! Together these are the deterministic building blocks any byte-equal NARG
//! simulator must compose, and they witness the core HVZK property at the
//! algebraic level (per-position simulation distributionally matches the real
//! encoding for `|S| ≤ N − M`).
//!
//! ## What's pending
//!
//! [`simulate_single`] and [`simulate_multi`] are **stubs**. A byte-equal NARG
//! simulator that produces transcripts accepted by a real Fiat–Shamir
//! [`spongefish::VerifierState`] requires either
//!
//! 1. **Programmable spongefish** — the simulator dictates each verifier
//!    challenge to fit messages it has already sampled (i.e. random-oracle
//!    programming, which spongefish does not currently expose), or
//! 2. **A parallel interactive-prover code path** — drop FS for a synthetic
//!    verifier whose challenges come from an external stream, so the
//!    simulator can be byte-compared against an equally interactive real
//!    prover.
//!
//! Both are substantial refactors. The protocol IS honest-verifier
//! zero-knowledge by Lemma 11 of the paper because (i) Constructions 6.3 /
//! 7.2 / 9.7 are active end-to-end in [`super::protocol`] and (ii) the
//! parameter sizing `|S| ≤ n* · Q_max ≤ N − M` is met by [`Params::n`] (see
//! [`lemma_11_hypothesis_holds`]). The end-to-end round-trip tests at
//! `n_star ∈ {1, 7, 31, 127, 1023}` provide empirical confirmation that the
//! composed protocol accepts.

// The simulator's helpers are exercised only by the in-module test suite and
// by the byte-level NARG simulator that will land in a follow-up. Until that
// arrives there are no production callers, so suppress the dead-code lint at
// the module level rather than per-item.
#![allow(dead_code)]

use rand::{CryptoRng, RngCore};

use super::encoding::ZkEncoding;
use crate::field::Goldilocks4;
use crate::params::Params;

/// Public input to the simulator — everything an outside observer sees about
/// one Jevil signature, minus the WHIR proof itself.
pub(crate) struct PublicClaim {
	/// The 32-byte zk-WHIR commitment root (= `PublicKey::root`).
	pub root: [u8; 32],
	/// The message being signed.
	pub msg: Vec<u8>,
	/// The `K` revealed evaluations `(y_t) = (f(x_t))`.
	pub y_values: Vec<Goldilocks4>,
}

/// Simulated transcript: the NARG bytes a real signer would write to the
/// spongefish transcript for the given `(public_claim, challenge_seed)`.
pub(crate) struct SimulatedTranscript {
	pub narg_bytes: Vec<u8>,
}

/// Bound on the cumulative WHIR queries to the main public-key codeword across
/// every signature in the multi-opening regime: `n* · Q_max`.
///
/// This is the upper bound on `|S|` that any honest multi-opening transcript
/// can hit. `Params::nu_prime` is sized so that `N − M ≥ n* · Q_max` holds for
/// every deployable `n_star`, which gives Lemma 11's hypothesis by
/// construction.
pub(crate) fn joint_query_set_budget(params: Params) -> u64 {
	params.n_star as u64 * Params::Q_MAX
}

/// Check Lemma 11's hypothesis: the joint query set must (i) fit within
/// `n* · Q_max` and (ii) fit within `N − M = |r_zk|`. Both hold by
/// construction of `Params::nu_prime` for any `n_star` accepted by
/// [`Params::new`]; this runtime check is provided so any future byte-level
/// simulator can panic loudly on misconfiguration rather than silently
/// producing an unsound transcript.
pub(crate) fn lemma_11_hypothesis_holds(params: Params, joint_query_set_size: usize) -> bool {
	let size = joint_query_set_size as u64;
	let budget = joint_query_set_budget(params);
	let n_minus_m = (params.n() - params.m()) as u64;
	size <= budget && size <= n_minus_m
}

/// Derive a deterministic joint challenge seed from `n*` per-opening seeds.
/// Used by [`simulate_joint_codeword`] to seed the single `Sim_C` invocation
/// per multi-opening run.
pub(crate) fn derive_joint_seed(seeds: &[[u8; 32]]) -> [u8; 32] {
	use crate::hash::{Family, JV_RZK, hash};
	let inputs: Vec<&[u8]> = seeds.iter().map(|s| s.as_slice()).collect();
	let bytes = hash(Family::Xof, JV_RZK, &inputs, 32);
	let mut out = [0u8; 32];
	out.copy_from_slice(&bytes);
	out
}

/// Sample the joint codeword's values at the joint query set — the
/// distribution any HVZK simulator must match for the main code's encoding
/// under Lemma 11.
///
/// Uses [`ZkEncoding::simulate`] (Prop 3.19) to draw `|joint_query_set|`
/// uniformly random field elements from `rng`. For Reed–Solomon
/// (`ζ_C = 0`) the resulting distribution is *identical* to the joint
/// distribution of `Enc_C(c, r_zk)[joint_query_set]` over uniform `r_zk`
/// whenever the Lemma 11 hypothesis holds — i.e. perfect HVZK at the
/// per-position level.
///
/// Callers should seed `rng` with [`derive_joint_seed`] applied to the
/// per-opening challenge seeds to keep `Sim_C` deterministic in the
/// challenge stream.
///
/// Panics if [`lemma_11_hypothesis_holds`] would return `false`.
pub(crate) fn simulate_joint_codeword<R: RngCore + CryptoRng>(
	params: Params,
	joint_query_set: &[usize],
	rng: &mut R,
) -> Vec<Goldilocks4> {
	assert!(
		lemma_11_hypothesis_holds(params, joint_query_set.len()),
		"Lemma 11 hypothesis violated: |S|={} but n_star·Q_max={}, N-M={}",
		joint_query_set.len(),
		joint_query_set_budget(params),
		params.n() - params.m()
	);
	let zk_enc = ZkEncoding::new(params.m(), params.n() - params.m());
	zk_enc.simulate(joint_query_set, rng)
}

/// Single-opening NARG-byte simulator — stub. See the module docstring for
/// what blocks a real implementation. Callers needing the *algebraic*
/// per-position simulator should use [`simulate_joint_codeword`].
pub(crate) fn simulate_single(
	_params: Params,
	_claim: &PublicClaim,
	_challenge_seed: [u8; 32],
) -> SimulatedTranscript {
	SimulatedTranscript {
		narg_bytes: Vec::new(),
	}
}

/// Multi-opening NARG-byte simulator — stub orchestrating per-claim calls
/// to [`simulate_single`]. Once the byte-level simulator lands, this routine
/// also performs the Lemma 11 joint-query dispatch:
///
/// 1. enumerate each opening's query set `S_i` via a verifier-without-oracles
///    pass over the per-opening challenge stream,
/// 2. draw a single `ũ ← Sim_C(S)` from [`simulate_joint_codeword`] for
///    `S = ⋃_i S_i`,
/// 3. drive each per-opening simulator with `ũ|_{S_i}`.
///
/// Panics if `claims.len()` exceeds `params.n_star` (the budget) or if
/// `claims.len() ≠ challenge_seeds.len()`.
pub(crate) fn simulate_multi(
	params: Params,
	claims: &[PublicClaim],
	challenge_seeds: &[[u8; 32]],
) -> Vec<SimulatedTranscript> {
	assert_eq!(
		claims.len(),
		challenge_seeds.len(),
		"simulate_multi: |claims| = {} != |challenge_seeds| = {}",
		claims.len(),
		challenge_seeds.len()
	);
	assert!(
		claims.len() <= params.n_star as usize,
		"simulate_multi: |claims| = {} > n_star = {}",
		claims.len(),
		params.n_star
	);
	claims
		.iter()
		.zip(challenge_seeds)
		.map(|(claim, seed)| simulate_single(params, claim, *seed))
		.collect()
}

#[cfg(test)]
mod tests {
	use super::*;
	use rand::SeedableRng;
	use rand_chacha::ChaCha20Rng;

	fn rng_from_seeds(seeds: &[[u8; 32]]) -> ChaCha20Rng {
		ChaCha20Rng::from_seed(derive_joint_seed(seeds))
	}

	#[test]
	fn joint_seed_is_deterministic() {
		let seeds = [[1u8; 32], [2u8; 32], [3u8; 32]];
		assert_eq!(derive_joint_seed(&seeds), derive_joint_seed(&seeds));
	}

	#[test]
	fn joint_seed_distinguishes_inputs() {
		let seeds_a = [[1u8; 32], [2u8; 32]];
		let seeds_b = [[1u8; 32], [3u8; 32]];
		assert_ne!(derive_joint_seed(&seeds_a), derive_joint_seed(&seeds_b));
	}

	#[test]
	fn lemma_11_holds_at_full_budget() {
		for n_star in [1u32, 3, 7, 15, 31, 63, 127, 255, 511, 1023] {
			let params = Params::new(n_star);
			let budget = joint_query_set_budget(params) as usize;
			assert!(
				lemma_11_hypothesis_holds(params, budget),
				"Lemma 11 must hold at full n*·Q_max budget for n_star={n_star}"
			);
			assert!(
				budget <= params.n() - params.m(),
				"n*·Q_max must fit in N-M for n_star={n_star} (got budget={budget}, N-M={})",
				params.n() - params.m()
			);
		}
	}

	#[test]
	fn lemma_11_rejects_over_n_minus_m() {
		let params = Params::new(3);
		let over = params.n() - params.m() + 1;
		assert!(!lemma_11_hypothesis_holds(params, over));
	}

	#[test]
	fn joint_codeword_simulator_returns_query_set_length() {
		let params = Params::new(7);
		let positions: Vec<usize> = (0..32).collect();
		let mut rng = rng_from_seeds(&vec![[5u8; 32]; params.n_star as usize]);
		let sample = simulate_joint_codeword(params, &positions, &mut rng);
		assert_eq!(sample.len(), positions.len());
	}

	#[test]
	fn joint_codeword_simulator_is_deterministic_in_seeds() {
		let params = Params::new(3);
		let positions: Vec<usize> = (0..16).collect();
		let seeds = vec![[9u8; 32]; params.n_star as usize];
		let s1 = simulate_joint_codeword(params, &positions, &mut rng_from_seeds(&seeds));
		let s2 = simulate_joint_codeword(params, &positions, &mut rng_from_seeds(&seeds));
		assert_eq!(s1, s2);
	}

	#[test]
	fn joint_codeword_simulator_differs_across_seeds() {
		let params = Params::new(3);
		let positions: Vec<usize> = (0..16).collect();
		let seeds_a = vec![[0u8; 32]; params.n_star as usize];
		let seeds_b = vec![[1u8; 32]; params.n_star as usize];
		let s1 = simulate_joint_codeword(params, &positions, &mut rng_from_seeds(&seeds_a));
		let s2 = simulate_joint_codeword(params, &positions, &mut rng_from_seeds(&seeds_b));
		assert_ne!(s1, s2);
	}

	#[test]
	#[should_panic(expected = "Lemma 11 hypothesis violated")]
	fn simulate_joint_codeword_panics_outside_budget() {
		let params = Params::new(3);
		let over = params.n() - params.m() + 1;
		let positions: Vec<usize> = (0..over).collect();
		let mut rng = rng_from_seeds(&vec![[0u8; 32]; params.n_star as usize]);
		let _ = simulate_joint_codeword(params, &positions, &mut rng);
	}

	#[test]
	fn multi_opening_dispatch_respects_budget() {
		let params = Params::new(3);
		let claims: Vec<PublicClaim> = (0..3)
			.map(|i| PublicClaim {
				root: [i as u8; 32],
				msg: vec![i as u8],
				y_values: Vec::new(),
			})
			.collect();
		let seeds: Vec<[u8; 32]> = (0..3).map(|i| [i as u8; 32]).collect();
		let sims = simulate_multi(params, &claims, &seeds);
		assert_eq!(sims.len(), claims.len());
	}

	#[test]
	#[should_panic(expected = "n_star")]
	fn multi_opening_rejects_over_budget() {
		let params = Params::new(1);
		let claims: Vec<PublicClaim> = (0..2)
			.map(|i| PublicClaim {
				root: [i as u8; 32],
				msg: vec![i as u8],
				y_values: Vec::new(),
			})
			.collect();
		let seeds = vec![[0u8; 32]; 2];
		let _ = simulate_multi(params, &claims, &seeds);
	}
}
