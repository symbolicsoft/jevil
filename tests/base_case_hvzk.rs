//! Lemma 7.3 of ePrint 2026/391: the base-case prover transcript is identically
//! distributed to the simulator's output for a fixed `(witness, public claim,
//! challenge seed)`, with statistical distance `ζ_C + n · ζ_{C_zk} = 0` for
//! the Reed–Solomon instantiation.
//!
//! Strong simulator-equivalence is asserted in `src/whir/simulator.rs` via the
//! multi-opening flagship test `tests/multi_opening_hvzk.rs` once the simulator
//! is fully wired (plan Task 14 / 17). Integration-level coverage here asserts
//! the base case round-trips through the public API at the smallest deployable
//! `n_star`.

use jevil::{Params, keygen, sign, verify};
use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;

#[test]
fn base_case_round_trip_at_minimal_n_star() {
	let params = Params::new(1);
	let mut rng = ChaCha20Rng::seed_from_u64(0);
	let (pk, sk, cache) = keygen(&mut rng, params);
	let sig = sign(&sk, &pk, &cache, params, b"base-case-test");
	verify(&pk, params, b"base-case-test", &sig).expect("must verify");
}

#[test]
fn base_case_round_trip_at_n_star_7() {
	// n_star = 7 → at least one full code-switching round exercises the
	// base case after a non-trivial WHIR recursion.
	let params = Params::new(7);
	let mut rng = ChaCha20Rng::seed_from_u64(0);
	let (pk, sk, cache) = keygen(&mut rng, params);
	let sig = sign(&sk, &pk, &cache, params, b"base-case-after-recursion");
	verify(&pk, params, b"base-case-after-recursion", &sig).expect("must verify");
}
