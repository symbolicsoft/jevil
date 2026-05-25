//! HVZK sumcheck (Construction 6.3, Lemma 6.4 of ePrint 2026/391) — integration
//! tests exercising the masked sumcheck path end-to-end via the public Jevil
//! API. Round-trip + determinism at `n_star = 3` (the smallest `n_star` whose
//! `n_star + 1` is a power of two and which Jevil's [`Params::new`] accepts).
//!
//! The strong simulator-vs-real byte-equality witness lives in
//! `tests/multi_opening_hvzk.rs` alongside the full Theorem 4.5 composed
//! simulator. This file is the cheap smoke test that the sumcheck IOR composes
//! correctly into the overall sign/verify pipeline.

use jevil::{Params, keygen, sign, verify};
use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;

#[test]
fn sumcheck_round_trip_at_n_star_3() {
	let params = Params::new(3);
	let mut rng = ChaCha20Rng::seed_from_u64(0);
	let (pk, sk, cache) = keygen(&mut rng, params);
	let sig = sign(&sk, &pk, &cache, params, b"sumcheck-hvzk");
	verify(&pk, params, b"sumcheck-hvzk", &sig).expect("honest sumcheck must verify");
}

#[test]
fn deterministic_signatures_imply_deterministic_sumcheck() {
	let params = Params::new(3);
	let mut rng = ChaCha20Rng::seed_from_u64(0);
	let (pk, sk, cache) = keygen(&mut rng, params);
	let sig_a = sign(&sk, &pk, &cache, params, b"determinism-1");
	let sig_b = sign(&sk, &pk, &cache, params, b"determinism-1");
	assert_eq!(
		sig_a.to_bytes(),
		sig_b.to_bytes(),
		"signing the same message under the same key must produce identical bytes"
	);
}
