//! Tests Prop 3.19 of ePrint 2026/391: for Reed–Solomon `C` with parameters
//! `(msg_len, rand_len)`, any query set `S` with `|S| ≤ rand_len` has the
//! property that `Enc_C(msg, r)[S]` is identically distributed (over uniform
//! `r`) to `|S|` uniformly random field elements.
//!
//! Strong simulator-vs-encoding distributional equivalence lives in
//! `src/whir/encoding.rs::tests` (which can reach `pub(crate)` `ZkEncoding`
//! directly). Integration-level coverage here asserts the property survives
//! through the public API: distinct seeds produce distinct commitments, and
//! the full sign/verify round-trip works against ZK-encoded commitments.

use jevil::{Params, keygen, sign, verify};
use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;

#[test]
fn round_trip_with_zk_encoding_at_n_star_1() {
	let params = Params::new(1);
	let mut rng = ChaCha20Rng::seed_from_u64(0);
	let (pk, sk, cache) = keygen(&mut rng, params);
	let sig = sign(&sk, &pk, &cache, params, b"jv-zk-encoding-prop319-test");
	verify(&pk, params, b"jv-zk-encoding-prop319-test", &sig).expect("should verify");
}

#[test]
fn distinct_seeds_produce_distinct_pk_roots_and_signatures() {
	// The encoding randomness `r_zk` should ensure that two distinct seeds
	// produce distinct commitments — Prop 3.19's perfect hiding requires
	// `r_zk` be uniformly sampled from the seed (deterministic, but distinct).
	let params = Params::new(3);
	let mut rng_a = ChaCha20Rng::seed_from_u64(42);
	let mut rng_b = ChaCha20Rng::seed_from_u64(99);
	let (pk_a, sk_a, cache_a) = keygen(&mut rng_a, params);
	let (pk_b, sk_b, cache_b) = keygen(&mut rng_b, params);
	assert_ne!(pk_a.root, pk_b.root);
	assert_ne!(sk_a.to_bytes(), sk_b.to_bytes());

	let sig_a = sign(&sk_a, &pk_a, &cache_a, params, b"x");
	let sig_b = sign(&sk_b, &pk_b, &cache_b, params, b"x");
	assert_ne!(sig_a.whir_proof, sig_b.whir_proof);
}
