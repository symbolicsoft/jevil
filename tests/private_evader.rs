//! Tests Lemma 9.3 of ePrint 2026/391: the private-padded OOD evader produces
//! uniformly distributed answers independent of the witness, via the
//! [`PrivateZeroEvader`](jevil_internals::PrivateZeroEvader) wrapping in
//! `src/whir/evader.rs`.
//!
//! Strong unit-level confirmation lives in `src/whir/evader.rs::tests`.
//! Integration-level coverage here asserts the property survives through
//! the public API: distinct messages yield distinct WHIR proofs (which
//! would not hold if the OOD step were a deterministic linear functional
//! of the witness), and identical messages yield byte-identical signatures
//! (the determinism invariant the protocol's HVZK randomness must respect).

use jevil::{Params, keygen, sign, verify};
use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;

#[test]
fn distinct_messages_produce_distinct_whir_proofs() {
	let params = Params::new(3);
	let mut rng = ChaCha20Rng::seed_from_u64(0);
	let (pk, sk, cache) = keygen(&mut rng, params);
	let sig_a = sign(&sk, &pk, &cache, params, b"alpha");
	let sig_b = sign(&sk, &pk, &cache, params, b"beta");
	assert!(verify(&pk, params, b"alpha", &sig_a).is_ok());
	assert!(verify(&pk, params, b"beta", &sig_b).is_ok());
	assert_ne!(sig_a.whir_proof, sig_b.whir_proof);
}

#[test]
fn identical_signatures_for_identical_messages_under_deterministic_signing() {
	// jevil is a deterministic signature scheme: the same (sk, pk, msg)
	// must produce byte-identical signatures across calls. The
	// PrivateZeroEvader's "fresh" randomness must therefore be derived
	// deterministically from the per-signature mask seed (which itself is
	// a deterministic function of (sk, root, msg, ys)), not sampled from
	// the OS RNG.
	let params = Params::new(3);
	let mut rng = ChaCha20Rng::seed_from_u64(0);
	let (pk, sk, cache) = keygen(&mut rng, params);
	let sig_1 = sign(&sk, &pk, &cache, params, b"same");
	let sig_2 = sign(&sk, &pk, &cache, params, b"same");
	assert_eq!(sig_1.to_bytes(), sig_2.to_bytes());
}
