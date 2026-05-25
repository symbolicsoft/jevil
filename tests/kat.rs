//! Known-answer test: pins the byte layout of `(pk, signature)` for a fixed
//! `(params, seed, msg)` triple. Catches accidental changes to hash framing,
//! Fiat–Shamir transcript layout, Goldilocks limb order, WHIR proof
//! serialisation, etc.
//!
//! When the protocol bytes are intentionally changed, regenerate fixtures
//! with `KAT_UPDATE=1 cargo test --test kat -- --nocapture`.

use jevil::{Params, keygen, sign};
use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;
use sha2::{Digest, Sha256};

#[test]
fn kat_n_star_3() {
	let params = Params::new(3);
	let mut rng = ChaCha20Rng::seed_from_u64(0);
	let (pk, sk, cache) = keygen(&mut rng, params);
	let sig = sign(&sk, &pk, &cache, params, b"jevil-kat-fixture");

	let pk_hex = hex::encode(pk.to_bytes());
	let sig_hash = hex::encode(Sha256::digest(sig.to_bytes()));

	if std::env::var("KAT_UPDATE").is_ok() {
		eprintln!("KAT pk:         {pk_hex}");
		eprintln!("KAT sig SHA256: {sig_hash}");
		return;
	}

	assert_eq!(
		pk_hex,
		"1fe9c44d840692615413d382b1d5f75d72410e190ed9d77a386f7e98bee9454e03000000"
	);
	assert_eq!(
		sig_hash,
		"d2342e8a75308f18f482f536faa886ee1d5a3d6e189b976cb3150d8e3f1f5918"
	);
}
