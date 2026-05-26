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
		"f29166e4db9da85fb8e3a03cbb312b5453a83080494cabf1d4aa1b5109c5cca4c9f92020d15fa3fade124e3c70e38d3e9cf66a64af1f3a6d98b18da2c9cef49703000000"
	);
	assert_eq!(
		sig_hash,
		"c628f3f89c5370a56d2ee5cc28e7fe047d6d2f4367a69fe94f09783b435f6159"
	);
}
