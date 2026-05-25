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
		"dcc49c81854fad472af87a621312de0c402500d388ab9ecb485e6eddfb36be608f286c98daeeb93ad81c07bbd5027f059b0bd0f6465995ac6e867d77c5a9398103000000"
	);
	assert_eq!(
		sig_hash,
		"071dd05343ab9ac1395f21cb2a4fcc622c49c77fa515a47f794ca9dda53279ad"
	);
}
