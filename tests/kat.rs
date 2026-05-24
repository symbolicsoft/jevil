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
		"e24d20f1b35319ac0b411fc90dfcb5f2dcf1983af95ec43a8e89735e4e63bc7d03000000"
	);
	assert_eq!(
		sig_hash,
		"f48ca354185b3945d1195f83085a2f099e8e211824cc4894e6e2505c987e9cf0"
	);
}
