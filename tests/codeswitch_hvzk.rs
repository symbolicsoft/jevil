//! HVZK code-switching (Construction 9.7, Lemma 9.8 of ePrint 2026/391) —
//! integration tests exercising the full codeswitch path: per-round padding
//! mask commit, privacy-padded OOD answers, mask-stack threading.
//!
//! Round-trip at `n_star = 7` exercises ≥ 1 codeswitch round; `n_star = 31`
//! exercises ≥ 2 codeswitch rounds (and so the IOR's accumulator across
//! multiple recursive folds). The strong simulator-vs-real witness is the
//! Theorem 4.5 composed simulator in `tests/multi_opening_hvzk.rs`.

use jevil::{Params, keygen, sign, verify};
use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;

#[test]
fn codeswitch_round_trip_at_n_star_7() {
	let params = Params::new(7);
	let mut rng = ChaCha20Rng::seed_from_u64(0);
	let (pk, sk, cache) = keygen(&mut rng, params);
	let sig = sign(&sk, &pk, &cache, params, b"codeswitch-hvzk");
	verify(&pk, params, b"codeswitch-hvzk", &sig).expect("honest codeswitch must verify");
}

#[test]
fn codeswitch_round_trip_at_n_star_31() {
	let params = Params::new(31);
	let mut rng = ChaCha20Rng::seed_from_u64(0);
	let (pk, sk, cache) = keygen(&mut rng, params);
	let sig = sign(&sk, &pk, &cache, params, b"codeswitch-deep");
	verify(&pk, params, b"codeswitch-deep", &sig).expect("deep codeswitch must verify");
}
