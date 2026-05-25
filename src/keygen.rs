//! Key generation — paper §4.1.

use std::sync::Arc;

use rand::{CryptoRng, RngCore};

use crate::field::Goldilocks4;
use crate::hash::{Family, JV_SEED, hash};
use crate::params::Params;
use crate::whir::code::{InterleavedCode, ReedSolomon};
use crate::whir::commitment::CodeCommitment;
use crate::whir::vc::MerkleVc;
use crate::{PublicKey, SecretKey};

/// Cached signer state held in memory after [`keygen`] for fast signing.
///
/// The cache stores the coefficient vector `c`. The WHIR commit-and-open path
/// constructs its own internal state fresh on each [`crate::sign`] call from
/// `c`. A signer that has lost the cache can rebuild it from the
/// [`SecretKey`] alone via [`SignerCache::from_secret`] — at the cost of
/// re-running the deterministic-derivation step.
///
/// `c` is *secret* material (it is `f`'s coefficient vector), so on drop we
/// zeroize the vector contents via [`zeroize::Zeroize`].
pub struct SignerCache {
	pub(crate) c: Vec<Goldilocks4>,
}

impl Drop for SignerCache {
	fn drop(&mut self) {
		use zeroize::Zeroize;
		self.c.zeroize();
	}
}

impl SignerCache {
	/// Rebuild the cache from the secret seed and the public-key parameters.
	/// Produces the exact same `c` (and therefore the same `PublicKey`) as
	/// the original [`keygen`] call.
	pub fn from_secret(sk: &SecretKey, params: Params) -> Self {
		Self {
			c: derive_c(sk.seed(), params),
		}
	}
}

/// Generate a fresh `(PublicKey, SecretKey, SignerCache)` triple from a
/// cryptographically-strong RNG.
///
/// `rng` is consumed only to draw a single 32-byte uniform `σ`; all
/// subsequent randomness is derived deterministically from `σ` via the
/// `JV-SEED` SHAKE256 stream. The same `(rng-state, params)` always produces
/// the same public key — useful for testing, but consequential for
/// production: re-seeding the RNG identically will re-derive the same
/// signing key.
pub fn keygen<R: RngCore + CryptoRng>(
	rng: &mut R,
	params: Params,
) -> (PublicKey, SecretKey, SignerCache) {
	let mut sigma = [0u8; SecretKey::BYTES];
	rng.fill_bytes(&mut sigma);

	let c = derive_c(&sigma, params);
	let root = commit_c_root(&c, params);

	let pk = PublicKey {
		root,
		n_star: params.n_star,
	};
	let sk = SecretKey::from_bytes(sigma);
	let cache = SignerCache { c };
	(pk, sk, cache)
}

/// Derive the coefficient vector `c = (c_0, …, c_{M-1})` from `σ`.
pub(crate) fn derive_c(sigma: &[u8; SecretKey::BYTES], params: Params) -> Vec<Goldilocks4> {
	derive_field_elements(sigma, JV_SEED, params.m())
}

/// Pull `count` uniform `Goldilocks4` elements from the `SHAKE256(tag ‖ σ)`
/// stream with per-limb rejection sampling.
fn derive_field_elements(
	sigma: &[u8; SecretKey::BYTES],
	tag: [u8; 8],
	count: usize,
) -> Vec<Goldilocks4> {
	if count == 0 {
		return Vec::new();
	}
	let mut buffer_size = count * 32 * 2 + 32;
	let mut refill_tag = 0u64;
	loop {
		let extra = refill_tag.to_le_bytes();
		let stream = if refill_tag == 0 {
			hash(Family::Xof, tag, &[sigma], buffer_size)
		} else {
			hash(Family::Xof, tag, &[sigma, &extra], buffer_size)
		};
		let mut out = Vec::with_capacity(count);
		let mut cursor = 0usize;
		while out.len() < count && cursor + 32 <= stream.len() {
			let chunk = &stream[cursor..cursor + 32];
			cursor += 32;
			if let Some(g) = Goldilocks4::from_bytes(chunk) {
				out.push(g);
			}
		}
		if out.len() == count {
			return out;
		}
		buffer_size *= 2;
		refill_tag += 1;
	}
}

/// Run WHIR's commit-only path on `c` and return just the 32-byte root.
pub(crate) fn commit_c_root(c: &[Goldilocks4], params: Params) -> [u8; 32] {
	let m = params.m();
	assert_eq!(c.len(), m);
	const INTERLEAVING: usize = 4;
	const RATE_INV: usize = 4;
	let inner_msg_len = m / INTERLEAVING;
	let inner = ReedSolomon::<Goldilocks4>::new(inner_msg_len);
	let code = Arc::new(InterleavedCode::new(inner, INTERLEAVING));
	let vc = Arc::new(MerkleVc::new(inner_msg_len * RATE_INV));
	let cc = CodeCommitment::new(code, vc);
	let (root, _state) = cc.commit_only(c.to_vec());
	root
}

#[cfg(test)]
mod tests {
	use super::*;
	use rand::SeedableRng;
	use rand_chacha::ChaCha20Rng;

	#[test]
	fn derive_c_is_deterministic() {
		let params = Params::new(1);
		let sigma = [42u8; 32];
		let a = derive_c(&sigma, params);
		let b = derive_c(&sigma, params);
		assert_eq!(a, b);
		assert_eq!(a.len(), params.m());
	}

	#[test]
	fn keygen_is_deterministic_under_seeded_rng() {
		let params = Params::new(1);
		let mut a = ChaCha20Rng::seed_from_u64(0);
		let mut b = ChaCha20Rng::seed_from_u64(0);
		let (pk_a, _, _) = keygen(&mut a, params);
		let (pk_b, _, _) = keygen(&mut b, params);
		assert_eq!(pk_a.root, pk_b.root);
		assert_eq!(pk_a.n_star, pk_b.n_star);
	}

	#[test]
	fn signer_cache_from_secret_matches_keygen() {
		let params = Params::new(1);
		let mut rng = ChaCha20Rng::seed_from_u64(1);
		let (pk, sk, cache) = keygen(&mut rng, params);
		let cache2 = SignerCache::from_secret(&sk, params);
		assert_eq!(cache.c, cache2.c);
		assert_eq!(pk.root, commit_c_root(&cache2.c, params));
	}
}
