//! Key generation — paper §4.1.

use std::sync::Arc;

use rand::{CryptoRng, RngCore};

use crate::field::Goldilocks4;
use crate::hash::{Family, JV_MASK, JV_SEED, hash};
use crate::params::Params;
use crate::whir::code::{InterleavedCode, ReedSolomon};
use crate::whir::commitment::CodeCommitment;
use crate::whir::vc::MerkleVc;
use crate::{PublicKey, SecretKey};

/// Cached signer state held in memory after [`keygen`] for fast signing.
///
/// The cache stores the padded coefficient vector `c^pad`. The WHIR
/// commit-and-open path constructs its own internal state fresh on each
/// [`crate::sign`] call from `c^pad`. A signer that has lost the cache can
/// rebuild it from the [`SecretKey`] alone via [`SignerCache::from_secret`] —
/// at the cost of re-running the deterministic-derivation step.
///
/// `c^pad` is *secret* material (its first `M` entries are `f`'s coefficients),
/// so on drop we zeroize the vector contents via [`zeroize::Zeroize`].
pub struct SignerCache {
	pub(crate) c_pad: Vec<Goldilocks4>,
}

impl Drop for SignerCache {
	fn drop(&mut self) {
		use zeroize::Zeroize;
		self.c_pad.zeroize();
	}
}

impl SignerCache {
	/// Rebuild the cache from the secret seed and the public-key parameters.
	/// Produces the exact same `c^pad` (and therefore the same `PublicKey`)
	/// as the original [`keygen`] call.
	pub fn from_secret(sk: &SecretKey, params: Params) -> Self {
		Self {
			c_pad: derive_c_pad(sk.seed(), params),
		}
	}
}

/// Generate a fresh `(PublicKey, SecretKey, SignerCache)` triple from a
/// cryptographically-strong RNG.
///
/// `rng` is consumed only to draw a single 32-byte uniform `σ`; all
/// subsequent randomness is derived deterministically from `σ` via the
/// `JV-SEED` / `JV-MASK` SHAKE256 streams. The same `(rng-state, params)`
/// always produces the same public key — useful for testing, but consequential
/// for production: re-seeding the RNG identically will re-derive the same
/// signing key.
pub fn keygen<R: RngCore + CryptoRng>(
	rng: &mut R,
	params: Params,
) -> (PublicKey, SecretKey, SignerCache) {
	let mut sigma = [0u8; SecretKey::BYTES];
	rng.fill_bytes(&mut sigma);

	let c_pad = derive_c_pad(&sigma, params);
	let root = commit_c_pad_root(&c_pad, params);

	let pk = PublicKey {
		root,
		n_star: params.n_star,
	};
	let sk = SecretKey::from_bytes(sigma);
	let cache = SignerCache { c_pad };
	(pk, sk, cache)
}

/// Compose the padded coefficient vector
/// `c^pad = (c_0, …, c_{M-1}, r_1, …, r_{N-M})` from `σ`.
pub(crate) fn derive_c_pad(sigma: &[u8; SecretKey::BYTES], params: Params) -> Vec<Goldilocks4> {
	let m = params.m();
	let n = params.n();
	let mut c_pad = Vec::with_capacity(n);
	c_pad.extend(derive_field_elements(sigma, JV_SEED, m));
	c_pad.extend(derive_field_elements(sigma, JV_MASK, n - m));
	c_pad
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
		buffer_size = (buffer_size * 2).max(64);
		refill_tag += 1;
	}
}

/// Run WHIR's commit-only path on `c^pad` and return just the 32-byte root.
pub(crate) fn commit_c_pad_root(c_pad: &[Goldilocks4], params: Params) -> [u8; 32] {
	let n = params.n();
	assert_eq!(c_pad.len(), n);
	const INTERLEAVING: usize = 4;
	const RATE_INV: usize = 4;
	let inner_msg_len = n / INTERLEAVING;
	let inner = ReedSolomon::<Goldilocks4>::new(inner_msg_len);
	let code = Arc::new(InterleavedCode::new(inner, INTERLEAVING));
	let vc = Arc::new(MerkleVc::new(inner_msg_len * RATE_INV));
	let cc = CodeCommitment::new(code, vc);
	let (root, _state) = cc.commit_only(c_pad.to_vec());
	root
}

#[cfg(test)]
mod tests {
	use super::*;
	use rand::SeedableRng;
	use rand_chacha::ChaCha20Rng;

	#[test]
	fn derive_c_pad_is_deterministic() {
		let params = Params::new(1);
		let sigma = [42u8; 32];
		let a = derive_c_pad(&sigma, params);
		let b = derive_c_pad(&sigma, params);
		assert_eq!(a, b);
		assert_eq!(a.len(), params.n());
	}

	#[test]
	fn derive_c_pad_seed_and_mask_streams_differ() {
		let params = Params::new(1);
		let sigma = [99u8; 32];
		let c_pad = derive_c_pad(&sigma, params);
		let m = params.m();
		assert_ne!(&c_pad[..1], &c_pad[m..m + 1]);
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
		assert_eq!(cache.c_pad, cache2.c_pad);
		assert_eq!(pk.root, commit_c_pad_root(&cache2.c_pad, params));
	}
}
