//! Key generation — paper §4.1.

use std::sync::Arc;

use rand::{CryptoRng, RngCore};

use crate::field::Goldilocks4;
use crate::hash::{Family, JV_RZK, JV_SEED, hash};
use crate::params::Params;
use crate::whir::code::{InterleavedCode, ReedSolomon};
use crate::whir::commitment::CodeCommitment;
use crate::whir::vc::MerkleVc;
use crate::{PublicKey, SecretKey};

/// Cached signer state held in memory after [`keygen`] for fast signing.
///
/// The cache stores the ZK-encoded commitment vector `m = (c, r_zk) ∈ F^N`,
/// where `c` is `f`'s coefficient vector (length `M`) and `r_zk` is the
/// length-`(N − M)` Prop. 3.19 encoding randomness. The WHIR commit-and-open
/// path constructs its own internal state fresh on each [`crate::sign`] call
/// from `m`. A signer that has lost the cache can rebuild it from the
/// [`SecretKey`] alone via [`SignerCache::from_secret`].
///
/// `m` is *secret* material, so on drop we zeroize the vector contents.
pub struct SignerCache {
	pub(crate) m: Vec<Goldilocks4>,
}

impl Drop for SignerCache {
	fn drop(&mut self) {
		use zeroize::Zeroize;
		self.m.zeroize();
	}
}

impl SignerCache {
	/// Rebuild the cache from the secret seed and the public-key parameters.
	pub fn from_secret(sk: &SecretKey, params: Params) -> Self {
		Self {
			m: derive_commit_vector(sk.seed(), params),
		}
	}
}

/// Generate a fresh `(PublicKey, SecretKey, SignerCache)` triple from a CSPRNG.
/// Realizes `Jevil.KeyGen` of the paper (`§3.3, Construction 4`).
///
/// `rng` is consumed only to draw a 32-byte uniform `σ`; both `c` (from
/// `JV-SEED`) and `r_zk` (from `JV-RZK`) are derived deterministically from
/// `σ`. The same `(rng-state, params)` always produces the same public key.
pub fn keygen<R: RngCore + CryptoRng>(
	rng: &mut R,
	params: Params,
) -> (PublicKey, SecretKey, SignerCache) {
	let mut sigma = [0u8; SecretKey::BYTES];
	rng.fill_bytes(&mut sigma);

	let m = derive_commit_vector(&sigma, params);
	let root = commit_root(&m, params);

	let pk = PublicKey {
		root,
		n_star: params.n_star,
	};
	let sk = SecretKey::from_bytes(sigma);
	let cache = SignerCache { m };
	(pk, sk, cache)
}

/// Derive the ZK-encoded commitment vector `m = (c, r_zk)` of length `N`
/// from a 32-byte seed. The first `M` entries are `f`'s coefficients
/// (`JV-SEED` stream); the next `N − M` are Prop. 3.19 encoding randomness
/// (`JV-RZK` stream).
pub(crate) fn derive_commit_vector(
	sigma: &[u8; SecretKey::BYTES],
	params: Params,
) -> Vec<Goldilocks4> {
	let m = params.m();
	let n = params.n();
	let mut out = Vec::with_capacity(n);
	out.extend(derive_field_elements(sigma, JV_SEED, m));
	out.extend(derive_field_elements(sigma, JV_RZK, n - m));
	out
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

/// Run WHIR's commit-only path on `m` and return the 32-byte root.
pub(crate) fn commit_root(m: &[Goldilocks4], params: Params) -> [u8; 32] {
	let n = params.n();
	assert_eq!(m.len(), n);
	const INTERLEAVING: usize = 4;
	const RATE_INV: usize = 4;
	let inner_msg_len = n / INTERLEAVING;
	let inner = ReedSolomon::<Goldilocks4>::new(inner_msg_len);
	let code = Arc::new(InterleavedCode::new(inner, INTERLEAVING));
	let vc = Arc::new(MerkleVc::new(inner_msg_len * RATE_INV));
	let cc = CodeCommitment::new(code, vc);
	let (root, _state) = cc.commit_only(m.to_vec());
	root
}

#[cfg(test)]
mod tests {
	use super::*;
	use rand::SeedableRng;
	use rand_chacha::ChaCha20Rng;

	#[test]
	fn derive_commit_vector_is_deterministic() {
		let params = Params::new(1);
		let sigma = [42u8; 32];
		let a = derive_commit_vector(&sigma, params);
		let b = derive_commit_vector(&sigma, params);
		assert_eq!(a, b);
		assert_eq!(a.len(), params.n());
	}

	#[test]
	fn commit_vector_seed_and_rzk_streams_differ() {
		let params = Params::new(3);
		let sigma = [99u8; 32];
		let m = derive_commit_vector(&sigma, params);
		let split = params.m();
		assert_ne!(&m[..1], &m[split..split + 1]);
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
		assert_eq!(cache.m, cache2.m);
		assert_eq!(pk.root, commit_root(&cache2.m, params));
	}
}
