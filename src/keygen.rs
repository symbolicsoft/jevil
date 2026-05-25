//! Key generation — paper §4.1.

use rand::{CryptoRng, RngCore};

use crate::field::Goldilocks4;
use crate::hash::{Family, JV_SEED, hash};
use crate::params::Params;
use crate::whir::{ConcreteWhirProtocol, WhirSignerState};
use crate::{PublicKey, SecretKey};

/// Cached signer state held in memory after [`keygen`] for fast signing.
///
/// The cache stores the coefficient vector `c ∈ F^M` plus the opaque
/// [`WhirSignerState`] produced by `WHIR.Commit`. The state holds whatever
/// the WHIR primitive needs to drive `WHIR.Open` (the Prop. 3.19 encoding
/// randomness, internalised inside the primitive); jevil never names it.
///
/// A signer that has lost the cache can rebuild it from the
/// [`SecretKey`] alone via [`SignerCache::from_secret`].
///
/// `c` is *secret* material, so on drop we zeroize the vector contents
/// (the internal WHIR state is also zeroised on drop).
pub struct SignerCache {
	pub(crate) c: Vec<Goldilocks4>,
	pub(crate) whir_state: WhirSignerState,
}

impl Drop for SignerCache {
	fn drop(&mut self) {
		use zeroize::Zeroize;
		self.c.zeroize();
		self.whir_state.internal.zeroize();
	}
}

impl SignerCache {
	/// Rebuild the cache from the secret seed and the public-key parameters.
	pub fn from_secret(sk: &SecretKey, params: Params) -> Self {
		let c = derive_coefficient_vector(sk.seed(), params);
		let whir = build_whir_protocol(params);
		let (_root, state) = whir.commit(&c, sk.seed());
		Self {
			c,
			whir_state: state,
		}
	}
}

/// Generate a fresh `(PublicKey, SecretKey, SignerCache)` triple from a CSPRNG.
/// Realizes `Jevil.KeyGen` of the paper (`§3.3, Construction 4`).
///
/// `rng` is consumed only to draw a 32-byte uniform `σ`. The polynomial
/// coefficients `c` are derived from `σ` via `JV-SEED`, and `WHIR.Commit`
/// uses the same `σ` to deterministically derive its internal Prop. 3.19
/// encoding randomness via `JV-RZK`. The same `(rng-state, params)`
/// always produces the same public key.
pub fn keygen<R: RngCore + CryptoRng>(
	rng: &mut R,
	params: Params,
) -> (PublicKey, SecretKey, SignerCache) {
	let mut sigma = [0u8; SecretKey::BYTES];
	rng.fill_bytes(&mut sigma);

	let c = derive_coefficient_vector(&sigma, params);
	let whir = build_whir_protocol(params);
	let (root, whir_state) = whir.commit(&c, &sigma);

	let pk = PublicKey {
		root,
		n_star: params.n_star,
	};
	let sk = SecretKey::from_bytes(sigma);
	let cache = SignerCache { c, whir_state };
	(pk, sk, cache)
}

/// Derive the polynomial coefficient vector `c ∈ F^M` from a 32-byte
/// seed via `JV-SEED`.
pub(crate) fn derive_coefficient_vector(
	sigma: &[u8; SecretKey::BYTES],
	params: Params,
) -> Vec<Goldilocks4> {
	derive_field_elements(sigma, JV_SEED, params.m())
}

pub(crate) fn build_whir_protocol(params: Params) -> ConcreteWhirProtocol {
	let hvzk_budget = params.n() - params.m();
	ConcreteWhirProtocol::build(params.m(), hvzk_budget, 32, 64)
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

#[cfg(test)]
mod tests {
	use super::*;
	use rand::SeedableRng;
	use rand_chacha::ChaCha20Rng;

	#[test]
	fn derive_coefficient_vector_is_deterministic() {
		let params = Params::new(1);
		let sigma = [42u8; 32];
		let a = derive_coefficient_vector(&sigma, params);
		let b = derive_coefficient_vector(&sigma, params);
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
		assert_eq!(cache.whir_state.internal, cache2.whir_state.internal);
		let whir = build_whir_protocol(params);
		let (rebuilt_root, _) = whir.commit(&cache2.c, sk.seed());
		assert_eq!(pk.root, rebuilt_root);
	}
}
