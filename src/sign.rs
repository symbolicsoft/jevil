//! Signing — paper §4.2.

use spongefish::domain_separator;

use crate::field::{Goldilocks4, psi};
use crate::keygen::SignerCache;
use crate::params::Params;
use crate::positions::derive_positions;
use crate::transcript::{derive_betas, prefix_bytes};
use crate::whir::{ConcreteWhirProtocol, LinearForm};
use crate::{Error, PublicKey, SecretKey};

/// A Jevil signature.
///
/// Layout: the `K` revealed evaluations `y_1, …, y_K` (32 bytes each)
/// followed by the inline WHIR proof.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Signature {
	/// The revealed `y_t = f(x_t)` values, one per sampled position.
	pub y_values: Vec<Goldilocks4>,
	/// The WHIR linear-form proof. Opaque blob — its layout is an
	/// implementation detail of the WHIR sub-protocol.
	pub whir_proof: Vec<u8>,
}

impl Signature {
	/// Serialise as `y_1.to_bytes() ‖ … ‖ y_K.to_bytes() ‖ whir_proof`.
	pub fn to_bytes(&self) -> Vec<u8> {
		let mut out = Vec::with_capacity(self.y_values.len() * 32 + self.whir_proof.len());
		for y in &self.y_values {
			out.extend_from_slice(&y.to_bytes());
		}
		out.extend_from_slice(&self.whir_proof);
		out
	}

	/// Parse exactly `K = Params::K` field elements followed by the WHIR proof.
	///
	/// Returns [`Error::InvalidLength`] if the input is too short, or
	/// [`Error::NonCanonicalField`] if any `y_t` chunk is not a canonical
	/// extension-field element.
	pub fn from_bytes(b: &[u8]) -> Result<Self, Error> {
		let k = Params::K as usize;
		let head_len = k * 32;
		if b.len() < head_len {
			return Err(Error::InvalidLength);
		}
		let mut y_values = Vec::with_capacity(k);
		for chunk in b[..head_len].chunks_exact(32) {
			let g = Goldilocks4::from_bytes(chunk).ok_or(Error::NonCanonicalField)?;
			y_values.push(g);
		}
		let whir_proof = b[head_len..].to_vec();
		Ok(Self {
			y_values,
			whir_proof,
		})
	}
}

/// Produce a Jevil signature on `msg` under `(pk, sk)`. Realizes
/// `Jevil.Sign` of the paper (`§3.3, Construction 5`).
///
/// `cache` must be the [`SignerCache`] returned by [`crate::keygen`] (or
/// rebuilt via [`SignerCache::from_secret`]) for this `sk`; passing a cache
/// from a different signer produces undefined output.
///
/// `sk` is consumed to derive the per-signature prover-randomness seed
/// `ρ = H_xof(JV-OPRD, s, root, msg, y_1, …, y_K; 32)` per paper §2.2 /
/// Construction 5 step 6: that seed deterministically drives every
/// internally-sampled value inside `WHIR.Open` (sumcheck masks,
/// code-switching masks, OOD answers), so `Sign` is a pure function of
/// `(sk, pk, msg)`.
pub fn sign(
	sk: &SecretKey,
	pk: &PublicKey,
	cache: &SignerCache,
	params: Params,
	msg: &[u8],
) -> Signature {
	let k = Params::K as usize;
	let t = params.t();
	let m = params.m();
	let n = params.n();

	// 1. Positions and their ψ-images.
	let positions = derive_positions(&pk.root, msg, k, t);
	let xs: Vec<Goldilocks4> = positions.iter().map(|&i| psi(i as u64, t as u64)).collect();

	// 2. y_t = f(x_t) via Horner over the M coefficients held in the first
	//    M slots of `cache.m`.
	let coeffs = &cache.m[..m];
	let ys: Vec<Goldilocks4> = xs.iter().map(|&x| horner(coeffs, x)).collect();

	// 3. β challenges (verifier re-derives the same vector from
	//    `(root, msg, ys)`).
	let betas = derive_betas(&pk.root, msg, &ys);

	// 4. The lift `α = Σ_t β_t · u(x_t)` lives in `F^M`; we materialise it
	//    here, embedded into the WHIR primitive's length-`N` wire format
	//    (zero-padded over the WHIR encoding-randomness slots — see
	//    `lift.rs` module docs for the embedding). Computed via a parallel
	//    Horner pass over the K positions to avoid the K transient
	//    length-N allocations the old MonomialLift path used.
	let mut alpha = vec![Goldilocks4::ZERO; n];
	let mut x_powers = vec![Goldilocks4::ONE; xs.len()];
	for slot in alpha.iter_mut().take(m) {
		let mut sum = Goldilocks4::ZERO;
		for t in 0..xs.len() {
			sum += betas[t] * x_powers[t];
		}
		*slot = sum;
		for t in 0..xs.len() {
			x_powers[t] *= xs[t];
		}
	}
	// alpha[m..n] stays zero by initialisation (WHIR-pad region).

	// 5. Build the Fiat–Shamir transcript with the deterministic prefix
	//    binding (params, root, msg, ys) into its instance bytes — then run
	//    WHIR's prover on top.
	let prefix = prefix_bytes(params, &pk.root, msg, &ys);
	let domain = domain_separator!("jevil-v1")
		.without_session()
		.instance(&prefix);
	let mut transcript = domain.std_prover();

	// Derive the per-signature prover-randomness seed ρ from
	// (sk_seed, root, msg, y_1, …, y_K) under the JV-OPRD domain tag
	// (paper §2.2 / Construction 5 step 6). Deterministic but unique
	// per signature; the seed drives every internally-sampled value
	// inside `WHIR.Open` (sumcheck masks via Construction 6.3,
	// code-switching masks via Construction 9.7, OOD answers via
	// Lemma 9.3).
	let mask_seed = derive_prover_randomness_seed(sk, &pk.root, msg, &ys);

	let whir = ConcreteWhirProtocol::build(n, 32, 64);
	whir.prove_to_transcript(
		&mut transcript,
		cache.m.clone(),
		LinearForm::new(alpha),
		&mask_seed,
	);

	Signature {
		y_values: ys,
		whir_proof: transcript.narg_string().to_vec(),
	}
}

/// Derive the per-signature prover-randomness seed `ρ` from
/// `(sk_seed, root, msg, y_1, …, y_K)` per paper §2.2 / Construction 5
/// step 6.
///
/// The hash inputs match the spec literally: the seed `s`, root, message,
/// and each revealed evaluation `y_t` are passed as **separate**
/// length-prefixed inputs to the length-prefix framing of [`crate::hash`],
/// so the framing is `tag ‖ len_8(s) ‖ s ‖ len_8(root) ‖ root ‖
/// len_8(msg) ‖ msg ‖ len_8(y_1) ‖ y_1 ‖ … ‖ len_8(y_K) ‖ y_K`. This
/// matches the spec's `H_xof(JV-OPRD, s, root, M, y_1, …, y_K; ∞)`
/// exactly, except we truncate the XOF stream at 32 bytes and re-expand
/// downstream via `derive_field_vec` (which is itself JV-OPRD-tagged) —
/// security-equivalent in the random-oracle model.
///
/// Deterministic — the same `(sk, pk, msg)` produces the same seed —
/// but unique across messages and `y_t` tuples, which is what HVZK
/// requires per signature.
fn derive_prover_randomness_seed(
	sk: &SecretKey,
	root: &[u8; 32],
	msg: &[u8],
	ys: &[Goldilocks4],
) -> [u8; 32] {
	use crate::hash::{Family, JV_OPRD, hash};
	let y_bytes: Vec<[u8; 32]> = ys.iter().map(|y| y.to_bytes()).collect();
	let mut inputs: Vec<&[u8]> = Vec::with_capacity(3 + ys.len());
	inputs.push(sk.seed());
	inputs.push(root);
	inputs.push(msg);
	for yb in &y_bytes {
		inputs.push(yb);
	}
	let h = hash(Family::Xof, JV_OPRD, &inputs, 32);
	let mut out = [0u8; 32];
	out.copy_from_slice(&h);
	out
}

/// Horner evaluation of `Σ_k coeffs[k] · x^k`.
fn horner(coeffs: &[Goldilocks4], x: Goldilocks4) -> Goldilocks4 {
	let mut acc = Goldilocks4::ZERO;
	for c in coeffs.iter().rev() {
		acc = acc * x + *c;
	}
	acc
}

#[cfg(test)]
mod tests {
	use super::*;
	use crate::keygen::keygen;
	use rand::SeedableRng;
	use rand_chacha::ChaCha20Rng;

	#[test]
	fn signature_has_k_y_values_and_nonempty_proof() {
		let params = Params::new(1);
		let mut rng = ChaCha20Rng::seed_from_u64(0);
		let (pk, sk, cache) = keygen(&mut rng, params);
		let sig = sign(&sk, &pk, &cache, params, b"hello");
		assert_eq!(sig.y_values.len(), Params::K as usize);
		assert!(!sig.whir_proof.is_empty());
	}

	#[test]
	fn signature_serde_round_trip() {
		let params = Params::new(1);
		let mut rng = ChaCha20Rng::seed_from_u64(2);
		let (pk, sk, cache) = keygen(&mut rng, params);
		let sig = sign(&sk, &pk, &cache, params, b"round");
		let bytes = sig.to_bytes();
		let parsed = Signature::from_bytes(&bytes).unwrap();
		assert_eq!(parsed.y_values, sig.y_values);
		assert_eq!(parsed.whir_proof, sig.whir_proof);
	}
}
