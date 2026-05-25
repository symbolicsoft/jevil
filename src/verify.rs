//! Verification — paper §4.3.

use spongefish::domain_separator;

use crate::alpha::BatchedAlpha;
use crate::field::{Goldilocks4, psi};
use crate::keygen::derive_ood_point;
use crate::params::Params;
use crate::positions::derive_positions;
use crate::sign::Signature;
use crate::transcript::{derive_betas, prefix_bytes};
use crate::whir::ConcreteWhirVerifier;
use crate::{Error, PublicKey};

/// Verify a Jevil signature. Realizes `Jevil.Verify` of the paper
/// (`§3.3, Construction 6`).
///
/// Returns `Ok(())` iff the signature is valid for `(pk, params, msg)`.
/// Returns one of:
///
/// - [`Error::ParamsMismatch`] if `params.n_star != pk.n_star`,
/// - [`Error::InvalidLength`] if `sig.y_values.len() != Params::K`,
/// - [`Error::VerificationFailed`] for any cryptographic failure (tampered
///   y-values, malformed proof, wrong message, …).
///
/// [`Error::NonCanonicalField`] is **not** reachable from `verify`: a
/// `Goldilocks4` cannot exist non-canonically, so the only failure mode of
/// that kind is at signature *parse* time in [`Signature::from_bytes`].
///
/// `verify` does **not** distinguish *which* check failed — all
/// cryptographic-failure paths collapse to the single
/// [`Error::VerificationFailed`] variant to deny side channels.
pub fn verify(pk: &PublicKey, params: Params, msg: &[u8], sig: &Signature) -> Result<(), Error> {
	if pk.n_star != params.n_star {
		return Err(Error::ParamsMismatch);
	}
	let k = Params::K as usize;
	if sig.y_values.len() != k {
		return Err(Error::InvalidLength);
	}
	// Canonicality of `sig.y_values` is structural: a `Goldilocks4` cannot
	// exist with a non-canonical limb (parser enforces it in
	// `Signature::from_bytes`; `Goldilocks::new` reduces). No re-check here.

	// 1. Re-derive positions and x_t; re-derive the OOD point z from root
	//    (identical derivation to KeyGen — the verifier never receives z
	//    on the wire).
	let positions = derive_positions(&pk.root, msg, k, params.t());
	let mut xs: Vec<Goldilocks4> = positions
		.iter()
		.map(|&i| psi(i as u64, params.t() as u64))
		.collect();
	let z = derive_ood_point(&pk.root);
	xs.push(z);

	// 2. Re-derive K+1 betas (β_1..β_K for the message-derived points plus
	//    β_{K+1} for the OOD term) and compute the batched target
	//    v = Σ_{t≤K} β_t·y_t + β_{K+1}·w.
	let betas = derive_betas(&pk.root, msg, &sig.y_values);
	if betas.len() != xs.len() {
		return Err(Error::VerificationFailed);
	}
	let mut v: Goldilocks4 = betas
		.iter()
		.zip(sig.y_values.iter())
		.map(|(b, y)| *b * *y)
		.sum();
	v += betas[k] * pk.w;

	// 3. Reconstruct the spongefish transcript using the same instance bytes
	//    the signer used. The opaque whir_proof IS the narg_string in full.
	let prefix = prefix_bytes(params, &pk.root, &pk.w, msg, &sig.y_values);
	let domain = domain_separator!("jevil-v1")
		.without_session()
		.instance(&prefix);
	let mut transcript = domain.std_verifier(&sig.whir_proof);

	// 4. Build the symbolic length-`M` α handle from K+1 (x, β) pairs —
	//    α = Σ_{t≤K} β_t·u(x_t) + β_{K+1}·u(z). O(K · ν) verifier with no
	//    length-M alloc; the embedding into WHIR's length-`N` wire format
	//    happens inside `whir.verify`.
	let alpha = BatchedAlpha::new(&xs, betas, params.nu());

	// 5. Run WHIR's verifier on top.
	let hvzk_budget = params.n() - params.m();
	let whir = ConcreteWhirVerifier::build(params.m(), hvzk_budget, 64, 64);
	whir.verify(&mut transcript, alpha, v)
		.map_err(|_| Error::VerificationFailed)?;
	transcript
		.check_eof()
		.map_err(|_| Error::VerificationFailed)?;
	Ok(())
}

#[cfg(test)]
mod tests {
	use super::*;
	use crate::keygen::keygen;
	use crate::sign::sign;
	use rand::SeedableRng;
	use rand_chacha::ChaCha20Rng;

	#[test]
	fn honest_signature_verifies() {
		let params = Params::new(1);
		let mut rng = ChaCha20Rng::seed_from_u64(0);
		let (pk, sk, cache) = keygen(&mut rng, params);
		let sig = sign(&sk, &pk, &cache, params, b"hello");
		assert!(verify(&pk, params, b"hello", &sig).is_ok());
	}

	#[test]
	fn tampered_y_value_rejected() {
		let params = Params::new(1);
		let mut rng = ChaCha20Rng::seed_from_u64(1);
		let (pk, sk, cache) = keygen(&mut rng, params);
		let mut sig = sign(&sk, &pk, &cache, params, b"hi");
		sig.y_values[0] += Goldilocks4::ONE;
		assert_eq!(
			verify(&pk, params, b"hi", &sig),
			Err(Error::VerificationFailed)
		);
	}

	#[test]
	fn wrong_message_rejected() {
		let params = Params::new(1);
		let mut rng = ChaCha20Rng::seed_from_u64(2);
		let (pk, sk, cache) = keygen(&mut rng, params);
		let sig = sign(&sk, &pk, &cache, params, b"a");
		assert!(verify(&pk, params, b"b", &sig).is_err());
	}
}
