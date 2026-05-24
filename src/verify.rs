//! Verification — paper §4.3.

use spongefish::domain_separator;

use crate::alpha::BatchedAlpha;
use crate::field::{Goldilocks4, psi};
use crate::params::Params;
use crate::positions::derive_positions;
use crate::sign::Signature;
use crate::transcript::{derive_betas, prefix_bytes};
use crate::whir::{ConcreteWhirVerifier, LinearConstraint};
use crate::{Error, PublicKey};

/// Verify a Jevil signature.
///
/// Returns `Ok(())` iff the signature is valid for `(pk, params, msg)`.
/// Returns one of:
///
/// - [`Error::ParamsMismatch`] if `params.n_star != pk.n_star`,
/// - [`Error::InvalidLength`] if `sig.y_values.len() != Params::K`,
/// - [`Error::NonCanonicalField`] if a y-value isn't a canonical extension
///   element,
/// - [`Error::VerificationFailed`] for any cryptographic failure (tampered
///   y-values, malformed proof, wrong message, …).
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
	for y in &sig.y_values {
		if Goldilocks4::from_bytes(&y.to_bytes()).is_none() {
			return Err(Error::NonCanonicalField);
		}
	}

	// 1. Re-derive positions, x_t, β_t, and v = Σ β_t · y_t.
	let positions = derive_positions(&pk.root, msg, k, params.t());
	let xs: Vec<Goldilocks4> = positions
		.iter()
		.map(|&i| psi(i as u64, params.t() as u64))
		.collect();
	let betas = derive_betas(&pk.root, msg, &sig.y_values);
	let v: Goldilocks4 = betas
		.iter()
		.zip(sig.y_values.iter())
		.map(|(b, y)| *b * *y)
		.sum();

	// 2. Reconstruct the spongefish transcript using the same instance bytes
	//    the signer used. The opaque whir_proof IS the narg_string in full.
	let prefix = prefix_bytes(params, &pk.root, msg, &sig.y_values);
	let domain = domain_separator!("jevil-v1")
		.without_session()
		.instance(&prefix);
	let mut transcript = domain.std_verifier(&sig.whir_proof);

	// 3. Build the symbolic α handle (O(K · ν') verifier — no length-N alloc).
	let alpha = BatchedAlpha::new(&xs, betas, params.nu(), params.nu_prime());
	let constraint = LinearConstraint::new(alpha, v);

	// 4. Run WHIR's verifier on top.
	let whir = ConcreteWhirVerifier::build(params.n(), 32, 64);
	whir.verify_from_transcript(&mut transcript, constraint)
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
