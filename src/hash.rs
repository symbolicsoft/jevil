//! Domain-tagged hashing.
//!
//! Jevil uses two random-oracle-modelled primitives:
//!
//! - [`Family::Arith`]: **Poseidon2-Goldilocks-12** (state width 12, rate 8,
//!   capacity 4, S-box `x⁷`). Arithmetic-friendly, used *only* inside WHIR for
//!   its codeword vector commitment.
//! - [`Family::Xof`]: **SHAKE256** (extendable-output). Used for seed
//!   expansion, position derivation, Fiat–Shamir, and per-signature
//!   prover-randomness derivation.
//!
//! Every hash invocation is prefixed by an 8-byte ASCII *domain tag* that
//! separates the eight logical uses (paper §2.2). The tags are exposed as
//! module-level constants below.
//!
//! ## Serialisation format
//!
//! `hash(family, tag, [x₁, x₂, …, xₖ]; L)` returns the first `L` output bytes
//! of `family` applied to
//!
//! ```text
//! tag ‖ len_8(x₁) ‖ x₁ ‖ len_8(x₂) ‖ x₂ ‖ … ‖ len_8(xₖ) ‖ xₖ
//! ```
//!
//! where `len_8(x)` is the byte length of `x` encoded as an 8-byte
//! little-endian unsigned integer. The length prefix is what makes the
//! framing *injective*: concatenating two inputs of different lengths can
//! never produce the same serialised buffer as concatenating two of any other
//! lengths.

use p3_field::{PrimeCharacteristicRing, PrimeField64};
use p3_goldilocks::{Goldilocks, default_goldilocks_poseidon2_12};
use p3_symmetric::Permutation;
use shake::{ExtendableOutput, Shake256, Update, XofReader};

// ---------------------------------------------------------------------------
// Domain tags (paper §2.2)
// ---------------------------------------------------------------------------

/// Domain tag for the seed-derived polynomial coefficients (XOF).
pub(crate) const JV_SEED: [u8; 8] = *b"JV-SEED ";
/// Domain tag for the seed-derived ZK encoding randomness (Prop. 3.19 of
/// eprint 2026/391). Used **only** at `KeyGen` to extend the committed
/// message from `M` to `N` before NTT encoding so that any subset of
/// ≤ `N − M` codeword positions reveals nothing about the honest
/// coefficients. Per-signature prover-side randomness uses [`JV_OPRD`].
pub(crate) const JV_RZK: [u8; 8] = *b"JV-RZK  ";
/// Domain tag for per-message position derivation (XOF).
pub(crate) const JV_POSN: [u8; 8] = *b"JV-POSN ";
/// Domain tag for the Fiat–Shamir batching challenges (XOF).
pub(crate) const JV_FSCH: [u8; 8] = *b"JV-FSCH ";
/// Domain tag for WHIR's internal vector commitment (Poseidon2).
pub(crate) const JV_WHIR: [u8; 8] = *b"JV-WHIR ";
/// Domain tag for the WHIR transcript's instance-bytes prefix
/// ([`crate::transcript::prefix_bytes`]). Not used as input to a hash
/// invocation — consumed only by the Fiat–Shamir layer as the leading
/// 8 bytes of the binding prefix per paper §4.2 step 4 / §4.3 step 3.
pub(crate) const JV_OPEN: [u8; 8] = *b"JV-OPEN ";
/// Domain tag for the per-signature prover-randomness derivation (XOF).
/// Used by [`crate::sign`] to derive a deterministic seed from
/// `(s, root, msg, y_1, …, y_K)` per paper §2.2 / §4.2 step 6; the resulting
/// seed feeds **all** internal randomness consumed by `WHIR.Open` — sumcheck
/// round-polynomial masks (Construction 6.3), code-switching mask oracles
/// (Construction 9.7), and OOD answers (Lemma 9.3) — so that `Sign` is a
/// pure function of `(sk, pk, msg)` rather than a sampler. The trailing two
/// `0x20` bytes pad the seven-character `JV-OPRD` ASCII tag out to the
/// fixed 8-byte slot.
pub(crate) const JV_OPRD: [u8; 8] = *b"JV-OPRD ";
/// Domain tag for the OOD binding point derivation (XOF). Used by
/// [`crate::keygen`] to derive the out-of-domain point `z ∈ F` from
/// `mathsf{root}` and by [`crate::verify`] to re-derive the same `z`
/// (the signer never transmits it). Bound by the OOD value `w = f(z)`
/// stored in the public key, this collapses an outsider's cap-binding
/// extraction from "needs ≥ ⌈M/K⌉ accepting signatures" to "one accepting
/// signature suffices" (paper §5.1, Theorem 13). The trailing three
/// `0x20` bytes pad the five-character `JV-OOD` ASCII tag out to the
/// fixed 8-byte slot.
pub(crate) const JV_OOD: [u8; 8] = *b"JV-OOD  ";

// ---------------------------------------------------------------------------
// Hash family selector
// ---------------------------------------------------------------------------

/// Choice of hash primitive for [`hash`].
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub(crate) enum Family {
	/// Poseidon2-Goldilocks-12 sponge (arithmetic-friendly).
	Arith,
	/// SHAKE256 extendable-output function.
	Xof,
}

// ---------------------------------------------------------------------------
// Length-prefixed domain-tagged hash
// ---------------------------------------------------------------------------

/// Hash a sequence of byte-string inputs with the chosen `family`, prefixing
/// the canonical 8-byte domain `tag` and the length-prefix framing described
/// in the module docs. Returns exactly `out_len` output bytes.
pub(crate) fn hash(family: Family, tag: [u8; 8], inputs: &[&[u8]], out_len: usize) -> Vec<u8> {
	let total = 8 + inputs.iter().map(|x| 8 + x.len()).sum::<usize>();
	let mut buf = Vec::with_capacity(total);

	buf.extend_from_slice(&tag);
	for input in inputs {
		buf.extend_from_slice(&(input.len() as u64).to_le_bytes());
		buf.extend_from_slice(input);
	}

	match family {
		Family::Arith => poseidon2_hash(&buf, out_len),
		Family::Xof => shake256_hash(&buf, out_len),
	}
}

// ---------------------------------------------------------------------------
// SHAKE256
// ---------------------------------------------------------------------------

/// SHAKE256 of `input`, squeezing `out_len` bytes.
fn shake256_hash(input: &[u8], out_len: usize) -> Vec<u8> {
	let mut hasher = Shake256::default();
	hasher.update(input);
	let mut reader = hasher.finalize_xof();
	let mut out = vec![0u8; out_len];
	reader.read(&mut out);
	out
}

// ---------------------------------------------------------------------------
// Poseidon2-Goldilocks-12 sponge
// ---------------------------------------------------------------------------

/// Sponge state width (field elements).
const POSEIDON2_WIDTH: usize = 12;
/// Sponge rate in field elements (8 → 64 bytes per absorb).
const POSEIDON2_RATE: usize = 8;
/// Bytes per absorbed rate block.
const POSEIDON2_RATE_BYTES: usize = POSEIDON2_RATE * 8;

/// Sponge-hash `input` bytes with Poseidon2-Goldilocks-12, producing `out_len`
/// bytes.
///
/// Padding: `"10*"` — append `0x01` then zero-pad to a multiple of
/// `POSEIDON2_RATE_BYTES`.
fn poseidon2_hash(input: &[u8], out_len: usize) -> Vec<u8> {
	let perm = default_goldilocks_poseidon2_12();

	let padded_len = {
		let raw = input.len() + 1;
		if raw.is_multiple_of(POSEIDON2_RATE_BYTES) {
			raw
		} else {
			raw + (POSEIDON2_RATE_BYTES - raw % POSEIDON2_RATE_BYTES)
		}
	};
	let mut padded = vec![0u8; padded_len];
	padded[..input.len()].copy_from_slice(input);
	padded[input.len()] = 0x01;

	let mut state = [Goldilocks::ZERO; POSEIDON2_WIDTH];
	for chunk in padded.chunks_exact(POSEIDON2_RATE_BYTES) {
		for (i, elem_bytes) in chunk.chunks_exact(8).enumerate() {
			state[i] += bytes_to_goldilocks(elem_bytes);
		}
		perm.permute_mut(&mut state);
	}

	let mut output = Vec::with_capacity(out_len);
	loop {
		for elem in &state[..POSEIDON2_RATE] {
			output.extend_from_slice(&elem.as_canonical_u64().to_le_bytes());
			if output.len() >= out_len {
				output.truncate(out_len);
				return output;
			}
		}
		perm.permute_mut(&mut state);
	}
}

/// Deserialise 8 little-endian bytes into a canonical `Goldilocks` element.
#[inline]
fn bytes_to_goldilocks(chunk: &[u8]) -> Goldilocks {
	let raw = u64::from_le_bytes(chunk.try_into().unwrap());
	let v = if raw >= Goldilocks::ORDER_U64 {
		raw - Goldilocks::ORDER_U64
	} else {
		raw
	};
	Goldilocks::from_u64(v)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
	use super::*;

	#[test]
	fn arith_is_deterministic() {
		let a = hash(Family::Arith, JV_WHIR, &[b"hello"], 32);
		let b = hash(Family::Arith, JV_WHIR, &[b"hello"], 32);
		assert_eq!(a, b);
	}

	#[test]
	fn xof_is_deterministic() {
		let a = hash(Family::Xof, JV_POSN, &[b"hello"], 32);
		let b = hash(Family::Xof, JV_POSN, &[b"hello"], 32);
		assert_eq!(a, b);
	}

	#[test]
	fn domain_tags_separate() {
		let a = hash(Family::Xof, JV_SEED, &[b"x"], 32);
		let b = hash(Family::Xof, JV_POSN, &[b"x"], 32);
		assert_ne!(a, b);
	}

	#[test]
	fn families_separate() {
		let a = hash(Family::Arith, JV_WHIR, &[b"x"], 32);
		let b = hash(Family::Xof, JV_WHIR, &[b"x"], 32);
		assert_ne!(a, b);
	}

	#[test]
	fn variable_output_length() {
		let a = hash(Family::Xof, JV_POSN, &[b"x"], 16);
		let b = hash(Family::Xof, JV_POSN, &[b"x"], 64);
		assert_eq!(a.len(), 16);
		assert_eq!(b.len(), 64);
		assert_eq!(&a[..], &b[..16]);
	}

	#[test]
	fn length_prefix_is_injective_under_concat() {
		// Without a length prefix, ("abcd", "ef") and ("abc", "def") would
		// collide. With our framing they cannot.
		let a = hash(Family::Xof, JV_POSN, &[b"abcd", b"ef"], 32);
		let b = hash(Family::Xof, JV_POSN, &[b"abc", b"def"], 32);
		assert_ne!(a, b);
	}

	#[test]
	fn spec_tags_are_present() {
		// Per paper §2.2: 8-byte ASCII strings right-padded with 0x20 (space).
		assert_eq!(&JV_SEED, b"JV-SEED ");
		assert_eq!(&JV_RZK, b"JV-RZK  ");
		assert_eq!(&JV_POSN, b"JV-POSN ");
		assert_eq!(&JV_FSCH, b"JV-FSCH ");
		assert_eq!(&JV_WHIR, b"JV-WHIR ");
		assert_eq!(&JV_OPEN, b"JV-OPEN ");
		assert_eq!(&JV_OPRD, b"JV-OPRD ");
		assert_eq!(&JV_OOD, b"JV-OOD  ");
		// Pairwise distinct as 8-byte strings.
		let tags: [&[u8; 8]; 8] = [
			&JV_SEED, &JV_RZK, &JV_POSN, &JV_FSCH, &JV_WHIR, &JV_OPEN, &JV_OPRD, &JV_OOD,
		];
		for i in 0..tags.len() {
			for j in (i + 1)..tags.len() {
				assert_ne!(tags[i], tags[j]);
			}
		}
	}
}
