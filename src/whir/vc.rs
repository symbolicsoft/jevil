//! Vector commitment used inside WHIR.
//!
//! WHIR is parameterised over an abstract [`VectorCommitment`] trait, but
//! Jevil only ever instantiates it with [`MerkleVc`] — a Poseidon2-Goldilocks
//! Merkle tree where each *leaf* is the hash of a `Vec<Goldilocks4>`
//! (corresponding to one position of an interleaved codeword) and each
//! *internal node* is the hash of its two children. See [`crate::merkle`] for
//! the tree primitive itself.

use spongefish::{Encoding, NargDeserialize};

use crate::field::Goldilocks4;
use crate::hash::{Family, JV_WHIR, hash};
use crate::merkle::{MerkleTree, verify_path};

/// An opened-position bundle for a vector commitment.
#[derive(spongefish::Encoding, spongefish::NargDeserialize)]
pub(crate) struct Opening<VC: VectorCommitment> {
	/// The opened symbols at the queried positions, in order.
	pub openings: Vec<VC::Alphabet>,
	/// Whatever proof material the underlying VC needs (e.g. Merkle paths).
	pub vc_proof: VC::OpeningProof,
}

/// Abstract interface to a position-binding vector commitment.
pub(crate) trait VectorCommitment: Sized {
	/// Symbol type of the committed vector.
	type Alphabet: Clone + Encoding;
	/// Commitment digest type.
	type Commitment: Encoding + NargDeserialize;
	/// Proof type that accompanies an opening.
	type OpeningProof: Encoding;
	/// Prover-side state retained between [`Self::commit`] and [`Self::open`].
	type CommitState;

	/// Commit to `input`. Returns the commitment digest and prover state.
	fn commit(&self, input: &[Self::Alphabet]) -> (Self::Commitment, Self::CommitState);

	/// Open the commitment at the given positions, returning the claimed
	/// symbols together with a proof.
	fn open(&self, state: &Self::CommitState, indexes: &[usize]) -> Opening<Self>;

	/// Verify an opening against `commitment` at the given positions.
	fn verify(
		&self,
		commitment: &Self::Commitment,
		indexes: &[usize],
		proof: &Opening<Self>,
	) -> bool;
}

// ---------------------------------------------------------------------------
// MerkleVc
// ---------------------------------------------------------------------------

/// A binary Merkle vector commitment whose alphabet is `Vec<Goldilocks4>`
/// (one position of an interleaved codeword).
///
/// The leaf hash for position `i` is
///
/// ```text
/// H_arith(JV-WHIR, concat(positions[i].iter().flat_map(|g| g.to_bytes())); 32)
/// ```
///
/// where `H_arith` is the Poseidon2-Goldilocks sponge of [`crate::hash`].
pub(crate) struct MerkleVc {
	/// Number of committed positions.
	pub(crate) n: usize,
}

impl MerkleVc {
	/// Construct an `MerkleVc` for `n` positions.
	pub(crate) fn new(n: usize) -> Self {
		Self { n }
	}
}

impl VectorCommitment for MerkleVc {
	type Alphabet = Vec<Goldilocks4>;
	type Commitment = [u8; 32];
	type OpeningProof = Vec<[u8; 32]>;
	type CommitState = (MerkleTree, Vec<Vec<Goldilocks4>>);

	fn commit(&self, input: &[Self::Alphabet]) -> (Self::Commitment, Self::CommitState) {
		assert_eq!(
			input.len(),
			self.n,
			"MerkleVc::commit: input length mismatch"
		);
		let leaf_hashes: Vec<[u8; 32]> = input
			.iter()
			.map(|symbols| {
				let mut buf = Vec::with_capacity(symbols.len() * 32);
				for s in symbols {
					buf.extend_from_slice(&s.to_bytes());
				}
				let h = hash(Family::Arith, JV_WHIR, &[&buf], 32);
				h.try_into().unwrap()
			})
			.collect();
		let tree = MerkleTree::build_from_hashes(leaf_hashes);
		let root = tree.root();
		(root, (tree, input.to_vec()))
	}

	fn open(&self, state: &Self::CommitState, indexes: &[usize]) -> Opening<Self> {
		let (tree, input) = state;
		let openings: Vec<Vec<Goldilocks4>> = indexes.iter().map(|&i| input[i].clone()).collect();
		let mut vc_proof: Vec<[u8; 32]> = Vec::new();
		for &i in indexes {
			vc_proof.extend(tree.path(i));
		}
		Opening { openings, vc_proof }
	}

	fn verify(
		&self,
		commitment: &Self::Commitment,
		indexes: &[usize],
		proof: &Opening<Self>,
	) -> bool {
		if proof.openings.len() != indexes.len() {
			return false;
		}
		let path_len = self.n.next_power_of_two().trailing_zeros() as usize;
		if proof.vc_proof.len() != indexes.len() * path_len {
			return false;
		}
		for (k, &i) in indexes.iter().enumerate() {
			let mut buf = Vec::with_capacity(proof.openings[k].len() * 32);
			for s in &proof.openings[k] {
				buf.extend_from_slice(&s.to_bytes());
			}
			let leaf_hash: [u8; 32] = hash(Family::Arith, JV_WHIR, &[&buf], 32)
				.try_into()
				.unwrap();
			let path = &proof.vc_proof[k * path_len..(k + 1) * path_len];
			if !verify_path(*commitment, i, leaf_hash, path) {
				return false;
			}
		}
		true
	}
}

#[cfg(test)]
mod tests {
	use super::*;
	use crate::field::Goldilocks;

	fn dummy_symbols(i: usize, k: usize) -> Vec<Goldilocks4> {
		(0..k)
			.map(|j| {
				Goldilocks4::new([
					Goldilocks::new((i * 100 + j) as u64),
					Goldilocks::new(0),
					Goldilocks::new(0),
					Goldilocks::new(0),
				])
			})
			.collect()
	}

	#[test]
	fn round_trip() {
		let n = 64;
		let k = 4;
		let vc = MerkleVc::new(n);
		let input: Vec<Vec<Goldilocks4>> = (0..n).map(|i| dummy_symbols(i, k)).collect();
		let (commit, state) = vc.commit(&input);
		let indexes = [3, 17, 42, 50];
		let opening = vc.open(&state, &indexes);
		assert!(vc.verify(&commit, &indexes, &opening));
	}

	#[test]
	fn tampering_rejected() {
		let n = 32;
		let k = 4;
		let vc = MerkleVc::new(n);
		let input: Vec<Vec<Goldilocks4>> = (0..n).map(|i| dummy_symbols(i, k)).collect();
		let (commit, state) = vc.commit(&input);
		let indexes = [10];
		let mut opening = vc.open(&state, &indexes);
		opening.openings[0][0] += Goldilocks4::ONE;
		assert!(!vc.verify(&commit, &indexes, &opening));
	}
}
