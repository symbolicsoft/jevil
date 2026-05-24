//! Binary Merkle tree used internally by WHIR's vector commitment.
//!
//! This is a private utility — callers pass in pre-hashed leaves and get back
//! a root plus path-extraction / path-verification routines. The internal-node
//! hash uses [`crate::hash::JV_WHIR`] (Poseidon2-Goldilocks-12).
//!
//! ## Why the length-prefixed framing is safe to reuse
//!
//! The same domain tag (`JV-WHIR`) is used for both leaf hashing (called by
//! the vector-commitment layer in [`crate::whir::vc`]) and internal-node
//! hashing here. This is fine because [`crate::hash::hash`] prefixes every
//! input with its 8-byte length, so leaf inputs (variable length) cannot
//! collide with internal-node inputs (always two 32-byte sibling hashes).

use crate::hash::{Family, JV_WHIR, hash};

/// A binary Merkle tree over pre-hashed leaves.
///
/// `layers[0]` is the leaf-hash layer; `layers[k+1]` holds the parents of
/// `layers[k]`. The number of leaves is rounded up to the next power of two by
/// duplicating the last hash (no padding tag is used because the consumer
/// never accepts an opening for a padded index).
#[derive(Clone, Debug)]
pub(crate) struct MerkleTree {
	layers: Vec<Vec<[u8; 32]>>,
}

impl MerkleTree {
	/// Build a Merkle tree over the given leaf hashes. Pads up to a power of
	/// two by duplicating the last hash if needed.
	///
	/// Panics if `leaf_hashes` is empty.
	pub(crate) fn build_from_hashes(mut leaf_hashes: Vec<[u8; 32]>) -> Self {
		assert!(
			!leaf_hashes.is_empty(),
			"MerkleTree::build_from_hashes: empty"
		);
		let n = leaf_hashes.len().next_power_of_two();
		while leaf_hashes.len() < n {
			leaf_hashes.push(*leaf_hashes.last().unwrap());
		}
		let mut layers = vec![leaf_hashes];
		while layers.last().unwrap().len() > 1 {
			let prev = layers.last().unwrap();
			let next: Vec<[u8; 32]> = prev
				.chunks(2)
				.map(|pair| {
					let v = hash(Family::Arith, JV_WHIR, &[&pair[0], &pair[1]], 32);
					v.try_into().unwrap()
				})
				.collect();
			layers.push(next);
		}
		Self { layers }
	}

	/// The 32-byte root.
	pub(crate) fn root(&self) -> [u8; 32] {
		self.layers.last().unwrap()[0]
	}

	/// The Merkle authentication path (sibling hashes, leaf layer to root) for
	/// `leaf_idx`. The returned vector has length `log₂(num_leaves)`.
	pub(crate) fn path(&self, mut leaf_idx: usize) -> Vec<[u8; 32]> {
		assert!(leaf_idx < self.layers[0].len());
		let mut path = Vec::with_capacity(self.layers.len() - 1);
		for layer in &self.layers[..self.layers.len() - 1] {
			let sib = leaf_idx ^ 1;
			path.push(layer[sib]);
			leaf_idx >>= 1;
		}
		path
	}
}

/// Verify a Merkle authentication path against `root`. Returns `true` iff
/// recomputing the path from `(leaf_hash, path)` yields `root`.
pub(crate) fn verify_path(
	root: [u8; 32],
	mut leaf_idx: usize,
	leaf_hash: [u8; 32],
	path: &[[u8; 32]],
) -> bool {
	let mut cur = leaf_hash;
	for &sib in path {
		let (l, r) = if leaf_idx & 1 == 0 {
			(cur, sib)
		} else {
			(sib, cur)
		};
		cur = hash(Family::Arith, JV_WHIR, &[&l, &r], 32)
			.try_into()
			.unwrap();
		leaf_idx >>= 1;
	}
	cur == root
}

#[cfg(test)]
mod tests {
	use super::*;

	fn fake_hashes(n: usize) -> Vec<[u8; 32]> {
		(0..n)
			.map(|i| {
				let mut x = [0u8; 32];
				x[0] = i as u8;
				x
			})
			.collect()
	}

	#[test]
	fn build_and_verify_path() {
		let hashes = fake_hashes(1024);
		let tree = MerkleTree::build_from_hashes(hashes.clone());
		let root = tree.root();
		for (i, leaf) in hashes.iter().enumerate() {
			let path = tree.path(i);
			assert!(verify_path(root, i, *leaf, &path), "index {i}");
		}
	}

	#[test]
	fn tampering_rejected() {
		let hashes = fake_hashes(64);
		let tree = MerkleTree::build_from_hashes(hashes.clone());
		let root = tree.root();
		let mut path = tree.path(7);
		path[0][0] ^= 1;
		assert!(!verify_path(root, 7, hashes[7], &path));
	}

	#[test]
	fn pads_to_power_of_two() {
		let hashes = fake_hashes(13);
		let tree = MerkleTree::build_from_hashes(hashes);
		assert_eq!(tree.layers.len(), 5);
	}
}
