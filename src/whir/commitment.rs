//! Code-based polynomial commitments.
//!
//! `CodeCommitment` packages an [`AdditiveCode`] together with a
//! [`VectorCommitment`] into a single object that can `commit` to a message
//! (encode it as a codeword, then VC-commit) and later `open` arbitrary
//! positions.
//!
//! Two flavours of prover state exist:
//!
//! - [`CodeCommitmentProverState`]: the initial state after `commit`. Holds
//!   the original message symbols.
//! - [`FoldedCodeCommitmentProverState`]: the state after one or more
//!   sumcheck rounds have folded the underlying message. Holds the folded
//!   message; the original VC state is still inside the wrapped inner state.
//!
//! Verifier-side mirrors live in this file too:
//! [`ExplicitCodeCommitmentHandle`] and [`FoldedCodeCommitmentHandle`].

use std::sync::Arc;

use spongefish::{ProverState, VerificationError, VerificationResult};

use super::code::{AdditiveCode, InterleavedCode, LinearCode};
use super::linear_form::fold_evaluations;
use super::vc::{Opening, VectorCommitment};

// ---------------------------------------------------------------------------
// CodeCommitment (prover-side factory)
// ---------------------------------------------------------------------------

/// A reusable commit-key bundling a code and a vector commitment.
pub(crate) struct CodeCommitment<EC, VC> {
	pub(crate) code: Arc<EC>,
	pub(crate) vc: Arc<VC>,
}

impl<EC, VC> CodeCommitment<EC, VC> {
	/// Pair a code with a vector commitment.
	pub(crate) fn new(code: Arc<EC>, vc: Arc<VC>) -> Self {
		Self { code, vc }
	}
}

impl<EC, VC> CodeCommitment<EC, VC>
where
	EC: AdditiveCode,
	VC: VectorCommitment<Alphabet = EC::OutputAlphabet>,
{
	/// Commit to `msg` via `code.encode` then `vc.commit`, writing the
	/// resulting commitment digest into the Fiat–Shamir transcript.
	pub(crate) fn commit(
		&self,
		transcript: &mut ProverState,
		msg: Vec<EC::InputAlphabet>,
	) -> CodeCommitmentProverState<EC, VC> {
		assert_eq!(msg.len(), self.code.msg_len());

		let encoding = self.code.encode(&msg);
		assert_eq!(encoding.len(), self.code.codeword_len());

		let (commitment, vc_state) = self.vc.commit(&encoding);
		transcript.prover_message(&commitment);

		CodeCommitmentProverState {
			code: self.code.clone(),
			vc: self.vc.clone(),
			msg,
			vc_state,
		}
	}

	/// Commit to `msg` *without* writing anything to a transcript. Returns
	/// `(root, state)`. Used by [`crate::keygen`] to derive the WHIR root
	/// for the public key independently of any signature transcript.
	pub(crate) fn commit_only(
		&self,
		msg: Vec<EC::InputAlphabet>,
	) -> (VC::Commitment, CodeCommitmentProverState<EC, VC>) {
		assert_eq!(msg.len(), self.code.msg_len());
		let encoding = self.code.encode(&msg);
		assert_eq!(encoding.len(), self.code.codeword_len());
		let (commitment, vc_state) = self.vc.commit(&encoding);
		let state = CodeCommitmentProverState {
			code: self.code.clone(),
			vc: self.vc.clone(),
			msg,
			vc_state,
		};
		(commitment, state)
	}
}

// ---------------------------------------------------------------------------
// Prover-side state (initial + folded)
// ---------------------------------------------------------------------------

/// Prover-side state immediately after `commit`.
pub(crate) struct CodeCommitmentProverState<EC, VC>
where
	EC: AdditiveCode,
	VC: VectorCommitment,
{
	pub(crate) code: Arc<EC>,
	pub(crate) vc: Arc<VC>,
	pub(crate) msg: Vec<EC::InputAlphabet>,
	pub(crate) vc_state: VC::CommitState,
}

/// Prover-side state after sumcheck folding has reduced the inner message.
pub(crate) struct FoldedCodeCommitmentProverState<EC, VC>
where
	EC: LinearCode,
	VC: VectorCommitment<Alphabet = Vec<EC::Alphabet>>,
{
	pub(crate) inner: CodeCommitmentProverState<InterleavedCode<EC>, VC>,
	pub(crate) msg: Vec<EC::Alphabet>,
}

/// Common interface to both initial and folded prover states.
pub(crate) trait CodeCommitmentProverHandle {
	type Code: AdditiveCode;
	type VC: VectorCommitment;

	/// The (logical) code at this point in the protocol.
	fn code(&self) -> &Self::Code;

	/// The current (possibly folded) message.
	fn msg(&self) -> &[<Self::Code as AdditiveCode>::InputAlphabet];

	/// Length of the codeword. Defaults to delegating to the code.
	fn codeword_len(&self) -> usize {
		self.code().codeword_len()
	}

	/// Open the underlying VC at the given positions.
	fn open(&self, indexes: &[usize]) -> Opening<Self::VC>;
}

impl<EC, VC> CodeCommitmentProverHandle for CodeCommitmentProverState<EC, VC>
where
	EC: AdditiveCode,
	VC: VectorCommitment<Alphabet = EC::OutputAlphabet>,
{
	type Code = EC;
	type VC = VC;

	fn code(&self) -> &Self::Code {
		&self.code
	}

	fn msg(&self) -> &[<Self::Code as AdditiveCode>::InputAlphabet] {
		&self.msg
	}

	fn open(&self, indexes: &[usize]) -> Opening<Self::VC> {
		self.vc.open(&self.vc_state, indexes)
	}
}

impl<EC, VC> CodeCommitmentProverHandle for FoldedCodeCommitmentProverState<EC, VC>
where
	EC: LinearCode,
	VC: VectorCommitment<Alphabet = Vec<EC::Alphabet>>,
{
	type Code = EC;
	type VC = VC;

	fn code(&self) -> &Self::Code {
		self.inner.code.inner_code()
	}

	fn msg(&self) -> &[<Self::Code as AdditiveCode>::InputAlphabet] {
		&self.msg
	}

	fn open(&self, indexes: &[usize]) -> Opening<Self::VC> {
		self.inner.vc.open(&self.inner.vc_state, indexes)
	}
}

// ---------------------------------------------------------------------------
// Verifier-side handles
// ---------------------------------------------------------------------------

/// Common interface to verifier-side commitment handles.
pub(crate) trait CodeCommitmentHandle {
	type Code: AdditiveCode;
	type VC: VectorCommitment;

	/// The logical code at this point in the protocol.
	fn code(&self) -> &Self::Code;

	/// Message length on the code side.
	fn msg_len(&self) -> usize {
		self.code().msg_len()
	}

	/// Codeword length on the code side.
	fn codeword_len(&self) -> usize {
		self.code().codeword_len()
	}

	/// Verify the supplied opening at the given positions and return the
	/// opened codeword symbols.
	fn verify_openings(
		&self,
		indexes: &[usize],
		proof: &Opening<Self::VC>,
	) -> VerificationResult<Vec<<Self::Code as AdditiveCode>::OutputAlphabet>>;
}

/// Verifier handle for the initial (un-folded) commitment.
pub(crate) struct ExplicitCodeCommitmentHandle<EC, VC>
where
	EC: AdditiveCode,
	VC: VectorCommitment<Alphabet = EC::OutputAlphabet>,
{
	pub(crate) code: Arc<EC>,
	pub(crate) vc: Arc<VC>,
	pub(crate) commitment: VC::Commitment,
}

impl<EC, VC> ExplicitCodeCommitmentHandle<EC, VC>
where
	EC: AdditiveCode,
	VC: VectorCommitment<Alphabet = EC::OutputAlphabet>,
{
	pub(crate) fn new(code: Arc<EC>, vc: Arc<VC>, commitment: VC::Commitment) -> Self {
		Self {
			code,
			vc,
			commitment,
		}
	}
}

impl<EC, VC> CodeCommitmentHandle for ExplicitCodeCommitmentHandle<EC, VC>
where
	EC: AdditiveCode,
	VC: VectorCommitment<Alphabet = EC::OutputAlphabet>,
{
	type Code = EC;
	type VC = VC;

	fn code(&self) -> &Self::Code {
		&self.code
	}

	fn verify_openings(
		&self,
		indexes: &[usize],
		proof: &Opening<Self::VC>,
	) -> VerificationResult<Vec<<Self::Code as AdditiveCode>::OutputAlphabet>> {
		(proof.openings.len() == indexes.len() && self.vc.verify(&self.commitment, indexes, proof))
			.then(|| proof.openings.clone())
			.ok_or(VerificationError)
	}
}

/// Verifier handle for a sumcheck-folded commitment.
pub(crate) struct FoldedCodeCommitmentHandle<EC, VC>
where
	EC: LinearCode,
	VC: VectorCommitment<Alphabet = Vec<EC::Alphabet>>,
{
	pub(crate) inner: ExplicitCodeCommitmentHandle<InterleavedCode<EC>, VC>,
	pub(crate) rand: Vec<EC::Alphabet>,
}

impl<EC, VC> CodeCommitmentHandle for FoldedCodeCommitmentHandle<EC, VC>
where
	EC: LinearCode,
	VC: VectorCommitment<Alphabet = Vec<EC::Alphabet>>,
{
	type Code = EC;
	type VC = VC;

	fn code(&self) -> &Self::Code {
		self.inner.code.inner_code()
	}

	fn verify_openings(
		&self,
		indexes: &[usize],
		proof: &Opening<Self::VC>,
	) -> VerificationResult<Vec<<Self::Code as AdditiveCode>::OutputAlphabet>> {
		let n = self.inner.code.interleaving_factor();
		if n == 0 || !n.is_power_of_two() || self.rand.len() != n.ilog2() as usize {
			return Err(VerificationError);
		}

		self.inner
			.verify_openings(indexes, proof)?
			.into_iter()
			.map(|opening| {
				(opening.len() == n)
					.then(|| fold_evaluations(opening, &self.rand).pop().unwrap())
					.ok_or(VerificationError)
			})
			.collect()
	}
}
