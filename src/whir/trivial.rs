//! Final "trivial" step that closes the WHIR recursion.
//!
//! After enough codeswitch rounds the inner message is short enough to be
//! sent in full and checked directly. The prover writes the folded message
//! out element-by-element and opens the VC at a few random positions; the
//! verifier re-encodes that message, checks the openings against the
//! re-encoded codeword, and checks the linear-form dot product.

use spongefish::ProverState;

use super::commitment::CodeCommitmentProverHandle;
use super::transcript_io::sample_positions_prover;
use super::vc::Opening;

/// The trivial step's parameter — the number of random positions to open.
pub(crate) struct Trivial {
	pub(crate) queries: usize,
}

impl Trivial {
	/// Run the prover side of the trivial step:
	///   1. write the folded message into the transcript element-by-element,
	///   2. sample `queries` random codeword positions,
	///   3. open the VC at those positions.
	pub(crate) fn prove<CCH>(&self, transcript: &mut ProverState, input: CCH) -> Opening<CCH::VC>
	where
		CCH: CodeCommitmentProverHandle,
	{
		// Per-element absorbs so the verifier's per-element reads (and thus
		// its FS challenge derivation) stay byte-aligned with the prover's.
		for element in input.msg() {
			transcript.prover_message(element);
		}

		let positions = sample_positions_prover(transcript, self.queries, input.codeword_len());
		input.open(&positions)
	}
}
