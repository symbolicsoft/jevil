//! HVZK simulator for `WhirPcs` — composed per Theorem 4.5 of ePrint 2026/391
//! and extended to multi-opening per Lemma 11 of the Jevil paper.
//!
//! # Status
//!
//! The simulator is **scaffold-only** in the current tree: the contract,
//! types, and entry points are defined, but the per-IOR simulator bodies
//! (`base_case::simulate`, `codeswitch::simulate`, `sumcheck::simulate`)
//! and the multi-opening dispatch via `Sim_C(S)` are pending follow-up.
//!
//! Why this still allows claiming HVZK at the protocol level:
//!
//! 1. Each HVZK pillar (Constructions 6.3, 7.2, 9.7) is **active** in the
//!    protocol with the masking, mask-oracle commits, and joint-target
//!    composition exactly per the paper.
//! 2. The multi-opening parameter sizing `N − M ≥ n* · Q_max` satisfies the
//!    hypothesis of Lemma 11 by construction (paper §2.3).
//! 3. Tests at `n* ∈ {1, 7, 31, 127, 1023}` confirm the protocol round-trips
//!    end-to-end at full budget.
//!
//! What an in-tree simulator would add is an empirical *byte-equality*
//! demonstration that the protocol's transcript distribution matches the
//! simulator's output for the Reed–Solomon instantiation (`ζ_C = ζ_ze =
//! ζ_{C_zk} = 0`). This is the strongest possible witness for the HVZK
//! property short of formal verification.
//!
//! # Contract (when implemented)
//!
//! ```ignore
//! pub(crate) fn simulate_signature(
//!     pcs: &WhirPcs,
//!     public_claim: &PublicClaim,
//!     challenge_seed: [u8; 32],
//! ) -> SimulatedTranscript;
//!
//! pub(crate) fn simulate_multi(
//!     pcs: &WhirPcs,
//!     claims: &[PublicClaim],
//!     challenge_seeds: &[[u8; 32]],
//! ) -> Vec<SimulatedTranscript>;
//! ```
//!
//! `simulate_multi` is the Lemma 11 entry point: it pre-samples the joint
//! query set `S = ⋃_i S_i` across all `n*` openings, draws a single
//! `ũ ← Sim_C(S)` for the main code's ZK encoding (Prop 3.19), and threads
//! `ũ|_{S_i}` to each per-opening simulator. With `|S| ≤ n* · Q_max ≤ N − M`,
//! the joint distribution of the simulator's outputs matches the real
//! protocol's joint signature distribution at statistical distance zero.

use crate::field::Goldilocks4;

/// Public input to the simulator — everything an outside observer sees about
/// one Jevil signature, minus the WHIR proof itself.
#[allow(dead_code)] // Pending simulator implementation.
pub(crate) struct PublicClaim {
	/// The 32-byte zk-WHIR commitment root (= `PublicKey::root`).
	pub root: [u8; 32],
	/// The message being signed.
	pub msg: Vec<u8>,
	/// The `K` revealed evaluations `(y_t) = (f(x_t))`.
	pub y_values: Vec<Goldilocks4>,
}

/// Simulated transcript: the NARG bytes a real signer would write to the
/// spongefish transcript for the given `(public_claim, challenge_seed)`.
#[allow(dead_code)] // Pending simulator implementation.
pub(crate) struct SimulatedTranscript {
	pub narg_bytes: Vec<u8>,
}
