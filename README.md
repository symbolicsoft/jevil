# Jevil

> [!CAUTION]
> ## ⚠️ EXPERIMENTAL — DO NOT USE IN PRODUCTION ⚠️
>
> **This is a research-grade proof-of-concept implementation of a brand-new,
> completely novel cryptographic scheme.** Both the *scheme itself* and this
> *implementation* have received **close to zero peer review**.
>
> - The construction has **not** been vetted by the cryptographic community.
> - The security proofs have **not** been independently verified.
> - The code has **not** been audited.
> - There are almost certainly bugs, side channels, and possibly fundamental
>   design flaws that have not yet been discovered.
> - APIs, wire formats, and parameter choices may change without notice.
>
> Treat this repository as a **research artifact only**. Do not use it to
> protect anything you care about. Do not deploy it. Do not rely on it for
> any security property whatsoever.

> A stateless few-time signature scheme with a sharp cliff at the
> `(n* + 1)`-th signature.

Jevil is a post-quantum few-time signature scheme parameterised by a single
signing budget `n*`. Signatures `1..=n*` are existentially unforgeable; at
the `(n* + 1)`-th signature the secret signing key becomes **publicly
recoverable** by anyone observing the signatures — the cap is enforced not by
counters or hardware, but by the algebraic structure of a single committed
polynomial. `Params::new` accepts only `n_star` values for which `n_star + 1`
is a power of two (the paper's recommended regime), so `n_cliff = n_star + 1`
exactly for every deployment.

| | |
| --- | --- |
| **Public key** | ~36 bytes |
| **Secret key** | 32 bytes |
| **Signature** | ~30–50 KB |
| **Classical security** | ≥ 128 bits below the cliff |
| **Quantum security** | ≥ 85 bits at default capacity, raise to 128 bits with `c = 6` |

## When to use Jevil

Jevil is designed for **audit-budgeted credentials** — settings where
over-signing must be *self-exposing* rather than merely policy-forbidden:

- a firmware vendor capping its own release count,
- an operator binding themselves to a per-tenure attestation budget,
- an ephemeral session signer with a per-session cap,
- any audit-budgeted credential whose holder shouldn't be trusted to honour
  the budget unilaterally.

It is **not** a general-purpose signature scheme. For everyday signing use a
stateful or unlimited-use post-quantum scheme such as ML-DSA or Falcon.
Jevil's value is in the cliff.

## Installation

Add to your `Cargo.toml`:

```toml
[dependencies]
jevil = "0.1"
```

The crate is `#![forbid(unsafe_code)]` and exposes a single library target.

## Quick start

```rust
use jevil::{Params, keygen, sign, verify};
use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;

// Pick a signing budget. n_star = 7 means: up to 7 honest signatures;
// the cliff fires at the 8th. Params::new accepts only n_star values for
// which n_star + 1 is a power of two (1, 3, 7, 15, 31, …).
let params = Params::new(7);

// Generate a fresh key.
let mut rng = ChaCha20Rng::seed_from_u64(0);
let (pk, sk, cache) = keygen(&mut rng, params);

// Sign a message.
let signature = sign(&sk, &pk, &cache, params, b"firmware-image-v1.0.0");

// Anyone holding `pk` can verify.
assert!(verify(&pk, params, b"firmware-image-v1.0.0", &signature).is_ok());
```

Try the bundled examples:

```bash
cargo run --release --example basic    # minimal sign/verify
cargo run --release --example bench    # latencies across n_star
cargo run --release --example cliff -- 3  # public-key recovery demo
```

## How the cliff works

The secret is a univariate polynomial `f ∈ F[X]` of degree `D = M − 1` over
the quartic Goldilocks extension `F_{q₀^4}`, derived deterministically from
the 32-byte secret seed. The public key is a [WHIR][whir] commitment to the
length-`M` coefficient vector `c = (c_0, …, c_{M−1})`.

A signature on message `M` opens `f` at `K = 16` message-derived positions
`x_1, …, x_K`, revealing `(y_t)_{t=1..K} = (f(x_t))_{t=1..K}` together with a
single batched WHIR linear-form proof. After `n` honest signatures, an
observer holds at most `nK` distinct `(x, f(x))` pairs:

- **`n ≤ n*`**: at most `n*·K ≤ M − K ≤ D` evaluations. `f` is information-
  theoretically undetermined; the observer's prediction at any unseen point
  is correct with probability `1 / |F| ≈ 2^{−256}`.
- **`n ≥ n_cliff = ⌈M/K⌉`** (worst case): at least `M = D + 1` distinct
  evaluations have accumulated. Lagrange interpolation recovers `f` in
  `O(D²)` field operations, and the secret becomes publicly known.

The cliff is *intrinsic to the public key*. An adversarial signer cannot
extend it by committing to a non-codeword or by reusing proofs across
roots — see the paper for the precise cap-binding lemma.

The WHIR commitment is instantiated as **zk-WHIR** (honest-verifier
zero-knowledge WHIR): per-signature transcripts are perfectly simulable
from `(pk, msg, (y_t))` alone. Below the cliff, no observer holding any
number of signatures learns anything about `f` beyond the `K` explicitly
opened `y_t` values per signature. See [`docs/zkwhir-spec-compliance.md`](docs/zkwhir-spec-compliance.md)
for the per-construction status against the [zk-WHIR paper][zkwhir].

## Parameter selection

Jevil takes a single deployment-time integer `n_star`, which must satisfy
**`n_star + 1` is a power of two**:

```text
n_star ∈ {1, 3, 7, 15, 31, 63, 127, 255, 511, 1023, 2047, 4095, 8191, 16_383}
```

`Params::new` panics on any other value. This restriction pins the paper's
*recommended regime*, in which `(n_star + 1) · K` is itself a power of two,
`M = (n_star + 1) · K` exactly, and the cliff fires precisely at signature
`n_star + 1`. Outside this regime `M` would round up to the next power of
two, leaving a gap between `n_star` and `n_cliff` that erodes the HORS
coverage margin (§5.3 of the paper); rejecting bad values at construction
removes that footgun.

From `n_star` every other quantity is derived:

- `M = (n_star + 1) · K` — the cliff dimension; degree bound is `D = M − 1`.
- `T = nextpow2(n_star · K · 2^{128/K})` — the HORS-coverage position space.
- `n_cliff = M / K = n_star + 1` — the first signature index at which the
  worst-case outsider recovers `f`.

Reference sizes at `K = 16`:

| `n_star` | `M`    | `T`    | KeyGen | Sig    |
|---------:|-------:|-------:|-------:|-------:|
| 127      | 2¹¹    | 2¹⁹    | 0.1 s  | 30 KB  |
| 1023     | 2¹⁴    | 2²²    | 1 s    | 35 KB  |
| 16,383   | 2¹⁸    | 2²⁶    | 30 s   | 45 KB  |

## API overview

| Item | Purpose |
| --- | --- |
| [`Params`][p] | Configuration (a single `n_star` field). |
| [`keygen`][k] | Generate `(PublicKey, SecretKey, SignerCache)` from an RNG. |
| [`sign`][s] | Produce a `Signature` for `(sk, pk, msg)`. |
| [`verify`][v] | Check a signature against `pk` and `msg`. |
| [`PublicKey`][pk] | 36-byte serialised public key. |
| [`SecretKey`][sk] | 32-byte signing seed. |
| [`SignerCache`][c] | Pre-derived signer state, for fast repeated signing. |
| [`Signature`][sg] | `K` y-values plus the inline WHIR proof. |
| [`Error`][e] | Verification / parsing failure modes. |
| [`Goldilocks4`][g] | Re-exported for inspecting `Signature::y_values`. |

[p]: https://docs.rs/jevil/latest/jevil/struct.Params.html
[k]: https://docs.rs/jevil/latest/jevil/fn.keygen.html
[s]: https://docs.rs/jevil/latest/jevil/fn.sign.html
[v]: https://docs.rs/jevil/latest/jevil/fn.verify.html
[pk]: https://docs.rs/jevil/latest/jevil/struct.PublicKey.html
[sk]: https://docs.rs/jevil/latest/jevil/struct.SecretKey.html
[c]: https://docs.rs/jevil/latest/jevil/struct.SignerCache.html
[sg]: https://docs.rs/jevil/latest/jevil/struct.Signature.html
[e]: https://docs.rs/jevil/latest/jevil/enum.Error.html
[g]: https://docs.rs/jevil/latest/jevil/struct.Goldilocks4.html

The public API is intentionally small. WHIR plumbing, hashing, field
arithmetic, the merkle tree, and the lift are all crate-private — only the
items above are exposed.

## Security

### Below the cliff

A PPT forger producing a fresh signature must either:

1. supply true `y_t` values at unseen `x_t`. Below the cliff this requires
   guessing `f` on undetermined points — success probability `|F|^{−K} ≈
   2^{−256K}` per attempt.
2. supply *false* `y_t` values. WHIR's knowledge soundness then implies a
   break of either the Reed–Solomon proximity-gap conjecture or the
   Poseidon2 collision resistance — each at ~2⁻¹²⁸ classical.
3. find a fresh message whose `K` hash-derived positions all land in the
   already-seen pool (the **HORS coverage path**). Bounded by `(nK/T)^K`,
   capped at 2⁻¹²⁸ by the choice of `T`.

The three paths together upper-bound classical forgery probability at ~2⁻¹²⁸.

### At and above the cliff

After `n_cliff = ⌈M/K⌉` honest signatures with all-distinct positions, the
outsider has ≥ `D + 1` distinct `(x, f(x))` pairs and recovers `f` in
`O(D²)` operations. From that point the outsider matches the legitimate
signer bit-for-bit.

### Hiding

Each signature reveals exactly `K` `(x_t, y_t)` pairs and one WHIR proof.
Below the cliff, the WHIR proof is honest-verifier zero-knowledge with
respect to a query-bounded distinguisher: the joint distribution of revealed
codeword symbols and sumcheck round polynomials is perfectly simulable from
`(pk, msg, (y_t))` alone.

### Post-quantum security

All primitives are post-quantum:

- Poseidon2-Goldilocks, capacity `c = 4` → ~2²⁵⁶ classical /
  ~2¹²⁸ quantum second-preimage; ~2¹²⁸ classical / ~2⁸⁵ quantum collision.
- SHAKE256 → same numbers at its 256-bit output.
- Reed–Solomon proximity-gap conjecture backs WHIR's soundness.

For ≥ 128-bit *quantum* security throughout, raise the Poseidon2 capacity to
`c = 6` (384-bit state). Classical security is ≥ 128 bits in all paths at
the default `c = 4`.

### Multi-target

With `N_keys` honest signers in the deployment, replace `2^{128/K}` with
`2^{(128 + log₂ N_keys)/K}` in the `T` formula above to preserve 128-bit
*population-level* security against any-key forgery.

## Testing

Unit and integration tests cover:

- Field arithmetic correctness (commutativity, distributivity, inverse, NTT).
- Hash domain separation (every tag combination is distinct).
- Position-derivation distinctness, sortedness, and rejection sampling bias.
- Lift / symbolic-α correctness across `(ν, ν', K, R)` sweeps.
- Signature round-trip across `n_star ∈ {1, 3, 7, 15, 31}`.
- Tamper rejection: y-value flip, proof byte flip, wrong root, wrong
  message, wrong `n_star`, non-canonical field element, truncated signature.
- Determinism (same seed → byte-identical pk / signature).
- A pinned known-answer test (KAT) for `n_star = 3`, `seed = 0`,
  `msg = "jevil-kat-fixture"`.
- The **cliff property**: at `n_cliff` signatures, Lagrange interpolation
  recovers `f` byte-for-byte from observed `(x, y)` pairs.

```bash
cargo test                                            # standard
cargo test --release --test slow -- --ignored         # n_star = 127, 1023
KAT_UPDATE=1 cargo test --test kat -- --nocapture     # regenerate fixtures
```

## License

Licensed under either of [Apache License, Version 2.0](LICENSE-APACHE) or
[MIT license](LICENSE-MIT) at your option.

[whir]: https://eprint.iacr.org/2024/1586
[zkwhir]: https://eprint.iacr.org/2026/391
