//! Demonstrates the cliff property: after `n_cliff` signatures an outsider
//! has enough `(x, y)` pairs to recover the secret polynomial `f` via
//! Lagrange interpolation, at which point the secret key becomes publicly
//! recoverable.
//!
//! Run with `cargo run --release --example cliff -- 3` (defaults to
//! `n_star = 3`). Recall that `Params::new` only accepts `n_star ∈ {1, 3, 7,
//! 15, 31, 63, 127, 255, 511, 1023, …}` (i.e. `n_star + 1` a power of two);
//! anything else will panic at construction.
//!
//! This example re-implements the position-derivation procedure inline,
//! mirroring the spec — it does not reach into Jevil's private modules.

use std::collections::HashMap;
use std::env;

use jevil::{Goldilocks4, Params, keygen, sign, verify};
use rand::rngs::OsRng;
use shake::{ExtendableOutput, Shake256, Update, XofReader};

fn main() {
	let n_star: u32 = env::args().nth(1).and_then(|s| s.parse().ok()).unwrap_or(3);
	if n_star == 0 || ((n_star + 1) & n_star) != 0 {
		eprintln!(
			"error: n_star + 1 must be a power of two; got n_star = {n_star}.\n\
			 valid values: 1, 3, 7, 15, 31, 63, 127, 255, 511, 1023, …"
		);
		std::process::exit(2);
	}

	let params = Params::new(n_star);
	println!("Cliff demo · n* = {n_star}, K = {}", Params::K);
	println!(
		"  M = {}, N = {}, T = {}, n_cliff (worst case) = {}",
		params.m(),
		params.n(),
		params.t(),
		params.n_cliff()
	);

	// Fresh entropy from the OS each run, so the recovered polynomial differs
	// every time you run this demo.
	let mut rng = OsRng;
	let (pk, sk, cache) = keygen(&mut rng, params);
	println!(
		"Keygen done. Root prefix: {:02x}{:02x}…{:02x}{:02x}",
		pk.root[0], pk.root[1], pk.root[30], pk.root[31]
	);

	let t = params.t();
	let m = params.m();
	let mut pairs: HashMap<[u8; 32], Goldilocks4> = HashMap::new();
	let mut signed = 0usize;
	while pairs.len() < m {
		let msg = format!("msg-{signed}");
		let sig = sign(&sk, &pk, &cache, params, msg.as_bytes());
		verify(&pk, params, msg.as_bytes(), &sig).expect("honest sig verifies");

		let positions = derive_positions(&pk.root, msg.as_bytes(), Params::K as usize, t);
		for (pos, y) in positions.iter().zip(&sig.y_values) {
			let x = psi(*pos as u64, t as u64);
			pairs.insert(x.to_bytes(), *y);
		}
		signed += 1;
		println!(
			"  signed msg {signed:>4}: outsider has {:>5} distinct (x, y) pairs (M = {m})",
			pairs.len()
		);
		if signed > 4 * params.n_cliff() {
			println!("\nAborting after too many signatures.");
			return;
		}
	}

	let pts: Vec<(Goldilocks4, Goldilocks4)> = pairs
		.into_iter()
		.take(m)
		.map(|(xb, y)| (Goldilocks4::from_bytes(&xb).unwrap(), y))
		.collect();
	let recovered = lagrange_interpolate(&pts);

	println!(
		"\nLagrange-interpolated degree-{} polynomial — the signing key is now publicly recoverable.",
		recovered.len() - 1
	);
	println!("(c_0 = {:?})", recovered[0]);
}

// -----------------------------------------------------------------------------
// Spec-side reference implementations (don't reach into Jevil internals).
// -----------------------------------------------------------------------------

fn lagrange_interpolate(pairs: &[(Goldilocks4, Goldilocks4)]) -> Vec<Goldilocks4> {
	let n = pairs.len();
	let mut coeffs = vec![Goldilocks4::ZERO; n];
	for i in 0..n {
		let mut num = vec![Goldilocks4::ONE];
		let mut den = Goldilocks4::ONE;
		for j in 0..n {
			if i == j {
				continue;
			}
			let mut next = vec![Goldilocks4::ZERO; num.len() + 1];
			for (k, c) in num.iter().enumerate() {
				next[k] -= *c * pairs[j].0;
				next[k + 1] += *c;
			}
			num = next;
			den *= pairs[i].0 - pairs[j].0;
		}
		let den_inv = den.try_inverse().expect("distinct points");
		for (k, c) in num.iter().enumerate() {
			coeffs[k] += pairs[i].1 * *c * den_inv;
		}
	}
	coeffs
}

const JV_POSN: [u8; 8] = *b"JV-POSN\0";

fn shake256_tagged(tag: [u8; 8], inputs: &[&[u8]], out_len: usize) -> Vec<u8> {
	let mut hasher = Shake256::default();
	hasher.update(&tag);
	for input in inputs {
		hasher.update(&(input.len() as u64).to_le_bytes());
		hasher.update(input);
	}
	let mut reader = hasher.finalize_xof();
	let mut out = vec![0u8; out_len];
	reader.read(&mut out);
	out
}

fn psi(i: u64, t: u64) -> Goldilocks4 {
	let log_t = t.trailing_zeros() as usize;
	Goldilocks4::two_adic_generator(log_t).pow(i)
}

fn derive_positions(root: &[u8; 32], msg: &[u8], k: usize, t: usize) -> Vec<usize> {
	let log_t = t.trailing_zeros() as usize;
	let b = log_t.div_ceil(8);
	let initial_bytes = 32 + k * b * 4;
	let mut stream = shake256_tagged(JV_POSN, &[root, msg], initial_bytes);
	let mut cursor = 0usize;
	let mut refill = 0u64;
	let mut pool: Vec<usize> = (0..t).collect();
	let mut indices = Vec::with_capacity(k);
	let mut pool_size = t;
	for _ in 0..k {
		let m = pool_size as u64;
		let cutoff = ((1u128 << (8 * b)) / m as u128) * m as u128;
		let mask: u128 = if 8 * b == 128 {
			u128::MAX
		} else {
			(1u128 << (8 * b)) - 1
		};
		let j = loop {
			if cursor + b > stream.len() {
				refill += 1;
				let tag = refill.to_le_bytes();
				stream = shake256_tagged(JV_POSN, &[root, msg, &tag], initial_bytes);
				cursor = 0;
			}
			let mut buf = [0u8; 16];
			buf[..b].copy_from_slice(&stream[cursor..cursor + b]);
			cursor += b;
			let raw = u128::from_le_bytes(buf) & mask;
			if raw < cutoff {
				break (raw % m as u128) as usize;
			}
		};
		indices.push(pool[j]);
		let last = pool_size - 1;
		pool.swap(j, last);
		pool_size -= 1;
	}
	indices.sort_unstable();
	indices
}
