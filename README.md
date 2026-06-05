# ternary-noether

> **Noether's theorem for discrete ternary symmetries.**

## What problem does this solve?

Emmy Noether proved that every continuous symmetry of a Lagrangian implies a conservation law. But what happens when the configuration space is not a manifold—it is the **finite ternary alphabet** $\{-1, 0, +1\}$? This crate implements the discrete analogue: it derives conserved quantities (energy, momentum, angular momentum) from the discrete symmetries (translation, rotation, reflection) that leave a ternary dynamical system invariant. It is the group-theoretic backbone of any ternary physics simulation.

## Mathematical foundations

### Noether's theorem (continuous)

For a Lagrangian $L(q, \dot{q})$ invariant under a one-parameter family of transformations $q \to q + \epsilon \, \delta q$, the quantity

$$Q = \sum_i p_i \, \delta q_i$$

is conserved along trajectories. In the discrete setting the group is finite, so "conservation" means the quantity is **constant across state sequences** that respect the symmetry.

### Discrete symmetries on {-1, 0, +1}

| Symmetry | Generator | Action on coordinates |
|----------|-----------|----------------------|
| **Translation** | $\delta \in \mathbb{Z}$ | $x_i \to \operatorname{clamp}(x_i + \delta)$ |
| **Rotation** | $90° \times k$ | $(x, y) \to (-y, x)$ for $k=1$, etc. |
| **Reflection(X)** | — | $(x, y) \to (-x, y)$ |
| **Reflection(Y)** | — | $(x, y) \to (x, -y)$ |

After every transformation values are clamped back to $\{-1, 0, +1\}$, ensuring the group action closes on the discrete state space.

### Conserved quantities

| Symmetry | Conserved quantity | Formula |
|----------|-------------------|---------|
| Time translation | Energy | $E = \sum_i \bigl(p_i^2/2 + x_i^2/2\bigr)$ |
| Space translation | Momentum | $P = \sum_i p_i$ |
| Rotation | Angular momentum | $L = \sum_i (x_i \, p_{y_i} - y_i \, p_{x_i})$ |
| Reflection | Momentum (parity) | $P = \sum_i p_i$ |

## Architecture

```text
┌─────────────────────────────────────────┐
│           Physical concept              │
├─────────────────────────────────────────┤
│  Translation symmetry                  │──►│ DiscreteSymmetry::Translation  │
│  Rotation symmetry                     │──►│ DiscreteSymmetry::Rotation     │
│  Reflection symmetry (X/Y)             │──►│ DiscreteSymmetry::Reflection   │
│  Conserved quantity Q                  │──►│ ConservedQuantity              │
│  State snapshot (q, p, t)              │──►│ TernaryState                   │
│  Verify E / P / L constant             │──►│ Verification                   │
└─────────────────────────────────────────┘   └────────────────────────────────┘
```

## Getting Started

Add to `Cargo.toml`:

```toml
[dependencies]
ternary-noether = "0.1"
```

Create a state, apply a rotation, and read off the conserved angular momentum:

```rust
use ternary_noether::{TernaryState, DiscreteSymmetry, ConservedQuantity, Verification};

fn main() {
    // Particle at (1, 0) with momentum (0, 1): pure rotation about origin
    let state = TernaryState::new(vec![1, 0], vec![0, 1]);

    // Rotation by 90° is a symmetry of the isotropic harmonic system
    let rot = DiscreteSymmetry::Rotation(90);
    let angular_mom = ConservedQuantity::derive_from_symmetry(&rot, &state);
    println!("{} = {}", angular_mom.name, angular_mom.value); // L = 1

    // Verify energy is unchanged if we duplicate the state
    let seq = vec![state.clone(), state.clone()];
    assert!(Verification::verify_energy_conservation(&seq));
    println!("Energy is conserved across the sequence.");
}
```

## Running the Tests

Run the 24-test suite with `cargo test`. Each test maps a symmetry to its conservation law:

| Test group | What it verifies |
|------------|------------------|
| `translation_*` (4 tests) | Translation adds $\delta$ and clamps; zero translation is the identity; clamping respects $\{-1,0,+1\}$ bounds. |
| `reflection_*` (4 tests) | X-reflection negates even indices, Y-reflection negates odd indices; applying any reflection twice returns the original state ($\sigma^2 = \mathbb{1}$). |
| `rotation_*` (4 tests) | 180° is point inversion; four 90° steps yield identity; 0° is identity. Group closure is respected modulo clamping. |
| `ternary_state_*` (3 tests) | Energy $E = \sum (p^2 + x^2)/2$ computed correctly; zero state has zero energy; out-of-range values are clamped. |
| `derive_*` (2 tests) | Translation symmetry yields total momentum; rotation symmetry yields $L = x \, p_y - y \, p_x$. |
| `time_translation_*` (1 test) | Time-translation symmetry maps to energy conservation. |
| `space_translation_*` (1 test) | Space-translation symmetry maps to linear momentum conservation. |
| `rotation_conserved_quantity_*` (1 test) | Rotational symmetry maps to angular momentum conservation for multi-pair systems. |
| `verify_*` (4 tests) | `Verification` detects constant and varying energy, momentum, and angular momentum; empty sequences are trivially conserved. |

## Related crates

- [`ternary-hamiltonian`](https://github.com/phoenix/ternary-hamiltonian) — Hamiltonian mechanics and symplectic integration on ternary phase space
- [`ternary-electromagnetism`](https://github.com/phoenix/ternary-electromagnetism) — Maxwell's equations and Yee-lattice wave propagation
- [`ternary-symmetry`](https://github.com/phoenix/ternary-symmetry) — Group actions and orbit enumeration on discrete alphabets
- [`ternary-dynamics`](https://github.com/phoenix/ternary-dynamics) — General discrete dynamical systems on {-1, 0, +1}

## License

MIT
