//! # ternary-noether
//!
//! Noether's theorem for discrete ternary systems.
//!
//! This crate derives conserved quantities from discrete symmetries acting on
//! states defined over the ternary alphabet {-1, 0, +1}. It is a discrete
//! analogue of Emmy Noether's 1915 theorem, which connects continuous
//! symmetries of a Lagrangian to conservation laws in classical mechanics.

/// Axis for reflection symmetries.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Axis {
    X,
    Y,
}

/// A discrete symmetry that can act on ternary coordinate vectors.
///
/// All coordinates live in {-1, 0, +1}. After any transformation, values are
/// clamped back into this set.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum DiscreteSymmetry {
    /// Shift every coordinate by `delta`, clamping to {-1, 0, +1}.
    Translation(i8),
    /// Rotate in 2D by 90, 180, or 270 degrees. Coordinates are treated as
    /// (x, y) pairs; odd-length slices leave the last element unchanged.
    Rotation(u8),
    /// Reflect through the given axis.
    Reflection(Axis),
}

/// Clamp a value into {-1, 0, +1}.
#[inline]
fn clamp_ternary(v: i8) -> i8 {
    v.max(-1).min(1)
}

impl DiscreteSymmetry {
    /// Apply this symmetry transformation to a slice of ternary coordinates.
    ///
    /// The output vector has the same length as the input slice. Every output
    /// value is guaranteed to lie in {-1, 0, +1}.
    pub fn apply(&self, coords: &[i8]) -> Vec<i8> {
        match self {
            DiscreteSymmetry::Translation(delta) => {
                coords.iter().map(|&c| clamp_ternary(c + delta)).collect()
            }

            DiscreteSymmetry::Rotation(degrees) => {
                // Normalise to one of four canonical angles.
                let steps = ((*degrees as u32 / 90) % 4) as u8;
                let mut result = coords.to_vec();
                // Process (x, y) pairs; ignore trailing element if odd length.
                let pairs = result.len() / 2;
                for i in 0..pairs {
                    let x = result[2 * i];
                    let y = result[2 * i + 1];
                    // Each step is a 90° CCW rotation: (x, y) → (-y, x).
                    let (nx, ny) = match steps {
                        0 => (x, y),
                        1 => (clamp_ternary(-y), clamp_ternary(x)),
                        2 => (clamp_ternary(-x), clamp_ternary(-y)),
                        3 => (clamp_ternary(y), clamp_ternary(-x)),
                        _ => unreachable!(),
                    };
                    result[2 * i] = nx;
                    result[2 * i + 1] = ny;
                }
                result
            }

            DiscreteSymmetry::Reflection(axis) => match axis {
                Axis::X => {
                    // Negate x-coordinates (indices 0, 2, 4, …) in (x,y) pairs.
                    let pairs = coords.len() / 2;
                    let mut result = coords.to_vec();
                    for i in 0..pairs {
                        result[2 * i] = clamp_ternary(-result[2 * i]);
                    }
                    // Odd trailing element: unchanged (not part of a pair).
                    result
                }
                Axis::Y => {
                    // Negate y-coordinates (indices 1, 3, 5, …) in (x,y) pairs.
                    let pairs = coords.len() / 2;
                    let mut result = coords.to_vec();
                    for i in 0..pairs {
                        result[2 * i + 1] = clamp_ternary(-result[2 * i + 1]);
                    }
                    result
                }
            },
        }
    }
}

/// A conserved quantity derived from a symmetry of the system.
#[derive(Debug, Clone, PartialEq)]
pub struct ConservedQuantity {
    /// Human-readable name, e.g. "energy" or "momentum".
    pub name: String,
    /// Numerical value of the conserved quantity.
    pub value: f64,
}

impl ConservedQuantity {
    /// Construct a conserved quantity directly.
    pub fn new(name: impl Into<String>, value: f64) -> Self {
        Self {
            name: name.into(),
            value,
        }
    }

    /// Derive a conserved quantity from a symmetry and a ternary state.
    ///
    /// The mapping follows the discrete analogue of Noether's theorem:
    ///
    /// | Symmetry                    | Conserved quantity  |
    /// |-----------------------------|---------------------|
    /// | `Translation(_)`            | total momentum      |
    /// | `Rotation(_)`               | angular momentum    |
    /// | `Reflection(Axis::X or Y)`  | total momentum      |
    pub fn derive_from_symmetry(symmetry: &DiscreteSymmetry, state: &TernaryState) -> Self {
        match symmetry {
            DiscreteSymmetry::Translation(_) => {
                let momentum: i32 = state.momenta.iter().map(|&m| m as i32).sum();
                ConservedQuantity::new("momentum", momentum as f64)
            }
            DiscreteSymmetry::Rotation(_) => {
                // L = Σ (x_i * p_y_i − y_i * p_x_i) over (x,y) position pairs
                // paired with (px, py) momentum pairs.
                let pos_pairs = state.positions.len() / 2;
                let mom_pairs = state.momenta.len() / 2;
                let pairs = pos_pairs.min(mom_pairs);
                let l: f64 = (0..pairs)
                    .map(|i| {
                        let x = state.positions[2 * i] as f64;
                        let y = state.positions[2 * i + 1] as f64;
                        let px = state.momenta[2 * i] as f64;
                        let py = state.momenta[2 * i + 1] as f64;
                        x * py - y * px
                    })
                    .sum();
                ConservedQuantity::new("angular_momentum", l)
            }
            DiscreteSymmetry::Reflection(_) => {
                let momentum: i32 = state.momenta.iter().map(|&m| m as i32).sum();
                ConservedQuantity::new("momentum", momentum as f64)
            }
        }
    }
}

/// A snapshot of a discrete ternary dynamical system.
///
/// All positions and momenta must lie in {-1, 0, +1}. Values outside this
/// range are silently clamped on construction.
#[derive(Debug, Clone, PartialEq)]
pub struct TernaryState {
    pub positions: Vec<i8>,
    pub momenta: Vec<i8>,
    pub time: f64,
}

impl TernaryState {
    /// Construct a new state, clamping all values to {-1, 0, +1}.
    pub fn new(positions: Vec<i8>, momenta: Vec<i8>) -> Self {
        Self {
            positions: positions.into_iter().map(clamp_ternary).collect(),
            momenta: momenta.into_iter().map(clamp_ternary).collect(),
            time: 0.0,
        }
    }

    /// Discrete kinetic + potential energy.
    ///
    /// E = Σ p²/2 + Σ x²/2, where each p and x is cast to f64 first.
    pub fn energy(&self) -> f64 {
        let kinetic: f64 = self.momenta.iter().map(|&p| (p as f64).powi(2) / 2.0).sum();
        let potential: f64 = self
            .positions
            .iter()
            .map(|&x| (x as f64).powi(2) / 2.0)
            .sum();
        kinetic + potential
    }

    /// Total (scalar) momentum: Σ p_i.
    pub fn total_momentum(&self) -> i32 {
        self.momenta.iter().map(|&m| m as i32).sum()
    }

    /// Angular momentum: Σ (x_i * p_y_i − y_i * p_x_i) over coordinate pairs.
    pub fn angular_momentum(&self) -> f64 {
        let pos_pairs = self.positions.len() / 2;
        let mom_pairs = self.momenta.len() / 2;
        let pairs = pos_pairs.min(mom_pairs);
        (0..pairs)
            .map(|i| {
                let x = self.positions[2 * i] as f64;
                let y = self.positions[2 * i + 1] as f64;
                let px = self.momenta[2 * i] as f64;
                let py = self.momenta[2 * i + 1] as f64;
                x * py - y * px
            })
            .sum()
    }
}

// ---------------------------------------------------------------------------
// Symmetry-specific helpers that mirror the classical Noether taxonomy.
// ---------------------------------------------------------------------------

/// Time-translation symmetry → energy conservation.
pub struct TimeTranslation;

impl TimeTranslation {
    /// Return the symmetry associated with time translation.
    ///
    /// Time translation is represented as a unit translation in a
    /// one-dimensional ternary "time" coordinate.
    pub fn symmetry() -> DiscreteSymmetry {
        DiscreteSymmetry::Translation(1)
    }

    /// Derive the energy (conserved quantity dual to time translation).
    pub fn conserved_quantity(state: &TernaryState) -> ConservedQuantity {
        ConservedQuantity::new("energy", state.energy())
    }
}

/// Space-translation symmetry → momentum conservation.
pub struct SpaceTranslation;

impl SpaceTranslation {
    /// Return the spatial translation symmetry for displacement `delta`.
    pub fn symmetry(delta: i8) -> DiscreteSymmetry {
        DiscreteSymmetry::Translation(delta)
    }

    /// Derive the total momentum (conserved quantity dual to space translation).
    pub fn conserved_quantity(state: &TernaryState) -> ConservedQuantity {
        let momentum: i32 = state.total_momentum();
        ConservedQuantity::new("momentum", momentum as f64)
    }
}

/// Rotational symmetry → angular momentum conservation.
pub struct Rotation;

impl Rotation {
    /// Return the rotational symmetry for `degrees` ∈ {90, 180, 270, 360}.
    pub fn symmetry(degrees: u8) -> DiscreteSymmetry {
        DiscreteSymmetry::Rotation(degrees)
    }

    /// Derive the angular momentum (conserved quantity dual to rotation).
    pub fn conserved_quantity(state: &TernaryState) -> ConservedQuantity {
        ConservedQuantity::new("angular_momentum", state.angular_momentum())
    }
}

/// Tools for verifying conservation laws across a sequence of states.
pub struct Verification;

impl Verification {
    /// Return `true` when energy is the same in every state (within ε = 0.01).
    pub fn verify_energy_conservation(states: &[TernaryState]) -> bool {
        if states.is_empty() {
            return true;
        }
        let e0 = states[0].energy();
        states.iter().all(|s| (s.energy() - e0).abs() < 0.01)
    }

    /// Return `true` when total (scalar) momentum is constant across all states.
    pub fn verify_momentum_conservation(states: &[TernaryState]) -> bool {
        if states.is_empty() {
            return true;
        }
        let p0 = states[0].total_momentum();
        states.iter().all(|s| s.total_momentum() == p0)
    }

    /// Return `true` when angular momentum is constant across all states
    /// (within ε = 0.01).
    pub fn verify_angular_momentum(states: &[TernaryState]) -> bool {
        if states.is_empty() {
            return true;
        }
        let l0 = states[0].angular_momentum();
        states.iter().all(|s| (s.angular_momentum() - l0).abs() < 0.01)
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // ------------------------------------------------------------------
    // DiscreteSymmetry::Translation
    // ------------------------------------------------------------------

    #[test]
    fn translation_applies_delta() {
        let sym = DiscreteSymmetry::Translation(1);
        let result = sym.apply(&[-1, 0, 1]);
        assert_eq!(result, vec![0, 1, 1]); // last clamped from 2 to 1
    }

    #[test]
    fn translation_clamps_to_ternary() {
        let sym = DiscreteSymmetry::Translation(1);
        // 1 + 1 = 2, clamped to 1
        assert_eq!(sym.apply(&[1]), vec![1]);
        // -1 + (-1) = -2, clamped to -1
        let sym_neg = DiscreteSymmetry::Translation(-1);
        assert_eq!(sym_neg.apply(&[-1]), vec![-1]);
    }

    #[test]
    fn translation_zero_is_identity() {
        let sym = DiscreteSymmetry::Translation(0);
        let coords = vec![-1, 0, 1];
        assert_eq!(sym.apply(&coords), coords);
    }

    // ------------------------------------------------------------------
    // DiscreteSymmetry::Reflection
    // ------------------------------------------------------------------

    #[test]
    fn reflection_x_negates_x_coords() {
        let sym = DiscreteSymmetry::Reflection(Axis::X);
        // pairs: (1, -1), (-1, 0)
        let result = sym.apply(&[1, -1, -1, 0]);
        // x-coords negated: -1, 1 ; y-coords unchanged: -1, 0
        assert_eq!(result, vec![-1, -1, 1, 0]);
    }

    #[test]
    fn reflection_y_negates_y_coords() {
        let sym = DiscreteSymmetry::Reflection(Axis::Y);
        // pairs: (1, -1), (0, 1)
        let result = sym.apply(&[1, -1, 0, 1]);
        // x-coords unchanged: 1, 0 ; y-coords negated: 1, -1
        assert_eq!(result, vec![1, 1, 0, -1]);
    }

    #[test]
    fn reflection_x_double_is_identity() {
        let sym = DiscreteSymmetry::Reflection(Axis::X);
        let coords = vec![1, -1, -1, 0, 0, 1];
        let once = sym.apply(&coords);
        let twice = sym.apply(&once);
        assert_eq!(twice, coords);
    }

    #[test]
    fn reflection_y_double_is_identity() {
        let sym = DiscreteSymmetry::Reflection(Axis::Y);
        let coords = vec![1, -1, -1, 0, 0, 1];
        let once = sym.apply(&coords);
        let twice = sym.apply(&once);
        assert_eq!(twice, coords);
    }

    // ------------------------------------------------------------------
    // DiscreteSymmetry::Rotation
    // ------------------------------------------------------------------

    #[test]
    fn rotation_180_is_double_negation() {
        let sym = DiscreteSymmetry::Rotation(180);
        // (x, y) → (-x, -y) for 180°
        let coords = vec![1, -1, 0, 1, -1, 0];
        let result = sym.apply(&coords);
        let expected: Vec<i8> = coords.iter().map(|&v| -v).collect();
        assert_eq!(result, expected);
    }

    #[test]
    fn rotation_90_applied_four_times_is_identity() {
        // Applying 90° four times must return to the original.
        let sym90 = DiscreteSymmetry::Rotation(90);
        let coords = vec![1, 0, -1, 1];
        let r1 = sym90.apply(&coords);
        let r2 = sym90.apply(&r1);
        let r3 = sym90.apply(&r2);
        let r4 = sym90.apply(&r3);
        assert_eq!(r4, coords);
    }

    #[test]
    fn rotation_0_is_identity() {
        // 0° (i.e. 0 mod 360) is the identity.
        let sym = DiscreteSymmetry::Rotation(0);
        let coords = vec![1, -1, 0, 1];
        assert_eq!(sym.apply(&coords), coords);
    }

    // ------------------------------------------------------------------
    // TernaryState
    // ------------------------------------------------------------------

    #[test]
    fn ternary_state_energy_calculation() {
        // positions = [1, -1], momenta = [0, 1]
        // E = (0²/2 + 1²/2) + (1²/2 + (-1)²/2)
        //   = (0 + 0.5) + (0.5 + 0.5) = 1.5
        let state = TernaryState::new(vec![1, -1], vec![0, 1]);
        let e = state.energy();
        assert!((e - 1.5).abs() < 1e-10, "Expected 1.5, got {e}");
    }

    #[test]
    fn ternary_state_zero_energy() {
        let state = TernaryState::new(vec![0, 0], vec![0, 0]);
        assert_eq!(state.energy(), 0.0);
    }

    #[test]
    fn ternary_state_clamps_values() {
        let state = TernaryState::new(vec![5, -5], vec![3, -3]);
        assert_eq!(state.positions, vec![1, -1]);
        assert_eq!(state.momenta, vec![1, -1]);
    }

    // ------------------------------------------------------------------
    // ConservedQuantity::derive_from_symmetry
    // ------------------------------------------------------------------

    #[test]
    fn derive_translation_gives_momentum() {
        let state = TernaryState::new(vec![1, -1], vec![-1, 1]);
        let sym = DiscreteSymmetry::Translation(1);
        let cq = ConservedQuantity::derive_from_symmetry(&sym, &state);
        assert_eq!(cq.name, "momentum");
        // total momentum = -1 + 1 = 0
        assert_eq!(cq.value, 0.0);
    }

    #[test]
    fn derive_rotation_gives_angular_momentum() {
        let state = TernaryState::new(vec![1, 0], vec![0, 1]);
        let sym = DiscreteSymmetry::Rotation(90);
        let cq = ConservedQuantity::derive_from_symmetry(&sym, &state);
        assert_eq!(cq.name, "angular_momentum");
        // L = x*py - y*px = 1*1 - 0*0 = 1
        assert!((cq.value - 1.0).abs() < 1e-10);
    }

    // ------------------------------------------------------------------
    // TimeTranslation / SpaceTranslation / Rotation helpers
    // ------------------------------------------------------------------

    #[test]
    fn time_translation_conserved_quantity_is_energy() {
        let state = TernaryState::new(vec![1, 0], vec![0, 1]);
        let cq = TimeTranslation::conserved_quantity(&state);
        assert_eq!(cq.name, "energy");
        assert!((cq.value - state.energy()).abs() < 1e-10);
    }

    #[test]
    fn space_translation_conserved_quantity_is_momentum() {
        let state = TernaryState::new(vec![1, -1], vec![1, -1]);
        let cq = SpaceTranslation::conserved_quantity(&state);
        assert_eq!(cq.name, "momentum");
        assert_eq!(cq.value, 0.0);
    }

    #[test]
    fn rotation_conserved_quantity_is_angular_momentum() {
        let state = TernaryState::new(vec![1, 0, 0, 1], vec![0, 1, -1, 0]);
        let cq = Rotation::conserved_quantity(&state);
        assert_eq!(cq.name, "angular_momentum");
        // pair 0: L = 1*1 - 0*0 = 1
        // pair 1: L = 0*0 - 1*(-1) = 1
        // total  = 2
        assert!((cq.value - 2.0).abs() < 1e-10, "Expected 2.0, got {}", cq.value);
    }

    // ------------------------------------------------------------------
    // Verification
    // ------------------------------------------------------------------

    #[test]
    fn verify_energy_conservation_constant_sequence() {
        // All states have the same positions and momenta → same energy.
        let s0 = TernaryState::new(vec![1, 0], vec![0, 1]);
        let s1 = TernaryState::new(vec![1, 0], vec![0, 1]);
        let s2 = TernaryState::new(vec![1, 0], vec![0, 1]);
        assert!(Verification::verify_energy_conservation(&[s0, s1, s2]));
    }

    #[test]
    fn verify_energy_conservation_detects_change() {
        let s0 = TernaryState::new(vec![1, 0], vec![0, 1]);
        // Energy changes: different momenta/positions.
        let s1 = TernaryState::new(vec![1, 1], vec![1, 1]);
        assert!(!Verification::verify_energy_conservation(&[s0, s1]));
    }

    #[test]
    fn verify_momentum_conservation_constant_sequence() {
        let s0 = TernaryState::new(vec![0], vec![1]);
        let s1 = TernaryState::new(vec![1], vec![1]);
        let s2 = TernaryState::new(vec![-1], vec![1]);
        // Total momentum = 1 in all states.
        assert!(Verification::verify_momentum_conservation(&[s0, s1, s2]));
    }

    #[test]
    fn verify_momentum_conservation_detects_change() {
        let s0 = TernaryState::new(vec![0], vec![1]);
        let s1 = TernaryState::new(vec![0], vec![-1]);
        assert!(!Verification::verify_momentum_conservation(&[s0, s1]));
    }

    #[test]
    fn verify_angular_momentum_constant_sequence() {
        // Same state repeated → constant angular momentum.
        let s0 = TernaryState::new(vec![1, 0], vec![0, 1]);
        let s1 = TernaryState::new(vec![1, 0], vec![0, 1]);
        assert!(Verification::verify_angular_momentum(&[s0, s1]));
    }

    #[test]
    fn verify_empty_slice_is_trivially_conserved() {
        let empty: &[TernaryState] = &[];
        assert!(Verification::verify_energy_conservation(empty));
        assert!(Verification::verify_momentum_conservation(empty));
        assert!(Verification::verify_angular_momentum(empty));
    }
}
