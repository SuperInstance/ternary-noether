# PLUG_AND_PLAY — Noether

> Noether's theorem for discrete ternary systems

## 🚀 Quick Start

Add to your `Cargo.toml`:

```toml
[dependencies]
ternary-noether = { git = "https://github.com/SuperInstance/ternary-noether" }
```

Use in your code:

```rust
use ternary_noether::{TernaryState, conserved_quantity};

let state = TernaryState::new(&[1.0, -1.0, 0.0]);
let conserved = conserved_quantity(&state);
```

## 🔗 Integration

This crate is part of the [SuperInstance ternary fleet](https://github.com/SuperInstance). It uses the canonical `Ternary` type from `ternary-types` for cross-crate compatibility.

## 📄 License

MIT
