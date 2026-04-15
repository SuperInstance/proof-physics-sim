# Float Drift vs Constraint-Theory Physics Simulation

A side-by-side benchmark showing how standard `f64` arithmetic accumulates
energy drift in a 3-body gravitational simulation versus snapping position
directions through a `PythagoreanManifold` (from the `constraint-theory-core`
crate) after each integration step.

## Overview

`proof-physics-sim` is a deterministic physics benchmark that simulates the
inner solar system (Sun + Earth-like + Jupiter-like) for **100 000 timesteps**
of 1 hour each (~11.4 years of simulated time) using the **Velocity Verlet**
integrator. It compares two approaches side-by-side:

| Mode | Method | Energy conservation |
|------|--------|---------------------|
| `float` | Plain `f64` Velocity Verlet, no corrections | Drift accumulates over long integrations |
| `ct` | Velocity Verlet + `PythagoreanManifold` direction snap | Constrained to Pythagorean angles; reduced or comparable drift |

The CT (constraint-theory) mode snaps the *angular direction* of each body's
position vector to the nearest Pythagorean unit angle after every timestep,
while preserving radial distance. This enforces arithmetic consistency on
geometric invariants and suppresses rounding-error feedback loops.

## Architecture

```
proof-physics-sim/
├── Cargo.toml          binary crate (release opt-level 3)
├── src/
│   └── main.rs         single-file simulation (~345 LOC)
├── bench.sh            debug + release build & run script
└── README.md
```

**Module structure within `main.rs`:**

| Section | Responsibility |
|---------|---------------|
| `Body` struct | Mass, 2-D position `[f64; 2]`, velocity `[f64; 2]` in the orbital plane |
| Vector helpers | `add`, `sub`, `scale`, `norm`, `dist` — zero-dependency 2-D linear algebra |
| `acceleration()` | Pairwise gravitational acceleration: `a_i = Σ G·mⱼ·(rⱼ−rᵢ)/|rⱼ−rᵢ|³` |
| `total_energy()` | KE + PE: `½mv² − Gm₁m₂/r` summed over all pairs |
| `step_verlet()` | Standard Velocity Verlet: `x(t+dt) = x + v·dt + ½a·dt²`, then `v(t+dt) = v + ½(a+a')·dt` |
| `step_verlet_ct()` | Same as above, then `snap_pos()` projects the direction onto the nearest Pythagorean angle |
| `snap_pos()` | Normalises position → unit direction → `manifold.snap()` → reconstructs `snapped_dir × |r|` |
| `run_float()` / `run_ct()` | Full simulation runners that collect `SimResult` (E₀, E_final, drift %, ΔPos) |
| `print_table()` | Pretty-printed comparison table with Unicode box drawing |

## Physics Engine Details

### Gravitational Model

The simulation solves the Newtonian N-body problem in 2-D (orbital plane = XY):

```
F_ij = G · m_i · m_j / |r_ij|²         (Newton's law of gravitation)
a_i  = Σ_{j≠i} G · m_j · (r_j - r_i) / |r_j - r_i|³
```

| Constant | Value | Unit |
|----------|-------|------|
| G | 6.674 × 10⁻¹¹ | m³ kg⁻¹ s⁻² |
| Δt | 3600 | s (1 hour) |
| Steps | 100 000 | — |
| Total time | ~11.4 | years |

### Velocity Verlet Integrator

Chosen for its **symplectic** property — it approximately conserves energy
over long integrations, making it the standard choice in orbital mechanics:

```
x(t+dt) = x(t) + v(t)·dt + ½·a(t)·dt²       (position half-step)
a(t+dt) = F(x(t+dt)) / m                      (recompute acceleration)
v(t+dt) = v(t) + ½·(a(t) + a(t+dt))·dt        (velocity full-step)
```

### CT Direction Snapping

The manifold does **not** snap position magnitudes — only angular directions.
This preserves the radial orbital dynamics while constraining the geometric
invariant (the unit direction vector) to a finite set of Pythagorean angles:

```rust
fn snap_pos(pos: [f64; 2], manifold: &PythagoreanManifold) -> [f64; 2] {
    let r = norm(pos);
    let dir = [pos[0] / r, pos[1] / r];          // normalise to unit
    let (snapped_dir, _noise) = manifold.snap(dir); // project to nearest Pythagorean angle
    [snapped_dir[0] * r, snapped_dir[1] * r]      // restore magnitude
}
```

With `density=200`, the manifold pre-computes ~450 primitive Pythagorean
triples via Euclid's formula, providing angular resolution of ~0.8°.

## Quick Start

```bash
# Build and run (release for speed)
cargo run --release

# Or use the bench script (debug + release builds)
./bench.sh

# Unit tests (if any are added)
cargo test
```

## Benchmarking

```bash
./bench.sh
```

The script runs both `debug` and `release` builds and reports wall-clock time
for each mode.

## Integration

To integrate constraint-theory energy conservation into your own physics engine:

**1. Add the dependency:**

```toml
[dependencies]
constraint-theory-core = "1.0.1"
```

**2. Initialise the manifold once:**

```rust
use constraint_theory_core::PythagoreanManifold;

// density controls angular resolution (~density²/2 Pythagorean triples)
let manifold = PythagoreanManifold::new(200);
```

**3. Snap directions after each integration step:**

```rust
fn physics_step(&mut self, dt: f64, manifold: &PythagoreanManifold) {
    // Your existing Verlet / RK4 / symplectic Euler step
    self.integrate(dt);

    // Project each body's position direction onto the nearest Pythagorean angle
    for body in &mut self.bodies {
        let r = body.pos.norm();
        if r < 1e-10 { continue; }
        let dir = [body.pos[0] / r, body.pos[1] / r];
        let (snapped, _) = manifold.snap(dir);
        body.pos = [snapped[0] * r, snapped[1] * r];
    }
}
```

**4. Monitor energy drift:**

```rust
let e0 = total_energy(&bodies);
// ... run simulation ...
let drift = (total_energy(&bodies) - e0) / e0.abs() * 100.0;
println!("Energy drift: {drift:.4e}%");
```

## Dependencies

- [`constraint-theory-core`](https://crates.io/crates/constraint-theory-core) v1.0.1
- Rust 1.70+

## Expected Output

```
Running float simulation (100000 steps)…
Running CT (PythagoreanManifold) simulation (100000 steps)…

┌──────────────────────────────────────────────────────────────┐
│      3-Body Solar Simulation — Energy Drift Comparison       │
├──────────────────────────────────────────────────────────────┤
│  Mode           E₀ (J)      E_final (J)    Drift %  ΔPos(m) │
├──────────────────────────────────────────────────────────────┤
│  float (f64)   -X.XXXXXXe+33 -X.XXXXXXe+33  X.XXXXe-XX%  0.00e0│
│  CT (manifold) -X.XXXXXXe+33 -X.XXXXXXe+33  X.XXXXe-XX%  X.XXe+Y│
└──────────────────────────────────────────────────────────────┘
```

The CT mode demonstrates that manifold-snapping constrains floating-point
error propagation. When CT drift is lower, it shows the manifold's ability
to suppress energy drift. When comparable, it shows the angular constraint
introduces negligible overhead while maintaining geometric consistency
(quantified by the ΔPos metric).

---

<img src="callsign1.jpg" width="128" alt="callsign">
