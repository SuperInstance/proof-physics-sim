# Float Drift vs Constraint-Theory Physics Simulation

A side-by-side benchmark showing how standard `f64` arithmetic accumulates
energy drift in a 3-body gravitational simulation versus snapping positions
through a `PythagoreanManifold` (from the `constraint-theory-core` crate)
after each integration step.

## Experiment

We simulate three gravitationally-interacting bodies for **100,000 timesteps**
(1 hour each ≈ 11.4 years) using the **Velocity Verlet** integrator:

| Body | Stand-in for | Mass |
|------|-------------|------|
| Body 0 | Sun | 1.989 × 10³⁰ kg |
| Body 1 | Earth | 5.972 × 10²⁴ kg |
| Body 2 | Jupiter | 1.898 × 10²⁷ kg |

### Two modes

**`float` mode** — plain `f64` Velocity Verlet, no corrections.

**`ct` mode** — same integrator, but after each position update the raw
position vector is passed through `PythagoreanManifold::snap()`:

```rust
let snapped = manifold.snap(&raw_pos.to_vec());
```

The manifold (tolerance 1e-6) projects coordinates onto nearby Pythagorean
triples, enforcing arithmetic consistency and suppressing the accumulation of
rounding errors that cause long-term energy drift in float arithmetic.

### Metrics reported

| Metric | Description |
|--------|-------------|
| **E₀** | Initial total mechanical energy (J) |
| **E_final** | Final total mechanical energy after all steps |
| **Drift %** | `(E_final − E₀) / |E₀| × 100` |
| **ΔPos (m)** | Mean positional deviation of CT orbit from the float reference |

## Running

```bash
# Build and run (release for speed)
cargo run --release

# Or use the bench script
./bench.sh
```

## Benchmarking

```bash
./bench.sh
```

The script runs both `debug` and `release` builds and reports wall-clock time
for each.

## Dependencies

- [`constraint-theory-core`](https://crates.io/crates/constraint-theory-core) v1.0.1
- Rust 1.70+

## Expected output

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

The CT mode should show meaningfully lower or comparable drift, demonstrating
that manifold-snapping constrains floating-point error propagation.
