/// Float Drift vs Constraint-Theory Physics Simulation
///
/// Demonstrates how standard f64 arithmetic accumulates energy drift in a
/// 3-body gravitational simulation versus snapping position directions to the
/// nearest Pythagorean angle (via PythagoreanManifold) after each integration
/// step.
///
/// API note: PythagoreanManifold::new(density: usize) pre-computes all
/// primitive Pythagorean triples up to parameter `density` via Euclid's
/// formula.  snap([f32; 2]) -> ([f32; 2], f32) normalises the input, finds
/// the nearest Pythagorean unit vector, and returns (snapped_unit, noise).

use constraint_theory_core::PythagoreanManifold;

const G: f64 = 6.674e-11; // gravitational constant (m³ kg⁻¹ s⁻²)
const TIMESTEP: f64 = 3600.0; // 1 hour in seconds
const STEPS: usize = 100_000;

// ─────────────────────────────────────────────────────────────────────────────
// 2-D body (orbital plane = XY plane)
// ─────────────────────────────────────────────────────────────────────────────

#[derive(Clone, Debug)]
struct Body {
    mass: f64,
    pos: [f64; 2],
    vel: [f64; 2],
}

impl Body {
    fn new(mass: f64, pos: [f64; 2], vel: [f64; 2]) -> Self {
        Body { mass, pos, vel }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Vector helpers (2-D)
// ─────────────────────────────────────────────────────────────────────────────

#[inline]
fn add(a: [f64; 2], b: [f64; 2]) -> [f64; 2] {
    [a[0] + b[0], a[1] + b[1]]
}

#[inline]
fn sub(a: [f64; 2], b: [f64; 2]) -> [f64; 2] {
    [a[0] - b[0], a[1] - b[1]]
}

#[inline]
fn scale(v: [f64; 2], s: f64) -> [f64; 2] {
    [v[0] * s, v[1] * s]
}

#[inline]
fn norm_sq(v: [f64; 2]) -> f64 {
    v[0] * v[0] + v[1] * v[1]
}

#[inline]
fn norm(v: [f64; 2]) -> f64 {
    norm_sq(v).sqrt()
}

#[inline]
fn dist(a: [f64; 2], b: [f64; 2]) -> f64 {
    norm(sub(a, b))
}

// ─────────────────────────────────────────────────────────────────────────────
// Physics
// ─────────────────────────────────────────────────────────────────────────────

fn acceleration(bodies: &[Body], i: usize) -> [f64; 2] {
    let mut acc = [0.0_f64; 2];
    for (j, other) in bodies.iter().enumerate() {
        if i == j {
            continue;
        }
        let r = sub(other.pos, bodies[i].pos);
        let d = norm(r);
        let mag = G * other.mass / (d * d * d);
        acc = add(acc, scale(r, mag));
    }
    acc
}

fn total_energy(bodies: &[Body]) -> f64 {
    let mut ke = 0.0_f64;
    let mut pe = 0.0_f64;

    for b in bodies {
        ke += 0.5 * b.mass * norm_sq(b.vel);
    }
    for i in 0..bodies.len() {
        for j in (i + 1)..bodies.len() {
            let d = dist(bodies[i].pos, bodies[j].pos);
            pe -= G * bodies[i].mass * bodies[j].mass / d;
        }
    }

    ke + pe
}

fn mean_pos_deviation(a: &[Body], b: &[Body]) -> f64 {
    let total: f64 = a
        .iter()
        .zip(b.iter())
        .map(|(ai, bi)| dist(ai.pos, bi.pos))
        .sum();
    total / (a.len() as f64)
}

// ─────────────────────────────────────────────────────────────────────────────
// Initial conditions — inner solar system (SI units, 2-D orbital plane)
// ─────────────────────────────────────────────────────────────────────────────

fn initial_bodies() -> Vec<Body> {
    let au = 1.496e11_f64; // 1 AU in metres
    let m_sun = 1.989e30_f64;
    vec![
        // Sun at origin
        Body::new(m_sun, [0.0, 0.0], [0.0, 0.0]),
        // Earth-like at 1 AU — circular orbit
        Body::new(5.972e24, [au, 0.0], [0.0, (G * m_sun / au).sqrt()]),
        // Jupiter-like at 5.2 AU — circular orbit
        Body::new(
            1.898e27,
            [5.2 * au, 0.0],
            [0.0, (G * m_sun / (5.2 * au)).sqrt()],
        ),
    ]
}

// ─────────────────────────────────────────────────────────────────────────────
// Velocity Verlet step — plain f64
// ─────────────────────────────────────────────────────────────────────────────

fn step_verlet(bodies: &mut Vec<Body>, accs: &mut Vec<[f64; 2]>) {
    let n = bodies.len();

    // x(t+dt) = x(t) + v(t)·dt + ½·a(t)·dt²
    for i in 0..n {
        bodies[i].pos = add(
            bodies[i].pos,
            add(scale(bodies[i].vel, TIMESTEP), scale(accs[i], 0.5 * TIMESTEP * TIMESTEP)),
        );
    }

    let new_accs: Vec<[f64; 2]> = (0..n).map(|i| acceleration(bodies, i)).collect();

    // v(t+dt) = v(t) + ½·(a(t)+a(t+dt))·dt
    for i in 0..n {
        bodies[i].vel = add(bodies[i].vel, scale(add(accs[i], new_accs[i]), 0.5 * TIMESTEP));
    }

    *accs = new_accs;
}

// ─────────────────────────────────────────────────────────────────────────────
// CT snap helper
//
// The manifold snaps a 2-D unit direction to the nearest Pythagorean angle.
// We preserve the radial distance (magnitude) and only constrain the *angle*,
// which is the geometric invariant constraint-theory enforces.
// ─────────────────────────────────────────────────────────────────────────────

fn snap_pos(pos: [f64; 2], manifold: &PythagoreanManifold) -> [f64; 2] {
    let r = norm(pos);
    if r < 1e-10 {
        return pos;
    }
    // Normalise to unit direction, convert to f32 for the manifold
    let dir_f32 = [(pos[0] / r) as f32, (pos[1] / r) as f32];
    let (snapped_dir, _noise) = manifold.snap(dir_f32);
    // Reconstruct position: snapped unit direction × original radius
    [snapped_dir[0] as f64 * r, snapped_dir[1] as f64 * r]
}

// ─────────────────────────────────────────────────────────────────────────────
// Velocity Verlet step — CT mode (snap positions after each integration step)
// ─────────────────────────────────────────────────────────────────────────────

fn step_verlet_ct(
    bodies: &mut Vec<Body>,
    accs: &mut Vec<[f64; 2]>,
    manifold: &PythagoreanManifold,
) {
    let n = bodies.len();

    // x(t+dt) = x(t) + v(t)·dt + ½·a(t)·dt²  → then snap
    for i in 0..n {
        let raw_pos = add(
            bodies[i].pos,
            add(scale(bodies[i].vel, TIMESTEP), scale(accs[i], 0.5 * TIMESTEP * TIMESTEP)),
        );
        bodies[i].pos = snap_pos(raw_pos, manifold);
    }

    let new_accs: Vec<[f64; 2]> = (0..n).map(|i| acceleration(bodies, i)).collect();

    // v(t+dt) = v(t) + ½·(a(t)+a(t+dt))·dt
    for i in 0..n {
        bodies[i].vel = add(bodies[i].vel, scale(add(accs[i], new_accs[i]), 0.5 * TIMESTEP));
    }

    *accs = new_accs;
}

// ─────────────────────────────────────────────────────────────────────────────
// Simulation runners
// ─────────────────────────────────────────────────────────────────────────────

struct SimResult {
    mode: &'static str,
    e0: f64,
    e_final: f64,
    drift_pct: f64,
    mean_pos_dev: f64,
}

fn run_float() -> (SimResult, Vec<Body>) {
    let mut bodies = initial_bodies();
    let e0 = total_energy(&bodies);
    let mut accs: Vec<[f64; 2]> = (0..bodies.len()).map(|i| acceleration(&bodies, i)).collect();

    for _ in 0..STEPS {
        step_verlet(&mut bodies, &mut accs);
    }

    let e_final = total_energy(&bodies);
    let drift_pct = (e_final - e0) / e0.abs() * 100.0;

    let result = SimResult {
        mode: "float (f64)",
        e0,
        e_final,
        drift_pct,
        mean_pos_dev: 0.0, // reference orbit
    };

    (result, bodies)
}

fn run_ct(reference_final: &[Body]) -> SimResult {
    let mut bodies = initial_bodies();
    let e0 = total_energy(&bodies);
    let mut accs: Vec<[f64; 2]> = (0..bodies.len()).map(|i| acceleration(&bodies, i)).collect();

    // density=200 gives ~450 Pythagorean states — sufficient angular resolution
    let manifold = PythagoreanManifold::new(200);

    for _ in 0..STEPS {
        step_verlet_ct(&mut bodies, &mut accs, &manifold);
    }

    let e_final = total_energy(&bodies);
    let drift_pct = (e_final - e0) / e0.abs() * 100.0;
    let mean_pos_dev = mean_pos_deviation(&bodies, reference_final);

    SimResult {
        mode: "CT (manifold)",
        e0,
        e_final,
        drift_pct,
        mean_pos_dev,
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Display
// ─────────────────────────────────────────────────────────────────────────────

fn print_table(results: &[SimResult]) {
    let sep_inner = "─".repeat(70);
    println!();
    println!("┌{}┐", sep_inner);
    println!("│{:^70}│", " 3-Body Solar Simulation — Energy & Position Drift ");
    println!("├{}┤", sep_inner);
    println!(
        "│  {:<14} {:>15} {:>13} {:>12}  {:>9} │",
        "Mode", "E₀ (J)", "E_final (J)", "Drift %", "ΔPos (m)"
    );
    println!("├{}┤", sep_inner);

    for r in results {
        println!(
            "│  {:<14} {:>15.6e} {:>13.6e} {:>11.4e}%  {:>9.2e} │",
            r.mode, r.e0, r.e_final, r.drift_pct, r.mean_pos_dev
        );
    }

    println!("└{}┘", sep_inner);

    println!();
    println!("Simulation parameters:");
    println!("  Steps       : {}", STEPS);
    println!("  Timestep    : {} s  ({:.1} hr)", TIMESTEP, TIMESTEP / 3600.0);
    println!(
        "  Total time  : {:.2} years",
        (STEPS as f64 * TIMESTEP) / (365.25 * 86400.0)
    );
    println!("  Bodies      : Sun + Earth-like + Jupiter-like (2-D orbital plane)");
    println!("  Integrator  : Velocity Verlet");
    println!("  CT density  : 200  (~450 Pythagorean states)");
    println!();

    let float_drift = results[0].drift_pct.abs();
    let ct_drift = results[1].drift_pct.abs();
    if ct_drift < float_drift {
        println!(
            "Energy drift reduced by {:.1}× in CT mode vs plain f64.",
            float_drift / ct_drift
        );
    } else if (float_drift - ct_drift).abs() < 1e-15 {
        println!("Both modes produced equivalent energy drift.");
    } else {
        // CT drift larger — expected when angular snapping disturbs orbits more
        // than raw floating-point noise; still a meaningful comparison.
        println!(
            "CT mode orbit constrained to Pythagorean angles (ΔPos = {:.2e} m).",
            results[1].mean_pos_dev
        );
        println!(
            "float drift: {:.4e}%   CT drift: {:.4e}%",
            float_drift, ct_drift
        );
    }
    println!();
}

// ─────────────────────────────────────────────────────────────────────────────
// Entry point
// ─────────────────────────────────────────────────────────────────────────────

fn main() {
    println!("Running float simulation ({} steps)…", STEPS);
    let (float_result, float_final) = run_float();

    println!("Running CT (PythagoreanManifold) simulation ({} steps)…", STEPS);
    let ct_result = run_ct(&float_final);

    print_table(&[float_result, ct_result]);
}
