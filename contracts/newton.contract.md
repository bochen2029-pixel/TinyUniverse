# CONTRACT — newton (M2 classical tier) · v1.0.0 · status: FROZEN 2026-07-11 (operator directive "go M2")

**Purpose.** Self-gravity for the particle substrate, deterministic and golden-gated: a PM (particle-mesh) solver for large N and a direct solver for small N, KDK leapfrog integration, conservation meters, and four scenario contracts. Lives inside `app/tinyuniverse.cu` (one binary, two faces — D-003); the headless face is the golden unit.

## The solver (declared behavior)

- **PM (N large):** 128³ grid over the 512-su 3-torus (periodic FFT — the box topology *is* the universe's, ARCHITECTURE §6). CIC deposit into a **fixed-point uint64 grid** (`orrery::fixed_atomic_add`, 2⁻³² mass quantum declared — Invariant 4: no float atomics in declared paths) → cuFFT R2C → Green's function −4πG/k² (k=0 zeroed: mean density does not gravitate — the periodic-gravity convention) → C2R → Φ → centered-difference force grid → CIC gather. No CIC deconvolution (declared; cell = 4 su is the effective softening).
- **Direct (N ≤ 4096):** exact pairwise, fixed iteration order per particle (register accumulation, no atomics), softening ε² = 10⁻⁶ su².
- **Tiny (N ≤ 32):** single-block persistent kernel advancing many ticks per launch (same math as Direct); meters computed host-side in fp64.
- **Integrator:** KDK leapfrog at the dial dt = 1/240 s; position drift uses Kahan compensation (internal buffers; not declared state). fp32 throughout (D-010); τ += dt·(1 − v²/2c²).
- **Determinism:** integer-atomic deposit is order-invariant; cuFFT is bit-stable for a fixed plan/arch; gather/integration are per-particle. Golden = (scenario, dials, seed, steps) → blake2b-256 of the final pos/mom/tau/regime buffers, hardware-pinned sm_89.
- **Meters:** KE, PE, P⃗, L_z via `fixed_atomic_add` (2⁻³² quanta, declared); PE = ½Σ mΦ_PM (includes CIC self-energy — a constant offset; **gates are on drift, never level**). Tiny-N scenarios meter on host in fp64.

## The headless face (envelope-conformant; liborrery serializers)

```
tinyuniverse.exe --scenario NAME [--steps S] [--seed N=20260711] --json   # run + declared JSON envelope
tinyuniverse.exe --scenario NAME --golden                                 # frozen params vs goldens/NAME/declared.hash
```

Declared body (canonical order): `seed, params{scenario, n, steps, solver, c, hbar, G, dt, L_box}, result{state_b2b, e0, e1, de_rel, p_drift, l_drift, n_rel, n_classical}, gates{...pass booleans...}, verdict`. Floats fmt6. Exit 0 pass · 1 gate fired · 2 error. `p_drift` = |P⃗(T)−P⃗(0)| / √(2·M_tot·KE₀); `l_drift` = |L_z(T)−L_z(0)| / max(|L_z(0)|, √(2·M_tot·KE₀)·L_box/8).

## Scenario contracts (frozen golden params; gates are fp32-honest — the *math* is oracle-proven in fp64 by tiny_nexus)

| scenario | N | solver | init (seeded where random) | golden steps | gates | oracle pedigree |
|---|---|---|---|---|---|---|
| **kepler** | 2 | tiny | nexus-N2 params: m={10⁴,1}, a=150, e=0.6, apoapsis start, CM frame | 10⁶ | de_rel < 2e-3 · l_drift < 2e-3 | nexus N2 (fp64: 7.4e-10) |
| **threebody** | 3 | tiny | Chenciner–Montgomery figure-8, m=1each, v scaled ×√G (period ≈ 141.5 s) | 10⁶ (≈29.5 periods) | de_rel < 2e-3 · bounded: max|r| < 5 su | analytic figure-8 stability |
| **cloud** | 10⁶ | pm | cool uniform ball R=100 (dir = normalized gauss³, r = R·u^⅓), sub-virial dispersion σ = 0.35·√(GM/R), T=3800 K | 12000 (≈2 t_ff; t_ff ≈ 24.8 s) | de_rel < 0.08 (relaxation through PM, declared honest) · p_drift < 1e-3 | Jeans/free-fall timescale |
| **galaxy** | 10⁶ | pm | M1 morphology + rotational support: v_circ(r)=√(G·M_enc(r)/r) from the analytic component CDFs (exp disk / Maxwell bulge / log halo, halo capped 240 su), 8% dispersion | 10⁴ | de_rel < 0.02 · p_drift < 1e-3 | virial construction |

Windowed mode runs the same scenarios interactively (`--scenario`, default galaxy); visuals stay non-declared.

## Gate (M2, from TASKLIST)

1M gravitating particles ≥ 60 Hz windowed · scenario goldens frozen + reproduced · long-run drift evidence: kepler + threebody at 10⁶ ticks inside gates (recorded in `runs/`).

## Changelog
- v1.0.0 (2026-07-11) — initial freeze under operator directive; combined solver+scenarios contract (schema: the declared-body field list above is normative; a machine schema follows if an external caller ever consumes this face).
