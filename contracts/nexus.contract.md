# CONTRACT — tiny_nexus · v1.0.0-rc1 · status: UNDER OPERATOR REVIEW

**Purpose.** The composition proof and standing oracle: demonstrate, in fp64 on CPU, that the four regimes (quantum, classical, relativistic, compact) *compose without contradiction* under one dial set, and freeze the dial defaults. Every later CUDA module is gated against nexus values at its declared tolerance. Pattern: `C:\ASTRA-7\proto\astra_nexus.cpp`.

**Tool shape.** Single file `nexus/tiny_nexus.cpp`, C++17, fp64 throughout, no dependencies beyond the standard library. Must compile with MSVC, clang++, and g++ (parity is part of the selftest story). Headless only. Deterministic: the only RNG is seeded `std::mt19937_64` (N6, N11); everything else is closed-form or fixed-step integration.

## CLI

```
tiny_nexus.exe [--dials PATH] [--seed N=20260711] [--json] [--selftest] [--golden] [--only ID]
```

`--dials` overrides the built-in v0 defaults (below). `--only N3` runs one test (dev convenience; never golden). `--selftest` = full battery, human-readable, exit 0/1. `--golden` = battery at frozen defaults, print blake2b of declared JSON, compare `goldens/nexus/`.

## Dials — v0 defaults (frozen by this contract; retuning = MINOR bump pre-golden, MAJOR post)

| dial | value | units |
|---|---|---|
| m_p | 1.0 | sm |
| c | 20.0 | su/s |
| ħ | 0.5 | sm·su²/s |
| G | 2.0e-3 | su³/(sm·s²) |
| dt | 1/240 | s (fixed tick) |
| L_box | 512.0 | su |
| k_B | 1.0 | (defines temperature units) |

## The battery (ordered; each row = declared gate; tolerance in fp64)

| id | name | asserts | oracle |
|---|---|---|---|
| **N1** | dials audit | derived groups computed & bounded: λ_dB(m_p, v=1) ∈ [0.1, 10] su · r_s(10⁵·m_p) = 2GM/c² ∈ [0.5, 5] su · r_s(10²·m_p) < 0.01 su · v_circ(10⁴·m_p @ r=10)/c ∈ [0.03, 0.3] · T_H(M=10⁵) and T_H(M=10³) computed, ratio = 100 exactly · all N-crossovers within 4 decades | closed-form definitions |
| **N2** | kepler | two-body leapfrog (KDK, dt) vs analytic ellipse, e=0.6, 100 orbits: |ΔE/E| < 1e-9 · |ΔL/L| < 1e-12 · Newtonian perihelion drift < 1e-4 rad/orbit | analytic Kepler |
| **N3** | pn1 | with 1PN correction on N2's orbit at compactness GM/(ac²) ≈ 5e-3: measured precession vs Δφ = 6πGM/(a(1−e²)c²) within 5% | analytic 1PN rate |
| **N4** | pw | Paczyński–Wiita Φ = −GM/(r−r_s): innermost stable circular orbit at r = 6GM/c² within 1%; a test orbit at 5.5 GM/c² is captured | P–W analytic ISCO |
| **N5** | qm | 1D split-step (FFT via bundled radix-2, or DFT — fp64): free packet σ(t) = σ₀√(1+(ħt/2mσ₀²)²) within 1% at t where σ doubles · SHO ground state E₀ = ħω/2 within 1e-3 (imaginary-time relaxation) | analytic QM |
| **N6** | ratchet | MC estimate of unwrite probability vs closed form min(1, p/((1−p)ρ))^R at (p=0.3, ρ=0.6, R∈{1,4,16}), 10⁷ trials seeded: within 1% | ORRERY `ratchet` closed form |
| **N7** | rapidity | relativistic kinematics: v(p) < c for p up to 10⁹ · γ via cosh(ω) matches 1/√(1−β²) to 1e-12 at γ ≤ 10³ · at γ = 10⁷ the cosh path is finite & correct while the naive path's cancellation error is *demonstrated and asserted* > 1e-3 (the test proves the discipline is load-bearing) | ASTRA-7 N3 pattern |
| **N8** | tau | proper time: SR circular motion dτ/dt = 1/γ exact to 1e-12 · weak-field static dτ/dt = √(1−2Φ/c²) vs Schwarzschild √(1−r_s/r) within 1e-6 at r ≥ 20 r_s · combined circular-orbit case vs exact √(1−3GM/rc²) within 1e-4 | analytic GR (weak field) |
| **N9** | t_emit | finite-c observation of a Keplerian body from a receding observer: apparent orbital phase computed at t_emit; asserts monotone phase-lag growth and correct classical recession rate (1 − v/c) | ASTRA-7 N4 pattern |
| **N10** | compose | one scene, all regimes, 10⁴ ticks: bound binary (N2 params) + photon in flight + free quantum packet (N5) + one P–W BH far away; asserts: no NaN/Inf · energy bookkeeping across regimes within 1e-6 · regime bitmask stable & correct per body · photon speed = c in all frames checked | all of the above, composed |
| **N11** | det | full battery twice in-process + documented 3× out-of-process: byte-identical declared JSON | envelope discipline |

## Declared output (schema: `nexus.schema.json`)

```
{ "tool":"tiny_nexus", "version":"1.0.0", "seed":20260711,
  "dials":{...},
  "derived":{ "lambda_dB_v1":0.5, "r_s_1e5":1.0, "v_circ_1e4_r10_over_c":0.0707,
              "T_H_1e5":..., "T_H_1e3":..., ... },
  "results":[ {"id":"N1","name":"dials","pass":true,"metrics":{...}}, ... ],
  "verdict":"pass|fail", "notes":"<non-declared>" }
```

Exit 0 = all pass · 1 = any gate fired · 2 = error. Hash domain excludes `notes`.

## Determinism clause

Same `--dials` + `--seed` ⇒ byte-identical declared JSON, cross-compiler *asserted for the integer/logic path*; floating outputs are compiler-pinned per platform in the golden (MSVC is the golden platform; clang/g++ runs must pass all *gates* but may differ in last-ulp metrics — the golden hash is MSVC's).

## Notes for the implementer

- N5's transform: bundle a ~40-line radix-2 FFT; do not link anything.
- N6: seed via splitmix64(seed ^ test_id) — liborrery idiom, hand-inlined (nexus predates the lift).
- Every gate prints measured vs expected on failure. Exit-1 output is a *result*, not a crash.
- Runtime budget: full battery < 30 s (I-selftest speed rule).
