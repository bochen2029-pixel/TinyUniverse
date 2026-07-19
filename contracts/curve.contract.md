# CONTRACT — curve (v2 N3: geometry curves) · v1.1.0 · status: FROZEN 2026-07-12 (N3 CLOSED — 3/3 goldens green: deflect `4e6c33ca` · precess `67272705` · shapiro `20bfd4d2`; D-024)

**Purpose.** The third rung of the v2 SUBSTRATE ladder (`ROADMAP.md` §3, N3): the **spatial** metric. N2 `lapse` curved *time* (the lapse α = √(1+2Φ/c²) → gravitational redshift, exact). N3 curves *space*: the substrate's potential Φ sources the full weak-field metric

```
ds² = −(1 + 2Φ/c²) c²dt² + (1 − 2Φ/c²)(dx² + dy² + dz²)
```

and **geodesics through that metric bend light and precess orbits at the exact General-Relativity rates.**

**The physics claim (the emergence), and why it's the cleanest possible test of spatial curvature — the factor of 2.** A photon passing a mass at impact parameter b is deflected by **Δφ = 4GM/(bc²)** (Eddington 1919). This splits into two equal halves:
- **half from time** (the lapse g_tt — exactly N2): a "Newtonian" photon falling through the potential deflects by 2GM/(bc²);
- **half from space** (the spatial factor 1−2Φ/c² — this rung, N3): the curvature of space adds the other 2GM/(bc²).

So N3's signature golden literally *measures the doubling*: integrate the same photon through the metric with (a) only the lapse curved → **2GM/bc²**, and (b) the full metric → **4GM/bc²**. The spatial curvature is exactly the difference. This is the famous 1919 result, decomposed into the two rungs of the substrate ladder. Orbits get the companion observable: perihelion precession **Δϖ = 6πGM/(c²a(1−e²))** per orbit (Einstein/Mercury 1915).

Single-file CPU fp64 oracle tool `substrate/curve_nexus.cpp` — the `_nexus` family idiom (blake2b golden hasher; envelope face `--scenario X --json|--golden|--selftest`; exit 0/1/2). **No GPU** (geodesic ODE integration is a few fp64 ODEs — deterministic, exact, runs under any card contention, like N0 `substrate_nexus`). The metric factors (1±2Φ/c²) are the substrate's Poisson potential; N2 already validated g_tt *on the lattice* (`redshiftPM`, 3.6%), so N3 establishes the **curvature observables** (deflection, precession) against exact GR via the geodesic oracle.

## Units & dials (v1 dials — c LIVE, as in N2)

| dial | value | role in N3 |
|---|---|---|
| c | 20 su/s | the metric's c; deflection 4GM/bc², precession 6πGM/c²a(1−e²); r_s = 2GM/c² |
| G | 2×10⁻³ | with M sets GM/c² = r_s/2 (the deflection/precession amplitude) |
| M | 10⁵ sm | the lensing/central mass → r_s = 1 su (weak field at the probed b, a) |
| ħ, m, dt, L_box | (frozen) | **inert in N3** (geodesic oracle; no field evolution) — declared, unused |

## Declared engine (the geodesic oracle)

**Light deflection — the refractive-index ray integral.** For the static isotropic weak-field metric, null geodesics obey Fermat's principle with an effective index of refraction **n(r) = √(g_ij/(−g_tt)) = √((1−2Φ/c²)/(1+2Φ/c²)) ≈ 1 − 2Φ/c² = 1 + 2GM/(rc²)** (Φ = −GM/r). The ray is integrated by the eikonal ODE in the orbital plane (arc-length s):

```
dx/ds = t̂ ,   dt̂/ds = [ ∇n − (t̂·∇n) t̂ ] / n        (t̂ stays unit; ∇n = −2GM/(r³c²)·(x,y))
```

RK4, a ray launched at (−X₀, b) with t̂=(1,0), integrated to x=+X₀ (X₀ ≫ b so the residual tail bending is negligible). The **deflection** = the angle the outgoing t̂ has turned. Two indices are run per impact parameter: **n_full = 1 + 2GM/(rc²)** (space+time) and **n_lapse = 1 + GM/(rc²)** (time only — the N2 lapse's contribution); their deflections are 4GM/bc² and 2GM/bc².

**Orbit precession — the relativistic orbit equation.** A bound timelike geodesic in the weak-field metric obeys, to 1PN, the orbit equation `d²u/dφ² + u = GM/L² + (3GM/c²)u²` (u=1/r, L the specific angular momentum). Integrated by RK4 in φ from perihelion to the next perihelion; the advance **Δϖ = φ_period − 2π** is the precession. The 3GM u²/c² term is the relativistic correction; its integral over one orbit is the coordinate-invariant 6πGM/(c²a(1−e²)).

**Determinism.** Pure fp64 RK4 with fixed step counts — same (scenario, dials, params) → byte-identical declared JSON → byte-identical blake2b. No RNG, no GPU, no atomics. **Golden = blake2b-256 of the declared JSON** (the measured deflections/precession *are* the declared state — the nexus/substrate_nexus CPU-oracle idiom). sm-independent (CPU).

## Scenarios & gates (the golden units)

| scenario | what runs | gate | oracle |
|---|---|---|---|
| **deflect** | photon rays at b ∈ {100, 200} su through n_full and n_lapse (M=10⁵ → r_s=1; X₀=5000, RK4) | **`|Δφ_full·b·c²/(4GM) − 1| < 0.02`** AND **`|Δφ_lapse·b·c²/(2GM) − 1| < 0.02`** (both b) — the full metric bends by exact GR 4GM/bc²; the lapse alone by exactly half; **space curvature is the difference** | **analytic GR** Δφ = 4GM/(bc²) (full), 2GM/(bc²) (lapse-only). EXACT weak-field. |
| **precess** | bound orbit a=1000, e=0.5 (p=a(1−e²)=750 ≫ r_s — the weak-field regime where the leading 1PN precession is the exact answer; 0.72°/orbit) integrated perihelion→perihelion | **`|Δϖ_meas/(6πGM/(c²·a(1−e²))) − 1| < 0.02`** | **analytic GR** perihelion precession 6πGM/(c²a(1−e²)). |
| **shapiro** | a photon ray past the mass at b=100 (X₀=5000); accumulate the excess coordinate light-travel time (1/c)∫(n−1)dl along the geodesic | **`|Δt_meas/((2GM/c³)·2·asinh(X₀/b)) − 1| < 0.02`** | **analytic GR** Shapiro delay Δt = (4GM/c³)·asinh(X₀/b) ≈ (4GM/c³)ln(2X₀/b) — the **4th classical test** (Shapiro 1964). |

**`--selftest` (flat):** GM=0 → n≡1 → the ray travels straight (Δφ < 1e-9) and the orbit closes (Δϖ < 1e-6). The no-gravity sanity + determinism smoke.

Both gates are exact against closed-form GR; the 2% tolerance absorbs the weak-field truncation (O((r_s/b)²), O((r_s/p)) — measured, not pre-tuned; D-018 discipline) plus RK4/finite-X₀ residuals. No scenario is gated against an unmeasured number.

## The headless face (envelope-conformant)

```
curve_nexus.exe --scenario deflect|precess [--seed N=20260711] --json    # declared JSON envelope
curve_nexus.exe --scenario NAME --golden                                 # vs goldens/curve_NAME/golden.hash
curve_nexus.exe --selftest                                               # flat: straight ray, closed orbit
```

Declared body (canonical order): `seed, params{scenario, c, G, M, r_s, b/a/e, X0, steps}, result{<defl_full_b50, defl_lapse_b50, ratio_full_b50, ratio_lapse_b50, … | precess_meas, precess_exact, precess_rel>, nan_free}, gates{primary, secondary, nan}, verdict`. Floats fmt6/fmt9. **Exit 0 pass · 1 gate fired · 2 error.** Build with MSVC `cl` (BUILD.md CPU path; clang/g++ parity owed like nexus).

## Oracle pedigree (Invariant 3)

- **Exact analytic GR:** deflection 4GM/bc² (Eddington) + 2GM/bc² lapse-half; precession 6πGM/c²a(1−e²) (Einstein). Oracle-gate failure ⇒ run-state SUSPECT.
- **The substrate tie:** the metric factors (1±2Φ/c²) are the substrate's Poisson Φ; g_tt is lattice-validated by N2 `redshiftPM` (3.6%). N3 adds the spatial factor and the geodesic observables. The v1 cross-check: v1 `photons` (0.83% deflection) and `keprel` (0.50% precession) got these by *pseudo-forces*; N3 gets them from the *metric* — the emergence, not a hard-coded factor-2.

## The honest boundary (permanent, printed)

1. **N3 is the geodesic/curvature ORACLE, not a dynamical GPU metric field.** It proves the static weak-field metric sourced by the substrate's potential bends light by the exact GR factor and precesses orbits. A **live GPU metric field a²(x) with dynamical back-reaction** (energy density sourcing the metric each step, the full Hamiltonian/momentum constraint) is the harder continuation toward N4/horizon — named, not faked.
2. **Weak field, 1PN.** The metric is the leading post-Newtonian form (1±2Φ/c²); the deflection/precession are the leading-order GR observables. Strong-field (near r_s) and higher PN orders are out of scope (N4+). The 2% gate reflects the honest weak-field truncation.
3. **Static, spherical, no frame-dragging.** Schwarzschild-class; Kerr/rotation is the v1 `gargantua` render, not N3.

## Build runbook

1. **Contract (this file).** ← *here.*
2. **`substrate/curve_nexus.cpp`** — blake2b + fmt + the ray/orbit integrators + envelope/golden/main. CPU fp64.
3. **Build** — `cl /std:c++17 /EHsc /O2 /W4 substrate\curve_nexus.cpp /Fe:build\curve_nexus.exe` (vcvars64; BUILD.md CPU path).
4. **Golden** — `goldens/curve_deflect/golden.hash`, `goldens/curve_precess/golden.hash`.
5. **Harness** — a `CURVE_GOLDENS` block runs with the CPU oracles (no preflight — GPU-independent).
6. **Register** — ARCHITECTURE §11, TASKLIST, DECISIONS (D-024), RUN_STATE, MODULE_curve.md.

## Changelog

- v1.1.0 (2026-07-12) — added the **shapiro** scenario: the Shapiro time delay, the 4th classical test of GR. Excess light-travel time (1/c)∫(n−1)dl along the ray vs (4GM/c³)·asinh(X₀/b), matched to 0.33% (`curve_shapiro` `20bfd4d2`). With N2's gravitational redshift, all four classical tests of GR now pass on the substrate. The JSON declared-format is unchanged, so the deflect/precess goldens stay byte-identical (backward-compatible minor bump).
- v1.0.0 (2026-07-12) — initial contract + build (operator directive: proceed N3+). The spatial metric: geodesic light deflection (4GM/bc², decomposed into the N2-lapse half + the N3-space half — the 1919 factor of 2) and perihelion precession (6πGM/c²a(1−e²)), vs exact GR. CPU fp64 geodesic oracle `substrate/curve_nexus.cpp`. Honest boundary: static weak-field geodesic oracle (the curvature observables), not a dynamical GPU metric field (that's N4+).
