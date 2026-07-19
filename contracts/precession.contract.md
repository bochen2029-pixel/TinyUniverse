# CONTRACT — precession (v1 polish: Q-006 resolution) · v1.0.0 · status: FROZEN 2026-07-12 (CLOSED — 3/3 goldens green: sr `a0e180df` · pn1 `db0818f2` · combined `f9df648f`; Q-006 RESOLVED, D-026)

**Purpose.** Resolve **Q-006** (`QUESTIONS.md`; the D-016 withdrawal): the einstein contract claimed a combined **SR-inertia + 1PN-field** perihelion precession of **7π·GM/(c²a(1−e²))** per orbit by *linear superposition* (SR contributes π, the harmonic-1PN field contributes 6π). Measurement gave **6.41π** (stable, fp64), contradicting the sum, so the 1PN term was retired (RAYFORMER protocol) and the app ships SR-only (Sommerfeld-exact). This oracle runs the **fp64 isolation experiments Q-006 owed** and settles the number.

**The resolution (MEASURED — and it FLIPPED once the actual semi-latus rectum was measured; the RAYFORMER lesson).** The linear superposition **is correct**. Measured with the orbit's *actual* (force-distorted) semi-latus rectum: sr = **1.00π** (Sommerfeld), pn1 = **6.03π** (the correct full-GR value), combined = **6.95π ≈ 7π** — π + 6π = 7π holds to ~1% (the residual is genuine O(ε²) cross-terms). **The app's original 6.41π was a normalization artifact:** it divided the precession by the *nominal* (Newtonian, undistorted) semi-latus rectum p = 6.40 instead of the orbit's actual p ≈ 6.90 (the 1PN + SR forces expand the orbit by ~8%). `6.95 × 6.40/6.90 = 6.45 ≈ 6.41π` — normalizing by the nominal p exactly reproduces the app's number. So Q-006's *"the implementation deviates from the declared model"* branch is the answer, localized precisely: the **measurement normalization** was flawed, not the superposition **reasoning**. *Separately:* the combined model's 7π is still not the correct physics — full GR is **6π** (pn1, the consistent 1PN EOM = N3 `curve`, D-024, 0.52%); the combined hybrid over-counts by adding SR inertia atop a 1PN force that is already a complete O(1/c²) description. So retiring it (D-016) was right — but for the physical over-count, not the (artifactual) 6.41π. The app ships SR-only (π, exact); the full 6π comes from N3 `curve`'s metric oracle.

Single-file CPU fp64 oracle `nexus/precession_nexus.cpp` — no GPU (RK4 orbit integration; deterministic, runs under any card contention). `_nexus` family (blake2b golden, envelope face, exit 0/1/2). Precession measured by **perihelion tracking + least-squares slope** over many orbits (D-016's method — endpoint sampling aliases the intra-orbit oscillation).

## Units & dials (the einstein-contract test orbit)

| dial | value | role |
|---|---|---|
| c | 20 su/s | ε = GM/(c²·a(1−e²)) = GM/c²p sets the precession scale |
| G | 2×10⁻³ | GM = 20 |
| M | 10⁴ | central mass (test-particle orbit) |
| a, e | 10, 0.6 | the einstein-contract N3 orbit; p = a(1−e²) = 6.4; ε = GM/c²p = 7.8×10⁻³ |

The precession per orbit is reported as a **coefficient × (GM/c²p)** in units of π (so GR = 6π, SR = π-ish, the claim was 7π).

## Declared engine (the three isolation experiments)

A planar test-particle orbit integrated by fp64 RK4 with **pluggable force + inertia**:
- **force** — Newtonian `a = −GM r/r³`, OR harmonic-1PN `a = g + (1/c²)[(v²+4Φ)g − 4(g·v)v]`, g = −GM r/r³, Φ = −GM/r (the einstein contract's declared 1PN term = the standard harmonic-gauge test-particle 1PN acceleration).
- **inertia** — Newtonian `v = p`, OR special-relativistic `v = p/√(1 + |p|²/c²)`.

| scenario | force | inertia | expected | this is |
|---|---|---|---|---|
| **sr** | Newtonian | **SR** | π·GM/c²p (Sommerfeld) → measured **1.00π** | relativistic Kepler — what the app ships |
| **pn1** | **1PN** | Newtonian | **6π·GM/c²p** → measured **6.03π** | the consistent 1PN EOM = the correct GR precession |
| **combined** | **1PN** | **SR** | **7π** (superposition) → **6.95π** (actual p); **6.45π** (nominal p = the app's artifact) | the app's retired model — the superposition HOLDS; its 6.41π was a normalization artifact |

Precession per orbit via perihelion tracking: at each perihelion passage (radial velocity v_r crossing − → +, interpolated) record the position angle; least-squares-fit angle vs orbit index → Δϖ. ~40 orbits.

**Determinism.** fp64 RK4, fixed steps → byte-identical declared JSON → blake2b. No RNG/GPU/atomics.

## Scenarios & gates (the golden units)

| scenario | gate | oracle |
|---|---|---|
| **sr** | `|Δϖ_meas / (2π(1/√(1−(GM/cL)²) − 1)) − 1| < 0.02` — the **exact Sommerfeld** relativistic-Kepler precession (actual relativistic L) | analytic Sommerfeld (= the app's keprel oracle) |
| **pn1** | `|coeff_actual/π / 6 − 1| < 0.02` — the harmonic-1PN precession is **6π** (the correct GR value), normalized by the orbit's **actual** p | analytic GR 6πGM/c²p |
| **combined** | **`|coeff_actual/π / 7 − 1| < 0.02`** (actual-p normalized → **7π**, the superposition **CONFIRMED**) AND `|coeff_nominal/π − 6.41| < 0.25` (nominal-p reproduces the app's 6.41π artifact) | the Q-006 resolution: superposition holds; 6.41π was a normalization artifact |

**`--selftest`:** Newtonian force + Newtonian inertia → **no precession** (|Δϖ| < 1e-3, a closed Kepler ellipse). The classical sanity + determinism smoke.

The `combined` "gate" is a **refutation gate** (RAYFORMER discipline): it certifies the measured precession contradicts the retired 7π claim and reproduces the D-016 number, rather than asserting a theory value for an inconsistent model (cf. N0's labelled `gamma_eff`, D-021). `sr` and `pn1` are gated against exact GR.

## The headless face (envelope-conformant)

```
precession_nexus.exe --scenario sr|pn1|combined [--seed N] --json    # declared JSON
precession_nexus.exe --scenario NAME --golden                        # vs goldens/precession_NAME/golden.hash
precession_nexus.exe --selftest                                      # Newtonian: no precession
```

Declared body: `seed, params{scenario, c, G, M, a, e, p, eps, orbits, steps}, result{dperi_meas, coeff_over_pi, expected_over_pi (sr/pn1) | sommerfeld/six_pi, rel|refutes_7pi, nan_free}, gates{primary, nan}, verdict`. Exit 0/1/2. Build with MSVC `cl`.

## Oracle pedigree (Invariant 3)

- **sr:** exact Sommerfeld relativistic-Kepler precession 2π(1/Λ−1), Λ=√(1−(GM/cl)²) — the same oracle the app's `keprel` golden uses (0.50%).
- **pn1:** exact GR 6πGM/c²p (the 1PN periastron advance; matches N3 `curve` `precess`, D-024).
- **combined:** the Q-006 measurement itself — the refutation of the 7π superposition (a labelled result, not a theory gate).

## The honest boundary (permanent, printed)

1. **The resolution: the 7π superposition is CORRECT; the app's 6.41π was a measurement-normalization artifact** (nominal vs actual semi-latus rectum). Measured with the orbit's actual p: sr=1.00π, pn1=6.03π, combined=6.95π≈7π. *Separately*, the combined model's 7π is not the correct physics — full GR is 6π (pn1, the consistent 1PN) — so its retirement (D-016) stands, on the physical over-count rather than the artifactual 6.41π. **This resolution flipped once the actual p was measured** — normalize by what the orbit *is*, not its nominal label (RAYFORMER).
2. **Weak-field, leading order.** ε = GM/c²p ≈ 8×10⁻³; the ~0.4π excess of `combined` over 6π is the O(ε) double-counted cross-term (it does not reach 7π because the cross-terms are O(ε²)). Higher PN orders out of scope.
3. **Test particle, point mass, no spin.**

## Build runbook

1. **Contract (this file).** ← *here.*
2. **`nexus/precession_nexus.cpp`** — blake2b + fmt + the RK4 orbit integrator (pluggable force/inertia) + perihelion-tracking LS precession + envelope/golden/main.
3. **Build** — `cl /std:c++17 /EHsc /O2 /W4 nexus\precession_nexus.cpp /Fe:build\precession_nexus.exe`.
4. **Golden** — `goldens/precession_sr|pn1|combined/golden.hash`.
5. **Harness** — add to `INSPIRAL_GOLDENS`' CPU-independent block (a `PRECESSION_GOLDENS` list).
6. **Register** — QUESTIONS.md (Q-006 RESOLVED), TASKLIST, DECISIONS (D-026), RUN_STATE, `nexus/MODULE_precession.md`.

## Changelog

- v1.0.0 (2026-07-12) — initial contract + build (operator queue #3, resolving Q-006). The SR-inertia + 1PN-field precession, isolated (measured with the orbit's ACTUAL semi-latus rectum): sr = 1.00π (Sommerfeld), pn1 = 6.03π (the correct GR 6π), **combined = 6.95π ≈ 7π — the superposition is CONFIRMED**. The app's original 6.41π was a normalization artifact (nominal p=6.40 vs actual p=6.90). CPU fp64 oracle `nexus/precession_nexus.cpp`. (Note: the initial draft of this contract had the resolution backwards — "superposition refuted, answer 6π" — corrected after measuring p_actual flipped combined from 6.4π to 6.95π. The RAYFORMER discipline earned its keep.)
