# CONTRACT — substrate_nexus (v2 N0: the substrate oracle) · v1.0.0 · status: FROZEN 2026-07-11 (operator directive "go substrate")

**Purpose.** The foundation stone of TINY UNIVERSE v2 "SUBSTRATE" (`docs/PROPOSAL_2026-07-12_v2_substrate.md`): a CPU fp64 spherically-symmetric Einstein–Klein–Gordon solver whose crown result is the **Choptuik critical-collapse mass-scaling exponent γ ≈ 0.374** — a number from the deepest end of gravitational physics, frozen as a golden **before any GPU substrate code exists**. It is the v2 analogue of M0 `tiny_nexus`: the standing fp64 oracle every future GPU substrate module (N1 field · N2 lapse · N3 horizon · N4 star) is gated against. Single file, C++17, fp64, stdlib only, no GPU (starts even while the shared card is busy elsewhere — the proposal's whole point).

## Units (declared, and why they differ from v1)

- **Geometric G = c = 1.** The Choptuik exponent is a *dimensionless, universal* number — it does not depend on the unit system. So N0 works in standard numerical-relativity units where the literature value γ ≈ 0.374 applies directly (Choptuik 1993; Gundlach γ = 0.3737). The v1 dials (c=20, ħ=0.5, G=2×10⁻³) are what make the *GPU* substrate (N1–N4) computable in-box — the proposal's thesis — not what this oracle needs. N0 owns no dials; it owns the equations.
- ħ does not appear (N0 is classical field + GR; the scalar is a classical KG field). The KG mass μ is a per-scenario parameter (μ=0 for Choptuik/massless; μ>0 for the massive-field bound test).

## Declared formulation (spherical polar-areal EKG, constrained evolution)

Metric `ds² = −α(t,r)²dt² + a(t,r)²dr² + r²dΩ²`. Real scalar φ, potential `V(φ) = ½μ²φ²`. First-order field variables `Φ ≡ ∂_rφ`, `Π ≡ (a/α)∂_tφ`.

**Evolution (method of lines; only the field evolves in time):**
- `∂_tφ = (α/a)Π`
- `∂_tΦ = ∂_r[(α/a)Π]`
- `∂_tΠ = (1/r²)∂_r[r²(α/a)Φ] − α a μ²φ`

**Constraints (a, α solved by outward radial integration each substep — the elegance of spherical symmetry):**
- `(ln a)_{,r} = (1−a²)/(2r) + 2πr(Π²+Φ²) + 4πr a² V`, with `a(t,0)=1`
- `(ln α)_{,r} = (a²−1)/(2r) + 2πr(Π²+Φ²) − 4πr a² V`, integrated from `α(0)=1` then rescaled so `α·a → 1` at `r_max` (proper time = coordinate time at infinity)

**Mass aspect** `m(t,r) = (r/2)(1 − 1/a²)`; ADM mass = `m(t,r_max)`. **Apparent horizon / black hole**: `2m/r → 1` (equivalently `a` diverges; polar slicing is singularity-avoiding, so collapse shows as `max_r(2m/r) → 1` together with **central lapse collapse** `α(0) → 0`).

**Numerics (declared, golden-pinned):** uniform radial grid `r_j = j·dr`, `j=0..N`; RK4 in time; centered 2nd-order finite differences; Kreiss–Oliger dissipation (ε_KO, declared) for stability; origin regularity by parity ghost zones (φ even, Φ odd, Π even); outgoing (Sommerfeld) radiation at `r_max`. CFL `dt = λ·dr` (λ declared). Fully deterministic fp64, single-threaded ⇒ `(params) → byte-identical declared JSON → blake2b` (tiny_nexus idiom; N-determinism gate runs the battery twice in-process).

## Battery (S1–S5; the crown is S4)

| id | scenario | gate | oracle |
|---|---|---|---|
| **S1** | flat-space wave — ε-amplitude Gaussian pulse (p=10⁻³), gravity negligible | ADM mass conserved `|M(t)−M(0)|/M < 1e-3`; pulse disperses (central \|φ\|→~0); `α(0)`, `a_max` stay ≈1; NaN-free | exact massless spherical wave + energy conservation |
| **S2** | subcritical dispersal — Gaussian below p* | `max_r,t(2m/r) < 1` (no horizon); `α(0)` dips then recovers `> 0.9`; field radiates off-grid; mass conserved to grid truncation | one side of the critical point (field disperses) |
| **S3** | supercritical collapse — Gaussian above p* | apparent horizon forms: `max_r(2m/r) > 0.99`; central lapse collapses `α(0) < 1e-2`; horizon mass `> 0` and `< M_ADM` (cosmic censorship sanity) | other side (black hole forms); Schwarzschild `M = r_AH/2` |
| **S4** | **Choptuik critical transition** — bisect p* (time-symmetric Gaussian family); verify the threshold brackets cleanly (below disperses, above collapses); measure `M_BH(p)` for a supercritical sequence and confirm the **Type-II signature**: `M_BH` falls continuously toward zero as `p → p*` (arbitrarily small black holes) | threshold bracketed · `M_BH` monotone-rising with p · `M_BH,min / M_BH,max < 0.35` (mass → 0 at p*) · a clean power-law fit exists (effective exponent + R² **reported as a grid-limited diagnostic, NOT gated to 0.374**) | the Choptuik *phenomenon* (Type-II transition); the precise universal γ ≈ 0.3737 is **deferred to the AMR contract** (D-021) |
| **S5** | massive-field stability — μ>0 Gaussian | massive-KG term (V in both field evolution and constraints) evolves **stably** (finite) and **conserves ADM mass** `< 2e-2`; no spurious collapse | verifies the mass term is correctly implemented; a bound oscillaton eigen-profile is **deferred to a later contract** (not claimed) |

Battery green ⇒ golden frozen (`goldens/substrate_nexus/golden.hash`), exit 0. A fired gate ⇒ exit 1. Error ⇒ exit 2 (never conflated). `--only S1..S5` runs one scenario; `--dev P` runs one evolution with diagnostics; `--json` prints the declared envelope; `--golden` freezes/checks the hash. Frozen config: `N=800, r_max=24` (~20 s).

## The crown, honestly (D-021 — RAYFORMER discipline)

The proposal named the precise Choptuik exponent γ ≈ 0.374 as N0's crown. **Measurement (this build) shows a uniform grid cannot resolve it**: critical collapse is self-similar, its central curvature growing without bound as p → p*, and a fixed grid caps that curvature (Choptuik needed adaptive mesh refinement). At the frozen resolution the *far-field* mass scaling is clean and reproducible (R² ≈ 0.998) but yields an effective exponent (≈ 0.24, grid-limited), **not** 0.374. Per RAYFORMER/D-016, N0 therefore ships what it can honestly prove — the **Type-II critical transition itself** (a sharp p* with M_BH → 0), the rigorous flat/subcritical/supercritical/massive oracles, and the effective exponent as a labelled diagnostic — and **defers the precise universal γ to the AMR contract** (N-later / the S3 `horizon` GPU stage may also reproduce it). Faking a tight γ gate against grid noise would poison every downstream v2 module gated on this oracle. The proposal (§6) anticipated exactly this: "fixed-grid horizon resolution bounds the smallest honest BH; declare it; AMR is a later contract."

## Gate (N0)

The spherical EKG evolution is trustworthy (flat-space mass conservation, subcritical dispersal, supercritical horizon formation, massive-KG stability — all green oracles a stranger runs cold), and the **Choptuik phenomenon is demonstrated** (Type-II transition, black holes forming with mass → 0 at a bracketed critical amplitude) — frozen as a golden, **before any GPU substrate code**. Precise γ: deferred, declared (D-021).

## Changelog

- v1.0.0 (2026-07-11) — initial freeze under operator directive "go substrate". CPU-only (the foundation stone by design). S4 crown rescoped at implementation to the Type-II *transition* (precise γ deferred to AMR — D-021, RAYFORMER/D-016 discipline); S5 rescoped to massive-KG stability+conservation (oscillaton eigen-profile deferred). Frozen config N=800/r_max=24, golden `13aa73e5`.
