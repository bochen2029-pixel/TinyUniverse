# CONTRACT — similarity_nexus (scalar-field DSS critical exponent γ, eigenvalue route) · v0.1.0 · status: **OPERATOR-APPROVED 2026-07-19 ("contract approved, proceed with the gamma crown") — BUILDING** (research scaffold first, C++ after Python validation, per the proven fluid pipeline)

**Purpose.** Measure the **Choptuik critical exponent of massless-scalar gravitational collapse** —
γ ≈ 0.374 — and the **echoing period** Δ ≈ 3.4453, by the eigenvalue route (D-028: γ = −1/Re λ₁ from
the discretely-self-similar critical solution's one growing perturbation mode), NOT from mass-scaling
fits near a wandering p\*. This is the γ-tournament's ruled path (CSS/DSS similarity coordinates ×
spectral collocation × eigenvalue observable), built ONLY NOW because the machine was first **proven
end-to-end on a known exponent**: the fluid-β warm-up (fluidcss_nexus v1.0.0, D-032 — β = 0.3557988
vs lit 0.35580192, with two analytic gauge controls). Primary source: **Gundlach gr-qc/9604019**
(`tournament/gamma/GUN96_src/chop2.tex`, on disk, gitignored; re-fetch:
`python C:\fetcher\fetch.py url "https://arxiv.org/e-print/gr-qc/9604019" --dest tournament\gamma\GUN96_src`).

**Relation to N0.** `substrate_nexus` (v2 N0) demonstrates the Type-II *transition* and DEFERS precise
γ (D-021, amended D-028). This tool supplies the deferred number by the coordinate route that deletes
D-021's mechanism (the self-similar solution becomes *stationary/periodic* on a fixed grid).

## The object (two stages, mirroring the proven fluid template)

**Stage A — the DSS critical background (Gundlach §2).** Fields Z = (a, g, X₊, X₋) on the cylinder
(ζ, τ) with τ-period Δ; coordinates τ = ln(t/r₀), ζ = ln(r/t) − ξ₀(τ) with the gauge function ξ₀(τ)
chosen so the **past self-similarity horizon (SSH) sits at ζ = 0 for all τ** (eq. coordcond — the DSS
analog of the fluid tool's sonic normalization). ζ-evolution (space-as-time), pseudo-Fourier in τ
(N_τ harmonics; the solution is smooth ⇒ rapid convergence). BCs: regular center (a = g = 1 at
ζ = −∞) + SSH analyticity at ζ = 0 (numerator regularity condition — Gundlach shows the single
condition enforces full analyticity, the exact analog of the fluid sonic condition). Reflection
symmetry imposed (κ = 0 boundedness of φ): a, g, ξ₀ even harmonics only; X± odd only. **Δ and ξ₀(τ)
are eigenvalues/eigenfunctions of the nonlinear BVP.** a is not evolved: reconstructed from the
constraint each ζ-step (Gundlach's own discipline).

**Stage B — the Floquet perturbation eigenvalue (Gundlach §4).** Ansatz δZ = e^{λτ} δᵢZ(ζ,τ) with
δᵢZ Δ-periodic; the linearized four-field system with the shifted operator (A + λB); BCs: one free
center function + SSH numerator-vanishing (+ linearized constraint) leaving free SSH data; two-sided
ζ-shoot matched at an intermediate ζ; λ = the eigenvalue. **In Gundlach's convention (t > 0, modes
grow toward the singularity τ = −∞): the relevant mode is Re λ < 0; λ₁ = −2.674 ± 0.009 ⇒
γ = −1/λ₁ = 0.374 ± 0.001.** Search strip: Im λ ∈ [0, π/Δ) (harmonic aliasing absorbs the rest).
Even/odd symmetry sectors decouple — search both, declare per-sector counts.

## Analytic controls (the D-032 lesson institutionalized)

1. **Gauge mode at λ = 0 exactly** (τ-translation of the background, δZ ∝ Z\*,τ — Gundlach reports
   it): the eigenvalue reader must resolve it to |λ_gauge| < gate. The direct analog of the fluid
   tool's origin-gauge control at κ = 1.
2. **The Δ/2 impostor check ("WHICH solution did you find")**: the Choptuik solution's φ → −φ
   half-period symmetry means a wrong-symmetry ansatz can land Δ/2 ≈ 1.72 or a φ-unbounded (κ ≠ 0)
   sibling. Declared: the κ-residual (must vanish), the harmonic-parity residuals, and Δ vs BOTH
   3.4453 and 2×1.7226 (must match the former).
3. **Redundant recovery across machines (D-028 gate)**: Δ from this tool vs the echoing period
   measured in `substrate_nexus` S4 near-critical evolutions (the N0 oracle farm) — agreement < 0.03.
   (γ itself has no in-house second machine until N4-GPU/AMR; its anchor is the literature value.)

## Declared state (hash domain; envelope idiom of the `_nexus` family)

Stage A: {N_τ, dζ, iteration counts; Δ; ‖ξ₀‖ harmonics; SSH-residual; center-residual; κ-residual;
parity residuals; Δ-convergence spread (N_τ↑, dζ/2)}.
Stage B: {λ₁ (Re, Im); γ = −1/Re λ₁; λ_gauge; per-sector growing-mode counts in the strip;
convergence spread; matching residual}. `%.9e` canonical serialization, notes excluded, blake2b-256.

## Gates (numbers PROPOSED — operator may tighten/loosen at review)

- **G-DELTA**: |Δ − 3.4453| < 5e-3 (lit ±0.0005; leave headroom for truncation systematics — see
  Honesty below). **G-DELTA-N0**: |Δ − Δ_N0-echo| < 0.03 (D-028 redundant recovery).
- **G-ANCHOR**: |γ − 0.374| < 4e-3 (mirrors the fluid G-ANCHOR; lit ±0.001).
- **G-GAUGE**: |λ_gauge| < 1e-3.
- **G-UNIQUE**: exactly ONE Re λ < 0 mode across both symmetry sectors in the strip.
- **G-CONVERGE**: Δ and λ₁ stationary under N_τ → N_τ+4 and dζ → dζ/2 (spread < 1e-3 each).
- **G-PARITY/G-KAPPA**: parity + κ residuals < 1e-10 (structural exactness of the symmetry ansatz).

## Determinism

Fixed N_τ, fixed ζ-step RK4, fixed Newton/bisection iteration counts, fp64 throughout, no
fast-math; single-file `substrate/similarity_nexus.cpp`, faces `--stageA | --stageB | --json |
--golden | --selftest`, exit 0/1/2. Goldens `similarity_stageA` / `similarity_stageB` frozen only
when every gate fires (D-016: no β/γ-bearing golden before the gates).

## Honest boundaries (named now, per D-016/D-021/D-032)

- **Truncation systematics are real**: Gundlach's own Δ moved 3.5σ (3.4439 → 3.4453) between his PRL
  and the long paper from a method inconsistency — the G-CONVERGE spread is a *floor*, not the full
  systematic; declare both.
- The Floquet strip search is bounded (strip + Re λ window [−10, 0)); modes outside are out of scope.
- The maximal-extension/curvature sections of the paper (§3, §Curvature) are NOT in scope — the tool
  measures Δ and γ only.
- If the nonlinear Stage-A BVP proves out-of-budget, the pre-registered fallback (D-028 ③) stands:
  the tournament's `reframe`/boson-star floor ships as N0's crown and γ defers to N4-GPU — nothing
  faked.

## Sanctioned prior art to port (verbatim-with-attribution, per house rules)

The fluid machine (proven, D-032): envelope/blake2b idiom, deterministic bisection shape, two-sided
shoot + analytic-control pattern (`substrate/fluidcss_nexus.cpp`); the Python scaffold pattern
(`tournament/gamma/phase4/nr_ec2.py`/`nr_lyap.py`/`nr_shoot_ec.py`) for the research-side validation
run that precedes the C++ port. Gundlach's appendices A–F are the implementation spec (BCs,
pseudo-Fourier operators, numerics, error estimates).

## Build

`cl /std:c++17 /EHsc /O2 /W4 substrate\similarity_nexus.cpp /Fe:build\similarity_nexus.exe` (CPU
fp64, no GPU — runs under any card contention, like N0/curve/fluidcss).

## Changelog

- v0.1.0 (2026-07-19) — **DRAFT for operator review** (autonomous session, post-D-032). Scope, both
  stages, analytic controls (λ=0 gauge; Δ/2-impostor check; N0-echo redundant recovery), gates,
  determinism, honest boundaries, fallback. **No source until this contract is operator-approved**
  (hard rule). Primary source fetched + digested (§2.2–2.3, §4.1 read verbatim; appendices A–F
  identified as the implementation spec).
