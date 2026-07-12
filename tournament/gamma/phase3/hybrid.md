# PHASE 3 — the ruling (rule every fork with a reason; never average)

**Operator directive:** *pick the most promising, go ambitious.* **Mode:** minimal-delta synthesis off the incumbent `substrate_nexus`; every conflict ruled with a reason; graveyard respected; every deletion carries a reinstatement trigger.

## The ruling in one sentence

**Resolve the scalar γ ≈ 0.374 via the CSS × spectral *similarity-collocation* path — the eigenvalue/perturbation route (γ = 1/Re λ₀) that the resolvability *and* oracle-honesty lenses both ranked #1 — validated first on a known critical exponent (the fluid-CSS warm-up oracle), with `reframe` (boson-star) as the pre-registered floor and AMR's γ deferred to N3.**

## D-021, amended (a real correction the tournament earned)

D-021 concluded "a uniform grid cannot resolve γ." **Evidence-grade correction (phase 2, `cheb_convergence` reproduced `b6ad2eba`): D-021 is a *gauge × discretization* artifact, not a grid law.** 2nd-order FD converges at order −1.977 and needs 9,722× the DOF of Chebyshev; polar-areal singularity-avoiding slicing crowds resolution into α(0)→0 exactly where the self-similar structure concentrates. What D-021 correctly killed is *uniform-2nd-order-FD-polar-areal*, not "a fixed deterministic grid." The honest residual it also proved: **γ<0.05 has historically required either refinement OR a coordinate change making the self-similar solution stationary/smooth-on-grid** — which is exactly the ruled path. (→ amends into `DECISIONS.md` on merge.)

## The forks, ruled

- **F1 formulation — RULED: similarity (CSS/DSS) coordinates.** Reason: they make the self-similar solution *stationary/periodic in τ*, so a fixed deterministic grid resolves it (kills D-021's mechanism by deletion, not refinement). Polar-areal (incumbent) and double-null both *evolve* a shrinking feature; similarity coordinates stop it shrinking.
- **F2 resolution — RULED: spectral (Chebyshev) collocation on the fixed similarity grid.** Reason: the *only* approach carrying **measured** evidence it escapes the FD ceiling (9,722× DOF, verified), and geometric convergence is what buys γ<0.05 at CPU budget. Deterministic by construction (no adaptivity, fixed differentiation matrices).
- **F3 observable — RULED: the eigenvalue γ = 1/Re λ₀ (primary), with mass-scaling / echo Δ as redundant-recovery cross-checks.** Reason: an eigenvalue *must converge* and is immune to grid-capped-curvature fit-noise; it is the honesty winner. Redundant recovery (≥2 observables agree) is a firing gate.
- **F4 crown — RULED: go for scalar γ (ambitious), `reframe` as the pre-registered floor.** Reason: operator directive + the resolvability verdict that γ is a *gauge/discretization* problem, not impossible. Reframe (boson-star M_max=0.633 to <1%, deterministic) is the reversibility fallback: **if the γ tool busts the ≤5-min budget or the DSS construction proves intractable in-scope, reframe ships as N0's crown and scalar-γ defers to N3** — nothing wasted.
- **F5 search — RULED: `autotune` level-crossing against a pre-registered target (deterministic, hash-cited), not hand-bisection.** Reason: proven evidence-grade across all five personas; kills the p\*-wander half of D-021.

## Graveyard (re-proposal needs a genuinely new argument)

1. **uniform-2nd-order-FD-polar-areal for γ** — D-021, now understood as gauge×discretization (not a grid law). Reinstatement: never for γ; it remains correct for the *transition* (N0 v1.0.0 ships it).
2. **double-null for a tight γ** — WOUNDED 2/3 lenses: its γ<0.05 overreaches its own exhibit (Garfinkle reported echoing, **no exponent**; the one null-γ-accuracy precedent needed *characteristic* AMR), and it carries no firing metamorphic gate ("the exact gap D-021 fell through"). Reinstatement trigger: a characteristic-AMR null contract at N3. *Steal kept: null coordinates for the light-history / horizon work.*
3. **AMR-in-N0** — determinism *cured* (constraint-solved polar-areal → no refluxing) but the γ golden is 6–20 min, busting ≤5 min. **Reinstatement: N3 (GPU), where AMR is affordable** — it is the published reference method and belongs there.

## Steals adopted (house discipline for the γ tool)

- **`G-CONVERGE`** — γ (the eigenvalue) MUST be stationary under grid refinement or the tool **exits 1**. The D-021 tooth made mechanical; non-negotiable.
- **`G-UNIQUE`** — exactly one growing perturbation mode (the DSS physics-honesty gate; Gundlach's single unstable mode).
- **Redundant recovery** — eigenvalue γ and mass-scaling γ (and/or echo Δ=3.44) must agree to < 0.03, else exit 1.
- **Warm-up oracle** — validate the eigenvalue machinery on a *known* critical exponent before claiming the scalar γ (the continuously-self-similar fluid case, β≈0.3558 — Koike–Hara–Adachi — a clean ODE eigenvalue; or an equivalently-anchored warm-up).
- **Hash domain = `declared_object`** (D-013: notes excluded, `%.6f`), not raw stdout, before any golden freezes.

## The chosen approach — pre-registered deciding experiment

**Tool** `similarity_nexus` (working name): single-file C++17 fp64, spectral collocation in similarity coordinates, constructs the self-similar critical solution + its linear perturbation spectrum, extracts γ = 1/Re λ₀. `autotune` drives the critical search. Deterministic golden over `declared_object`; ≤ 5 min CPU.

**PASS ladder (staged, ambitious, honest):**
1. **Warm-up (validates the loop on ground truth):** eigenvalue machinery recovers a *known* critical exponent to < 3% (fluid-CSS β≈0.3558, or an equivalently-anchored CSS test) — with `G-CONVERGE` + `G-UNIQUE` firing correctly. *This alone validates the tournament → ORRERY → build → measure loop on a real exponent.*
2. **The crown (ambitious target):** scalar DSS γ recovered with **|γ − 0.374| < 0.05**, `G-CONVERGE` (stationary under refinement), `G-UNIQUE` (one growing mode), redundant recovery (< 0.03), deterministic golden, ≤ 5 min.
3. **Floor (reversibility):** if (2) busts budget or the DSS construction is out-of-scope, `reframe`/`bosonstar_nexus` ships as N0's crown (M_max to <1%) and scalar-γ defers to N3 — pre-registered, not a failure.

**Honesty stance (D-016/D-021):** the scalar DSS γ is genuine NR frontier work on a CPU; the build reports **what actually lands** — a validated warm-up + real progress is a real result; a faked 0.374 is not. The register holds the doubt.
