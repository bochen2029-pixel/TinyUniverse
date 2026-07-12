# Stage B UPDATE — the perturbation operator is now CORRECT (2026-07-12)

Follow-up to `RESULTS_hka_beta.md`. The Stage-B blocker (the perturbation operator failing the
gauge-mode exactness gate) is **RESOLVED**. Isolating the physical eigenvalue κ is the remaining step.

## The fix (`hka_pert_derive.py`) — the breakthrough
Instead of transcribing HKA's precomputed coefficients (5.5)–(5.13) [subtly wrong as transcribed], the
operator is **DERIVED by directly linearizing the full background PDEs (3.6a,b,d,e)** with ∂/∂s → κ:
`A→A(1+ε·aa), N→N(1+ε·nn), om→om(1+ε·oo), V→V+ε·vv`, mode `~ ε·pert(x)·e^{κs}`.
Correct by construction: the gauge transformation is an EXACT symmetry of (3.6), so the pure-gauge mode
(5.20) is annihilated automatically. `hka_pert_derive.py` overwrites `hka_pert_L.pkl` (drop-in for
`hka_pert_core.Lnum`; run it once to regenerate the operator).

## Verification — ALL PASS
- **GAUGE-MODE GATE (decisive):** residuals **~1e-6…1e-9, INVARIANT across κ̄=0.357/1.0/2.81** (was
  O(1)–O(1e4) and κ̄-dependent). The gauge mode is now an exact solution for all κ̄, as required.
- **Ā row** reproduces the known-correct (G1,0,G3,G4) exactly.
- **Center indicial** eig(L)={−2,−1,0,0} (κ-independent): 2 regular (λ=0) + 2 excluded; λ=0 geometric
  multiplicity **2** (non-defective) — the 2 regular modes = gauge + physical.
- **Sonic residue** eig(R)={0,0,0,1−2κ} **exactly** (3 analytic + 1 non-analytic), matching HKA.

## Eigenvalue extraction — IN PROGRESS (not yet landed; NO β claimed)
Two-sided shoot (`hka_beta_solve.py`): center 2D regular subspace integrated to x_mid; sonic 3D analytic
Frobenius modes integrated outward to x_mid; eigenvalue = rank-drop of the stacked [U|W] 4×5 matrix
(σ₄→0, i.e. the physical mode joining the always-present gauge mode in the center-regular ∩ sonic-analytic
intersection). **STATUS:** σ₄ shows structure but no crisp, refinement-stable zero at κ≈2.81 — the
strongly-singular non-analytic mode t^{1−2κ} (exponent −4.6 at κ=2.81) makes the plain shooting match
stiff (σ₄ scales with integration refinement at all κ). **No β is reported until a refinement-stable zero
is isolated. No number faked.**

## Next (recommended)
Robust eigenvalue isolation on the (now-correct) operator via a method that handles the singular sonic
point cleanly: a **global Chebyshev/collocation discretization** of Ψ′=L(x;κ)Ψ on [x_center, x_sonic⁻]
with regular-center + analytic-sonic BCs as a nonlinear-in-κ generalized eigenproblem, **or** a
compound-matrix / Riccati two-sided determinant (immune to the exponential mode-amplification). Then
β=1/Re(κ), discard the gauge modes κ≈0.357/1 (HKA fn15), port to `fluidcss_nexus.cpp`, freeze the golden.

Files: `hka_pert_derive.py` (the correct operator), `hka_beta_solve.py` / `hka_beta_probe.py` (WIP shoot),
`hka_pert.py` (the gate — now passes with the derived operator).
