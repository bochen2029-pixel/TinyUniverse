# CONTRACT — fluidcss_nexus (radiation-fluid CSS critical exponent, eigenvalue route) · v0.9.1 · status: **Stage-A golden FROZEN** (Evans–Coleman background, `goldens/fluidcss_stageA/`, D-027); Stage-B/β an honest wall — see RESULTS

**Purpose.** A single-file fp64 CPU tool that measures the **critical exponent of radiation-fluid
gravitational collapse** `β = 1/Re κ₀` via the **eigenvalue route** of Koike–Hara–Adachi
(PRL 74 (1995) 5170, gr-qc/9503007): construct the continuously-self-similar (CSS) Evans–Coleman
critical solution as a shooting problem through the sonic point, then solve the linear perturbation
eigenvalue problem for the single relevant mode. The verified target is
**Re κ₀ = 2.81055255, β = 0.35580192** (KHA95 §5), with an independent evolution cross-check
β ≈ 0.36 (Evans–Coleman, gr-qc/9402041). This is the warm-up **oracle** of the tournament→build→
measure loop: landing this known number honestly validates the whole machine on ground truth.

**Status note (honesty, D-016/D-021).** As of v0.9.0 the load-bearing **Stage-A sonic-point regularity
condition is derived and independently verified three ways** (§ below). The **Stage-A background
construction and Stage-B eigenvalue shoot are NOT yet landing β**, because reconciling the primary-source
ODE system exposed a discrepancy (KHA's printed fluid equations vs. an independent covariant derivation)
that is diagnosed but not resolved — see `tournament/gamma/phase4/RESULTS_fluidcss.md`. Per the honesty
mandate, **β is not reported as measured** until Stage B converges with all gates firing. This contract
freezes the verified Stage-A object and the target/gate specification for the completed tool.

## Units & conventions (declared)

- Geometric `G = c = 1` (β is a dimensionless universal number; no dials).
- Equation of state `p = ρ/3` (radiation, adiabatic index γ_ad = 4/3).
- KHA similarity variables `s ≡ −ln(−t)`, `x ≡ ln(−r/t)` (t < 0); reduced fields
  `N ≡ α a⁻¹ e^{−x}`, `A ≡ a²`, `ω ≡ 4π r² a² ρ`, `V` = fluid 3-velocity.
- Sonic point gauge-fixed to `x = 0`. Center at `x → −∞` (BC: `A = 1, V = 0`). Dispersal at `x → +∞`.
- Perturbation ansatz `h(s,x) = H_ss(x) + ε h_p(x) e^{κs}`, `κ ∈ ℂ`; gauge `N_p(s,0) = 0`.
- `β = 1/Re κ₀` (KHA eq. 15), `κ₀` the unique relevant (`Re κ > 0`) mode.

## Verified Stage-A physics (the §1.6 derivation — DONE)

**Sonic locus (det = 0 condition).** Setting `∂_s → 0` in KHA eq. (18) reduces the fluid pair to a
2×2 linear system for `(ω_,x, V_,x)`. Its determinant vanishes on the sonic locus

> **`3 N² V² − N² + 4 N V − V² + 3 = 0`**  (equivalently `3 − V² − N² + 4NV + 3N²V² = 0`).

This is independently confirmed:
1. Symbolic elimination of `(A_,x, N_,x)` from the transcribed 2×2 reproduces it **exactly** (up to a
   harmless 4/3 factor) — matching the operator's pre-derived condition.
2. The sound-speed identity: `[(V+c_s)/(1+V c_s) + N]·[(V−c_s)/(1−V c_s) + N]` with `c_s = 1/√3`
   shares the same physical zero locus (product-of-two-sound-ray-conditions structure).
3. The two KHA **metric** slope equations `A_,x/A` and `N_,x/N` are verified term-by-term against the
   **independent** Evans–Coleman field equations (4),(5) (gr-qc/9402041) — exact algebraic identity.

**Regularity across the sonic point.** At `det = 0` the system is a Fredholm-singular point; a locally
analytic solution requires the Cramer numerators to vanish there (removable singularity), giving the
local expansion, with the single free parameter `V_ss(0)` fixed to a discrete set by the center BC.

## Declared formulation (the completed tool)

**Stage A — background CSS critical solution.** Integrate the autonomous ODE system (KHA eq. 18,
`∂_s = 0`): `A_,x, N_,x` from the constraint pair; `(ω_,x, V_,x)` from the 2×2 solve. Shoot the sonic
parameter `V_ss(0)` (equivalently the center density) so the solution is analytic through `x = 0` and
regular at the center (`A → 1, V → 0` as `x → −∞`), reaching the dispersive state (`V → 1, ω → 0`) as
`x → +∞`. Oracle anchors (Evans–Coleman): similarity exponent `n ≈ 1.1485`; near-field `a ≈ 1.07`
(`A ≈ 1.145`), `Ω ≈ 9.56×10⁻³`, `m/r → 0.0596`.

**Stage B — perturbation eigenvalue.** Linearize eq. (18) about `H_ss(x)` with `h_p(x) e^{κs}`;
homogeneous first-order linear ODEs for `(N_p, A_p, ω_p, V_p)`. Regularity at the center + analyticity
through the sonic point (now a regular singular point) leave a discrete `κ` set. Shoot complex `κ`
over the box `0 ≤ Re κ ≤ 15, |Im κ| ≤ 14`; the analyticity residual through the sonic point is the
objective whose zero locates each mode. Discard the **gauge mode at κ ≈ 0.35699**; the physical relevant
mode is `κ₀ ≈ 2.81055255`. Report `β = 1/Re κ₀`.

## CLI

```
fluidcss_nexus [--json] [--selftest] [--golden] [--stageA] [--scan]
               [--ngrid N] [--xmax X] [--tol T] [--kbox "0,15,-14,14"]
```
- `--stageA` : construct + report the background critical solution and its convergence table only.
- `--scan`   : run the Stage-B complex-κ scan and report all modes found in the box.
- `--json`   : one declared JSON envelope on stdout. `--selftest` : self-checks (< 30 s).
- `--golden` : freeze / verify the blake2b over the declared object.

## Declared output object (hash domain; notes excluded)

```
{ "tool":"fluidcss_nexus", "version":"0.9.0", "units":"G=c=1",
  "params":{ "ngrid":..., "xmax":..., "tol":..., "kbox":[...] },
  "sonic":{ "locus_coeffs":[3,-1,4,-1,3], "V0":..., "N0":..., "A0":..., "omega0":... },
  "background":{ "converged":bool, "n_similarity":..., "A_inf":..., "omega_inf":..., "resid":... },
  "beta":..., "re_kappa0":..., "im_kappa0":..., "n_relevant_modes":...,
  "converge_spread":..., "gauge_mode_kappa":...,
  "gates":{ "G_ANCHOR":bool, "G_CONVERGE":bool, "G_UNIQUE":bool }, "verdict":"pass|fail|blocked" }
```

## Gates

| gate | fires (→ exit 1) unless | meaning |
|---|---|---|
| **G-CONVERGE** | `β` is **stationary under grid/tolerance refinement**: spread of `β` across the refinement ladder `< tol_converge` (declared, e.g. 5×10⁻³) | the metamorphic tooth — a grid-noise "β" is caught |
| **G-UNIQUE** | **exactly one** relevant (`Re κ > 0`, non-gauge) mode in the searched box | KHA's uniqueness of the relevant mode; >1 or 0 ⇒ fail |
| **G-ANCHOR** | `|β − 0.35580192| < tol_anchor` (declared, e.g. 4×10⁻³, ~1%) | the ground-truth check; its PASSING validates the loop |

Battery green (all three fire correctly) ⇒ golden frozen (`goldens/fluidcss/golden.hash`), exit 0.
A fired gate ⇒ exit 1. Error ⇒ exit 2. Never conflated. **`verdict:"blocked"` + exit 1** is the honest
state while Stage B does not converge (a declared negative result, not an error).

## Determinism

Single-threaded fp64, fixed grids/tolerances, no RNG, no wall-clock in the declared path. Shooting uses
deterministic Newton/bisection with declared iteration caps. `(params) → byte-identical declared JSON →
blake2b` (tiny_nexus / substrate_nexus idiom). Build: `cl /std:c++17 /EHsc /O2 /W4`, ≤ 5 min.

## The one gap (§1.6) — status

The explicit sonic-point regularity equation and L'Hôpital local expansion (not printed in KHA) **are
derived and verified** (see above and RESULTS). The remaining blocker is *upstream* of the sonic point:
reconciling the primary-source **fluid** ODEs so the background integrates to a regular center. See
`RESULTS_fluidcss.md` for the precise fault localization (the `ω_,x` coefficients / `2N` source terms)
and the two candidate systems tested.

## Changelog

- v0.9.1 (2026-07-12) — **Stage-A LANDED + banked as a golden** (D-027). The tool now implements the Evans–Coleman background (HKA regular center 4.11–4.13 → sonic 4.5; constraint 4.2 as a first integral eliminating A), ported to deterministic C++ fp64 (RK4 center→sonic shoot on the central density oi so V_sonic=−c_s). **oi*=3/8, sonic point (3/2, 2/√3, 3/4, −1/√3) EXACT, 2m/r=1/3, invariants N=N∞e⁻ˣ & A=1+⅔ω** — all EMERGE (nothing tuned). Golden `fluidcss_stageA` `b4f4e463` frozen (declared Stage-A object: oi*, sonic point, 2m/r, invariant residuals, grid-convergence). Supersedes the walled v0.9.0 4D reduction (which lacked a regular center on the *ingoing* sound cone). **β/Stage-B remains BLOCKED** (perturbation ∂ₛ-coupling; D-016/D-021 — β not reported). Tool `substrate/fluidcss_nexus.cpp`; `substrate/MODULE_fluidcss.md`.
- v0.9.0 (2026-07-11) — DRAFT. Stage-A sonic-point condition derived + verified three ways; contract,
  target, gates, and declared object frozen. Stage-B eigenvalue not yet landing β (fluid-equation
  discrepancy, documented). β **not** reported as measured (D-016 honesty). No golden frozen for a
  β-bearing object until G-ANCHOR/G-CONVERGE/G-UNIQUE fire.
