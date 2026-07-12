# RESULTS — HKA radiation-fluid critical exponent (fluid-β track), fresh HKA Eq. 4.1 implementation

**Date:** 2026-07-12 · **Branch:** `substrate-gamma-tournament` · **Equation source:**
`HKA_beta_equations.md` (Hara–Koike–Adachi gr-qc/9607010 = PRD 59 104008; Maison gr-qc/9504008;
KHA gr-qc/9503007; Evans–Coleman gr-qc/9402041) · **Target:** κ = 2.81055255 ⇒ β = 1/Re κ = 0.35580192.

## BOTTOM LINE (honest)

- **STAGE A — LANDED.** The correct HKA self-similar **background** (the Evans–Coleman critical CSS
  solution) is obtained fresh from a **genuine regular center** on the **ingoing** sound cone, shooting
  to a real sonic point. All Stage-A invariants match analytic values to **machine precision**, with
  full convergence receipts. This **resolves the prior wall** (which had put the sonic point on the
  wrong branch V=+1/√3, lacked a regular center, and gave β≈0.99).
- **STAGE B — NOT LANDED (precisely diagnosed).** SUPERSEDED by `STAGE_B_UPDATE.md` (v1: operator
  CORRECTED via direct linearization → passes the gauge gate) and `STAGE_B_UPDATE_2.md`-content in
  `STAGE_B_UPDATE.md` (this pass: robust eigenvalue solvers built). **Current honest state:** the
  corrected, gauge-exact operator (`hka_pert_derive.py`) has correct center {−2,−1,0,0} and sonic
  {0,0,0,1−2κ} indicial structure, but a **compound-matrix Evans function** (amplification-immune,
  VALIDATED by cleanly finding the gauge mode at **κ=1**) shows it carries **NO physical eigenvalue near
  κ=2.81** (|E(2.81)| is δ-robust on the baseline; a non-convergent FD-BVP "ghost" near 2.79 fails every
  refinement test). Diagnosis: the operator is under-determined by **one** physical ∂_s coupling that the
  gauge mode + indicial structure cannot fix and the available equation transcriptions do not pin down.
  **No β emitted; none tuned.** The Stage-A section below remains valid and current.

| gate | status | evidence |
|---|---|---|
| Stage-A regular center (4.11)–(4.13) | **PASS** | A→1, V→0, N=N_inf·e^{−x}, NV→−2/(3γ)=−1/2 |
| Stage-A sonic = sound cone (4.5), constraint (4.2) | **PASS** | Dson=0, C(4.2)=0 identically |
| Stage-A oi* convergence | **PASS** | oi*→3/8 to ~1e-11 across tol, launch depth, 3 integrators |
| Stage-A sonic point vs HKA (4.7–4.9) | **PASS** | (A0,N0,om0,V0)=(3/2,2/√3,3/4,−1/√3) to 12 s.f. |
| Stage-B perturbation-operator gauge gate | **PASS** (v1 fix) | corrected operator: gauge residual ~1e-9, invariant in κ̄ (`hka_pert.py`, `hka_pert_derive.py`) |
| Stage-B operator indicial structure | **PASS** | center {−2,−1,0,0}, sonic {0,0,0,1−2κ} exact |
| Stage-B eigenvalue solver (Evans, validated) | **PASS on gauge / no physical** | κ=1 gauge found (\|E\| dip 4.4×); κ=2.81 NOT an eigenvalue (\|E\| δ-robust on baseline) |
| Stage-B physical κ / β | **NOT REPRODUCED** | operator missing 1 ∂_s coupling ⇒ no physical mode; no κ/β emitted, none tuned |

---

## STAGE A — Evans–Coleman background (HKA Eq. 4.1, γ=4/3)  ✅

### Method (fresh, correct)
The prior attempt walled because it put the sonic point at **V=+1/√3** (a from-scratch derivation that
gave β≈0.99) and lacked a regular center. Here the **ingoing** structure is used, exactly per HKA:

1. **Sonic data (4.7)–(4.9)** parameterised by V₀ — verified symbolically (`hka_verify_sonic.py`) to
   sit on the sound cone (4.5), satisfy the L'Hôpital condition (4.6) **and** the constraint (4.2)
   (all three residuals ≡ 0).
2. **Constraint reduction (the key stabiliser).** The constraint (4.2) obeys **dC/dx = −A·C**
   (`hka_constraint.py`) — a first integral. C=0 is only *marginally* stable, so integrating the full
   4-D flow toward the center amplifies tiny C-violations exponentially (this is the chaos the prior
   effort hit). **Cure:** eliminate A via (4.2), A = 1 + 2ω[1+(γ−1)V²]/(1−V²) + 2γNVω/(1−V²), and
   integrate the **reduced 3-D** system (N,ω,V) — exactly on C=0, no runaway.
3. **Shoot from the REGULAR CENTER outward** (the natural direction: the center has exactly ONE
   unstable direction). Launch from the center asymptotics (4.12)–(4.13) [verified against the ODE to
   ~1e-13], with gauge N_inf and central-density parameter **oi**. The critical oi is the one whose
   trajectory reaches the sonic point analytically (V at the Dson=0 crossing = −1/√3).

### Result (converged, machine precision — `hka_stageA_convergence.py`, `hka_stageA_verify.py`)
```
critical central density   oi* = 0.375000000 = 3/8   (gauge N_inf=1)
   convergence: oi*→3/8 as rtol→0 (|oi*−3/8| = 6e-12 at rtol=1e-13);
   robust over launch depth x0∈[−16,−8] and integrators RK45/DOP853/Radau (agree to ~1e-11).

sonic point (A0, N0, om0, V0) = (3/2, 2/√3, 3/4, −1/√3)   EXACTLY  (matches HKA 4.7–4.9 to 12 s.f.)
   Misner–Sharp  2m/r = 1 − 1/A0 = 1/3   exactly
   flow speed |V0| = 1/√3 = c_s  (sonic = sound cone);  Dson=0, constraint C(4.2)=0 identically.

EXACT invariants along the EC solution (new findings, verified to ~1e-13):
   N = N_inf · e^{−x}       (i.e. (lnN)' ≡ −1)
   A = 1 + (2/3) ω          (equivalently A = 1 + (2−γ)ω)
```
The center is genuinely regular (A→1, V→0, N→+∞, NV→−2/(3γ)=−1/2, ω→0). This is the first
Bogoyavlenskii/Evans–Coleman solution with V<0 (ingoing) throughout the center→sonic segment (0 zeros
of V there; the single V-zero of the EC solution lies in the outer region x>x_s).

### Provenance / cross-checks
- `hka_desing.py`: independent desingularised-flow derivation whose resolved RHS matches the direct
  `hka_ec.rhs3` **exactly** — two independent derivations of (4.1) agree.
- The sonic point (3/2, 2/√3, 3/4, 1/√3) is the same one the prior (mirror-branch) machinery used
  (`stageB_qr.py` `_precompute_ML`), corroborating the location.

---

## STAGE B — perturbation eigenvalue (β)  ❌ (diagnosed)

### The gate that failed
A **pure-gauge** perturbation (HKA 5.20), Ψ_g(x;κ̄) = ((lnA)'_x, (lnN)'_x + κ̄, (lnω)'_x, V'_x), is an
**exact** solution of Ψ' = L(x;κ̄)Ψ for **every** κ̄ (it is pure coordinate freedom). This is a
rigorous falsification test: a correct L gives |dΨ_g/dx − L Ψ_g| ≈ 0. The operator assembled from the
transcribed (5.5)–(5.13) gives **relative residuals O(1)–O(10)** (`hka_pert.py`) — so **L is not the
correct perturbation operator**, and no eigenvalue is trusted from it.

### What IS verified (the operator is close, and the wall is localized)
- **Abar row + gauge-mode form: CORRECT.** The auxiliary identity (5.14) holds **exactly** on Ψ_g
  (residual ~1e-11, `hka_check_514.py`), and independently the Abar-row `−κĀ_p + (5.14 RHS)` reproduces
  dĀ_p/dx on Ψ_g to ~1e-11 (`hka_check_rows.py`).
- **Reduced flow Jacobian J3: CORRECT (κ=0).** On the constraint surface the background flow Jacobian
  J3 reproduces the κ=0 gauge mode (the flow tangent) to ~1e-9 (`hka_pert_reduced.py`). Its indicial
  exponents are {−2,0,0} at the center (modes 1,1,e^{−2x}) and {−19.2, −0.71±0.76i} at the sonic point
  (one singular ~1/t mode + the analytic pair) — the correct regular-singular structure.
- **Sonic Frobenius structure: CORRECT.** L has a genuine 1/t pole at the sonic point with indicial
  exponents {0,0,0, μ} where μ=1+2κ (our t=x−x_s sign; ⇔ the md's 1−2κ under t→−t) — i.e. **3 analytic
  modes + 1 non-analytic mode**, exactly as HKA (`hka_pert_sonic.py`, `hka_frobenius.py`).

### Where it breaks (equation-level)
The **κ-coupling of the Nbar and (ω̄,V) rows** is wrong as transcribed. First-principles analysis
(gauge invariance): writing L = J3 + κK on the constraint surface, exactness of Ψ_g for all κ̄ **forces**
`K·e_lnN = 0` and `K·D3 = −J3[:,0]`, where D3 = (background log-derivatives). These are 6 constraints on
the 9 entries of K — **under-determined by 3**, and the transcribed HKA s-coupling matrix Ms (5.5),
assembled as Minv·Ms for the fluid block, does **not** satisfy even the determined part (rows 1,2 of
`K·D3 = −J3[:,0]` fail by O(0.3)–O(10), `hka_solveK.py`). So the (5.5)–(5.10) coefficients as written
do not supply the correct s-coupling; the missing 3 degrees of freedom are the genuine ∂_s content of
the HKA PDEs, not recoverable from the gauge mode alone.

Systematic falsification of the assembled operators: `hka_pert_symbols.py` / `hka_pert_symbols2.py`
(5.13 direct + 5.14 Abar), `hka_pert_jac.py` / `hka_pert_jac2.py` (flow-Jacobian + guessed K),
`hka_pert_red3.py` (reduced, all gauge-shift directions × fluid-sign) — **none** pass the gate for all κ̄.

Consistent with the gate failure, the two-sided eigenvalue matches built on this operator
(`hka_beta.py` det[Qc|Qs]; `hka_beta4.py` det[v_center|s1|s2|s3] with QR — the prior-verified
non-degenerate scheme) show **no eigenvalue near κ=2.81**, and the positive control (`hka_validate_gauge.py`:
the known gauge mode at κ=0.35699 must give rejection→0) **also fails** — precisely because the
operator is wrong. The QR/rejection machinery itself is scale-invariant and healthy (residuals O(1));
it is the operator that is not.

### What would close Stage B
Obtain the **correct perturbation κ-coupling**: either re-transcribe HKA (5.5)–(5.10)+(5.13) from the
primary TeX (the md itself flags 5.13 as "assemble from 5.5–5.10", i.e. reconstructed prose), OR
re-derive the linear perturbation of the p=ρ/3 self-similar system from the time-dependent
∇_aT^{ab}=0 + Einstein equations keeping ∂_s (the ∂_s coefficients are the missing 3 dof). Then the
gauge-mode gate (`hka_pert.py`) will pass (residual→0 at κ=0.357 and κ=1), after which the **already-built,
validated** two-sided QR shoot (`hka_beta4.py`) — with the gauge modes κ≈0.35699/1 DISCARDED per HKA
fn15 — extracts the physical κ. All downstream machinery (Frobenius modes, center-regular subspace,
match determinant) is in place and tested; only the operator is blocking.

---

## Deliverables (all in `tournament/gamma/phase4/`)

**Stage A (landed):**
- `hka_background.py` — clean public API (re-exports the verified solver; prints oi*=3/8 + exact sonic).
- `hka_ec.py` — the EC shoot (center→sonic, solve_ivp + Dson event, brentq on V_cross=−1/√3).
- `hka_stageA_verify.py`, `hka_stageA_convergence.py` — the verification + convergence receipts.
- `hka_verify_sonic.py`, `hka_constraint.py`, `hka_desing.py` — sonic/constraint/desingularisation checks.

**Stage B (diagnosed, not landed):**
- `hka_pert.py` — the gauge-mode VALIDATION GATE + status (the honest gate).
- `hka_pert_symbols.py`, `hka_pert_core.py` — the (5.5)–(5.13) operator assembly + numeric L.
- `hka_frobenius.py`, `hka_pert_sonic.py` — sonic Frobenius analytic modes + residue/indicial structure.
- `hka_beta.py`, `hka_beta4.py` — two-sided eigenvalue matches (QR-stabilised) on the operator.
- `hka_check_514.py`, `hka_check_rows.py`, `hka_pert_reduced.py`, `hka_solveK.py`,
  `hka_pert_symbols2.py`, `hka_pert_jac*.py`, `hka_pert_red3.py`, `hka_validate_gauge.py` —
  the row-by-row verification and the systematic operator falsification.

## Reproduce
```
python hka_background.py            # Stage A: oi*=3/8, sonic=(3/2,2/√3,3/4,−1/√3), 2m/r=1/3
python hka_stageA_convergence.py    # Stage-A convergence receipt (tol / depth / integrator)
python hka_verify_sonic.py          # sonic data on cone, L'Hôpital (4.6), constraint (4.2) all ≡0
python hka_pert.py                  # Stage-B gate: operator FAILS gauge-mode exactness (blocked)
```

## Doctrine note
This is the intended honest outcome. **Stage A is a real, converged landing** (the correct HKA
background with a genuine regular center and an exact sonic point — a substantive advance past the
prior wall). **Stage B is blocked at a precisely-localized equation** (the perturbation κ-coupling),
**with no β faked and none tuned toward 0.3558**. The eigenvalue machinery is built and validated
against a scale-invariant gate; it is armed and waiting on a correct perturbation operator.

### Exact resume/port state
- **Port-ready NOW:** the Stage-A background (`hka_ec.ec_background` / `hka_background.ec_background`)
  → `substrate/fluidcss_nexus.cpp` background sector + a Stage-A golden over
  (oi*, sonic point, 2m/r, invariants N=e^{−x}, A=1+2ω/3). β is NOT part of that golden.
- **Blocked:** freezing any β golden — withheld until the perturbation operator passes `hka_pert.py`'s
  gauge-mode gate. Next action: fix the perturbation κ-coupling (re-transcribe 5.5–5.10 from primary
  TeX, or re-derive ∂_s), re-run `hka_beta4.py`, discard κ≈0.357/1, report the physical κ.
```
κ, β: NOT MEASURED convergently. Corrected operator is gauge-exact + correct indicial structure, but the
      Evans function (validated on the κ=1 gauge mode) finds NO physical eigenvalue at κ=2.81 — the operator
      is missing one physically-undetermined ∂_s coupling. See STAGE_B_UPDATE.md for the full five-method
      diagnosis. Stage-A: oi*=3/8, sonic=(3/2,2/√3,3/4,−1/√3) EXACT. No β faked or tuned toward 0.35580192.
```
