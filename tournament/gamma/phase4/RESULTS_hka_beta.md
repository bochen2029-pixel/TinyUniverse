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
- **STAGE B — NOT LANDED (precisely diagnosed).** The linear **perturbation operator** L(x;κ),
  assembled from the transcribed HKA coefficients (5.5)–(5.13), **fails a rigorous gauge-mode
  exactness gate** that any correct operator must pass. Consequently **no eigenvalue and no β are
  reported from it** — per the honesty doctrine, a number that would be built on an unvalidated
  operator is not emitted. The failure is localized to the equation level (below).

| gate | status | evidence |
|---|---|---|
| Stage-A regular center (4.11)–(4.13) | **PASS** | A→1, V→0, N=N_inf·e^{−x}, NV→−2/(3γ)=−1/2 |
| Stage-A sonic = sound cone (4.5), constraint (4.2) | **PASS** | Dson=0, C(4.2)=0 identically |
| Stage-A oi* convergence | **PASS** | oi*→3/8 to ~1e-11 across tol, launch depth, 3 integrators |
| Stage-A sonic point vs HKA (4.7–4.9) | **PASS** | (A0,N0,om0,V0)=(3/2,2/√3,3/4,−1/√3) to 12 s.f. |
| Stage-B perturbation-operator validity | **FAIL** | gauge mode (5.20) is NOT an exact solution of L |
| Stage-B eigenvalue κ / β | **not evaluable** | operator not validated ⇒ no κ extracted, none faked |

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
κ, β: NOT MEASURED (Stage-B operator unvalidated). Stage-A: oi*=3/8, sonic=(3/2,2/√3,3/4,−1/√3) EXACT.
```

---

## UPDATE 2026-07-12 (session 2) — the linear-combination diagnosis + a primary-EOM operator (β still walled, sharpened)

**New structural finding — J and K used inequivalent combinations.** The prior reduced operator mixed two *inequivalent linear combinations* of the fluid conservation laws. The reduced flow Jacobian `J3` (`hka_pert_reduced.py`) was built from the code's **recombined** Eq. 4 (`hka_ec.fluid_slopes`: the (ω̄,V) row uses `Cx=(γ−1)(N+V), Dx=γ(1+NV)`), while the s-coupling `As` (5.5) was transcribed from the **primary** KHA95 Eq. 4 (`[4Vω_s+(4V+N+3NV²)ω_x]/ω + 4[(1+V²)V_s+(1+V²+2NV)V_x]/(1−V²)+…`, gr-qc/9503007 `eq:EOM`). For the *background* (∂ₛ=0) any independent pair gives the same solution — so **Stage-A is unaffected** — but the ∂ₛ coupling is **combination-dependent**, so J and K were inconsistent, which is exactly what the gauge-mode gate detects. Confirmed numerically (`hka_solveK.py`): the fluid rows fail `K·D3=−J3[:,0]` by O(0.3)–O(10), the miss **growing toward the sonic point**; the sign is `s=−1` (physically `−Ax⁻¹As`), best near the center.

**A consistent operator (`hka_pert_primary.py`, new).** Rebuilt `L(x;κ)` by linearizing the **primary KHA95 EOM directly** — both ∂ₓ and ∂ₛ from the *same* 4 equations — in the mixed vars (Ā,N̄,ω̄,V). The metric rows match the transcription exactly (Eq1→(5.9) G, Eq2→(5.10) H, verified algebraically). Gauge-mode gate, honest metric `|res|/|Ψ_g|` (the earlier finite-difference `|res|/|dP|` blew up spuriously near the center where the gauge mode is ~constant so `|dP|→1.7e-4` — a metric artifact, NOT the operator): the primary-EOM operator is **markedly better at the sonic point** (`|res|/|dP|` ≈ 0.3–0.8 vs the transcribed operator's ≈5–43) but **still O(1) overall** (`|res|/|Ψ_g|` ≈ 1–3): `L` maps the near-constant center gauge mode to O(1) instead of →0. Closer, not exact — the remaining error is in the near-center metric-row (Ā,N̄) coupling.

**Where the authoritative operator actually lives (the real blocker).** The on-disk primary — **KHA95 PRL gr-qc/9503007** (`KHA95_src/9503007.tex`) — prints the *background* EOM (`eq:EOM`) explicitly but for the perturbation only *states* (§5) that substitution "yields linear, homogeneous first-order ODEs"; it **does not print the perturbation operator**. So neither the transcription (`HKA_beta_equations.md`, itself flagged "5.13 = reconstructed prose") nor a PRL-based re-derivation has the authoritative coefficients. The explicit §V perturbation operator is in the **long paper HKA gr-qc/9607010 (PRD 59, 104008)** — which is **NOT on disk**.

**Cleanest unblock (next attempt):** (1) **fetch gr-qc/9607010 and transcribe §V** (the explicit perturbation coefficients) → drop into `hka_pert_primary`/`hka_pert_core` → the gate should pass; OR (2) finish the primary-EOM derivation's near-center (Ā,N̄) coupling (the residual localized above). Then the already-built, validated match determinant `hka_beta4.Delta` reads κ (discard the gauge modes κ≈0.35699/1 per HKA fn15) → β=1/Re κ.

**Honest state (D-016/D-021).** κ, β: STILL NOT MEASURED, none faked. Stage-A golden `fluidcss_stageA` stands. Two independent operators now fail the gauge gate (transcribed: worst at sonic; primary-EOM: worst at center) — the wall is real and localized, not a tuning problem. `β` is withheld until an operator passes `hka_pert.py`'s gauge-mode gate.

---

## UPDATE 2026-07-12 (session 2, cont.) — **THE OPERATOR IS CRACKED** (gauge gate → 1e-10)

**Fetched the authoritative long paper** (HKA gr-qc/9607010 = PRD 59 104008, `HKA99_src/rflanl.tex` §V `eq:EOM-var`/`eq:EOM-eigenmodes`). Coefficient-by-coefficient, the transcribed (5.5)–(5.10) were **already correct**. **The bug was never the coefficients — it was the ASSEMBLY.** The paper writes the eigenmode system as

> `M_x · Ψ' = (Gmat − κ M_s) · Ψ`,  where **`M_x` is NOT the identity**: `M_x = diag-block(1, 1, [[Ax,Bx],[Cx,Dx]])`.

So the correct operator is **`L = M_x⁻¹ (Gmat − κ M_s)`** — the fluid rows must be multiplied by `[[Ax,Bx],[Cx,Dx]]⁻¹`. The prior `hka_pert_symbols`/`hka_pert_core` assembly and my `hka_pert_primary` both mishandled this `M_x` inversion (treating the fluid ∂ₓ block as if pre-inverted / a linear-combination mismatch).

**Verified (`hka_pert_hka99.py`, the authoritative operator):** the gauge-mode exactness gate now passes to **machine precision — `|res|/|Ψ_g| ≈ 1e-10` for every κ̄ ∈ {0.357, 1.0, 2.81}**, at every x. The κ¹ gauge condition `E2 = As·ω̄_x + Bs·V_x`, `F2 = Cs·ω̄_x + Ds·V_x` holds to ~1e-15 on the background (it is the background Eq 3/4). **This is the wall broken** — the exact thing that blocked the overnight run and session-2's first pass.

**Remaining: the eigenvalue extraction (β = 1/Re κ₀, κ₀≈2.81055255).** The paper (Table I) confirms γ=4/3 → κ=2.81055255, β=0.35580192, unique relevant mode. Getting κ from the *correct* operator is now a standard-but-delicate singular-BVP shoot (regular-singular sonic point + stiff center; mode collapse under long integration). The prior `hka_beta4` match machinery (calibrated for the wrong operator) is too noisy; a robust BVP eigensolver (spectral collocation / orthonormalized shoot / `solve_bvp`) on the verified `hka_pert_hka99.Lnum` is the clean next step. Tools: `hka_pert_hka99.py` (operator, GATE-VERIFIED), `hka_beta_solve.py` (shoot scaffold). β still withheld until κ lands + G-UNIQUE/G-CONVERGE fire — the machine is armed on a correct operator for the first time.
