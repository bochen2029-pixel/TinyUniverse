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

**DONE (session 2, cont.):** unblock (1) executed — gr-qc/9607010 fetched, §V transcribed, and the operator CRACKED (M_x≠I assembly; gauge gate → 1e-10, `hka_pert_hka99.py`). See the section below for the eigenvalue-extraction state.

---

## UPDATE 2026-07-13 — eigenvalue extraction: correct machinery built, blocked at the sonic-point numerics (β not yet landed, none faked)

With the **verified** operator (`hka_pert_hka99.Lnum`, gauge gate 1e-10), the remaining task is to extract κ₀≈2.81055255 (→ β=0.35580192). This is the *known-hard* part of critical-collapse — the sonic point is a **regular singular point** and the extraction is delicate (HKA/Maison/Gundlach each engineered specific techniques). Eight formulations were built and run; all correctly recover the **gauge mode** but the **physical mode is masked**. Precise findings, all measured:

- **The algebraic identity (eq:alg-PP / eq:EOM-eigenmodes.2) is NOT preserved by the 4-row ODE** (measured: an identity-exact vector drifts to 24% off-identity by x≈−2). The paper uses it "as a check," but the physical solutions are the 3-D subspace satisfying it. ⇒ any 4-row shoot/collocation admits identity-violating modes and **misses the physical one**. Confirmed: 4-D `solve_bvp` and 4-D Chebyshev spectral find only {−1.67, −0.40, 0.994} — the **gauge mode at κ≈0.994** (in this background's gauge, where N̄_p(x_s)=0 ⇒ κ=−N̄′(x_s)≈1), never 2.81.
- **Fix — the identity-reduced 3-D system** (`Ā_p` eliminated via eq:alg-PP; `hka_beta_match.py::L3`): identity built in, no drift. Its center indicial exponents are correct ({−2,0,0}); **but near the sonic point it has a violent positive eigenvalue (+7.4 at x=−0.5, larger toward x_s)**. Backward integration (sonic→center) blows up ~1e57 (Frobenius truncation excites the +7 mode); forward integration (center→sonic) is stable in the bulk but the match determinant **explodes to −1e5 in the physical κ band**, swamping the eigenvalue zero. The reduced-3D forward 2×2 match determinant cleanly shows the **gauge sign-change at κ≈1.0**, but no resolvable zero near 2.81.
- **Correctly built + validated pieces:** the physical sonic leading coefficient `a₀=(Ā=1, N̄=0, ω̄=0.014, V=−0.76)` (3×3 solve: Ā₀=1 norm + gauge N̄₀=0 + identity=0, cond≈10, identity residual 2e-16, κ-dependent) — `hka_beta_solve.py::v_sonic`/`hka_beta_match.py`; the Frobenius analytic sonic modes (`ker(R)`, 3-D, indicial {0,0,0,−4.48}); the sonic solvability BC (fluid-block left-null `(Cx,−Ax)`).

**Precise remaining step (well-scoped, literature-standard):** bridge the sonic region by **high-order analytic continuation** — evaluate the physical sonic-analytic solution via its Frobenius series out to where the +7 instability has not yet dominated, match to the forward-integrated center-regular pair there (a matched two-sided shoot away from the singular point), OR reformulate in **double-null / CSS-similarity coordinates** (the tournament's F1 discussion) where the self-similar solution is smooth-on-grid and the sonic instability is tamed. Then the eigenvalue drops out; discard the gauge root at κ≈1 (this background's gauge).

**Honest state (D-016/D-021):** β STILL NOT MEASURED, none faked. **The operator is the hard scientific crux and it is DONE + verified.** The eigenvalue is a bounded numerical-methods task on a correct operator — the machinery (reduced-3D + identity + Frobenius) is built and characterized; only the sonic-point analytic-continuation step remains. Scaffolds: `hka_beta_match.py` (reduced-3D, best), `hka_beta_spec.py` (spectral), `hka_beta_bvp.py` (collocation), `hka_beta_solve.py`, `hka_beta_validate.py`, `hka_beta_eigfn.py`.

### 2026-07-13 (cont.) — the matched-shoot attempt + the definitive diagnosis (polar-areal Frobenius is a dead end for the new operator)

Attempted the precise matched shoot (`det[c1|c2|v_sonic]` at a point outside the singular region but inside the Frobenius radius). It **cannot work as built**, because the sonic **Frobenius machinery is inaccurate for the correct operator**:
- **The physical Frobenius series has O(1)–O(100) ODE-residual at every |t| (0.05→0.20), worsening with order** — it is not a usable solution of the ODE anywhere. Root: `hka_pert_sonic.bg_series_near_sonic` fits the background near the sonic point with a **low-order polyfit of sampled points**, and `L_laurent`'s DFT then inherits that error.
- **The residue R (indicial exponent) is not even pinned down:** DFT `L_laurent` → μ≈−4.48; an analytic `R = (1/d′)·adj(fluid-block)·(Gmat−κM_s)` → μ≈+10.5; the old-code "expectation" was 1−2κ=−4.62. Three different values ⇒ the sonic expansion must be **rebuilt from scratch** (exact background sonic Taylor series from the desingularized flow, exact Laurent) before any Frobenius/matched-shoot can be trusted — and even then the small radius (<0.3) + the +7/singular amplification make the polar-areal shoot delicate.

**Definitive conclusion:** the polar-areal ODE-shooting/collocation route to β is a **dead end without a full rebuild of the sonic analytic-expansion machinery**, and remains delicate even then. **The robust path is a different method, both standard in the literature:**
1. **The Lyapunov / time-evolution method** (HKA §V's *own* "second method", called "very regular"): evolve the linearized `(s,x)` PDE — fluid `(ω̄,V)` evolve in `s` via `[[As,Bs],[Cs,Ds]]⁻¹`, metric `(Ā,N̄)` from the spatial constraints (`∂_xĀ=G·Ψ`, `∂_xN̄=H·Ψ`) — from generic data; the dominant growth rate is Re κ ⇒ β. Avoids the sonic ODE eigenvalue entirely. **Recommended.**
2. **Double-null / CSS-similarity reformulation** (tournament F1): re-derive the operator in coordinates where the self-similar solution is smooth-on-grid.

Either is a focused, dedicated implementation on the **already-verified operator** (`hka_pert_hka99.py`). β is withheld; the operator crank stands.

**Lyapunov spectral method attempted (`hka_beta_lyap.py`) — also hits a numerical wall.** Built the s-evolution generator Q (regular at the sonic point, as predicted: metric rows → spatial constraints `∂_xĀ=G·Ψ`, `∂_xN̄=H·Ψ`; fluid rows evolve via `[[As,Bs],[Cs,Ds]]⁻¹`), discretized on a Chebyshev grid, `eig(Q)`. Result: the spectrum is **dominated by huge O(N²) spurious eigenvalues** (the advective `∂_x` term + the near-singular constraint-solve inversions), with only one clean O(1) real mode (κ≈−13.4); the physical κ=2.81 is not cleanly separable without careful spectral filtering (companion-form deflation), sonic-point-aware upwinding, and proper center/+∞ radiation BCs.

**FINAL (2026-07-13): nine distinct methods** — Frobenius-match (`hka_beta4`), rejection, growth shoot (4-D & reduced-3D), `solve_bvp` collocation, Chebyshev spectral (4-D & identity-augmented), reduced-3D forward match determinant, and the Lyapunov generator — **all extract the gauge mode but wall on the physical κ=2.81.** The blocker is uniformly the **regular-singular sonic point** (± its downstream effects: identity non-preservation, +7 instability, tiny/uncertain Frobenius radius, spurious spectral modes). This is the known-hard core of critical-collapse eigenvalue computation, which Choptuik/Maison/HKA/Gundlach each solved with dedicated, carefully-engineered schemes (AMR-in-null, compound-matrix Evans, or the fully-worked Lyapunov PDE with proper BCs). **The operator — the scientific crux — is DONE and verified (1e-10). Landing β is now a focused NR-numerics research task on that correct operator, not a brute-force sweep.** β NOT measured, none faked (D-016/D-021). All scaffolds committed for the dedicated next effort.

**Honest state (D-016/D-021).** κ, β: STILL NOT MEASURED, none faked. Stage-A golden `fluidcss_stageA` stands. Two independent operators now fail the gauge gate (transcribed: worst at sonic; primary-EOM: worst at center) — the wall is real and localized, not a tuning problem. `β` is withheld until an operator passes `hka_pert.py`'s gauge-mode gate.

---

## UPDATE 2026-07-12 (session 2, cont.) — **THE OPERATOR IS CRACKED** (gauge gate → 1e-10)

**Fetched the authoritative long paper** (HKA gr-qc/9607010 = PRD 59 104008, `HKA99_src/rflanl.tex` §V `eq:EOM-var`/`eq:EOM-eigenmodes`). Coefficient-by-coefficient, the transcribed (5.5)–(5.10) were **already correct**. **The bug was never the coefficients — it was the ASSEMBLY.** The paper writes the eigenmode system as

> `M_x · Ψ' = (Gmat − κ M_s) · Ψ`,  where **`M_x` is NOT the identity**: `M_x = diag-block(1, 1, [[Ax,Bx],[Cx,Dx]])`.

So the correct operator is **`L = M_x⁻¹ (Gmat − κ M_s)`** — the fluid rows must be multiplied by `[[Ax,Bx],[Cx,Dx]]⁻¹`. The prior `hka_pert_symbols`/`hka_pert_core` assembly and my `hka_pert_primary` both mishandled this `M_x` inversion (treating the fluid ∂ₓ block as if pre-inverted / a linear-combination mismatch).

**Verified (`hka_pert_hka99.py`, the authoritative operator):** the gauge-mode exactness gate now passes to **machine precision — `|res|/|Ψ_g| ≈ 1e-10` for every κ̄ ∈ {0.357, 1.0, 2.81}**, at every x. The κ¹ gauge condition `E2 = As·ω̄_x + Bs·V_x`, `F2 = Cs·ω̄_x + Ds·V_x` holds to ~1e-15 on the background (it is the background Eq 3/4). **This is the wall broken** — the exact thing that blocked the overnight run and session-2's first pass.

**Remaining: the eigenvalue extraction (β = 1/Re κ₀, κ₀≈2.81055255).** The paper (Table I) confirms γ=4/3 → κ=2.81055255, β=0.35580192, unique relevant mode. Getting κ from the *correct* operator is now a standard-but-delicate singular-BVP shoot (regular-singular sonic point + stiff center; mode collapse under long integration). The prior `hka_beta4` match machinery (calibrated for the wrong operator) is too noisy; a robust BVP eigensolver (spectral collocation / orthonormalized shoot / `solve_bvp`) on the verified `hka_pert_hka99.Lnum` is the clean next step. Tools: `hka_pert_hka99.py` (operator, GATE-VERIFIED), `hka_beta_solve.py` (shoot scaffold). β still withheld until κ lands + G-UNIQUE/G-CONVERGE fire — the machine is armed on a correct operator for the first time.

---

## UPDATE 2026-07-13 (session 3) — the sonic Frobenius is REBUILT + EXACT; **but κ=2.81 is ABSENT from the operator** (the wall moves from numerics to the operator)

Two results this session, one constructive and one a hard, triangulated negative that **redirects the whole β effort**.

### A · The sonic analytic machinery was rebuilt from scratch and is now EXACT (the §5a step-1 deliverable — DONE)

The prior Frobenius/Laurent path was poisoned by `hka_pert_sonic.bg_series_near_sonic`: a one-sided `np.polyfit` of the **singular resolved ODE** sampled right against the sonic point → the residue μ was never pinned (three disagreeing values −4.48/+10.5/−4.62). Rebuilt analytically:

- **`nr_sonic.py`** — EXACT background sonic Taylor series `(A,N,ω,V)(t)`, `t=x−x_s`, by an order-by-order power-series recursion. **Key structural finding:** the order-1 fluid solvability is **quadratic** in the null-coefficient α₁ (two analytic branches through the sonic point — `M₁·w₁` is a product of two α₁-linear factors); the Evans–Coleman branch is selected via the desingularized-flow eigenvector (`hka_desing`). Orders ≥2 are linear. **Validated reference-free**: ODE residual ~1e-15 at low order; matches the cleanly-integrated EC background to **3.7e-11 at t=−0.02**. (radius ≈0.12 → match inside |t|≤0.05.)
- **`nr_laurent.py`** — EXACT Laurent `L(t)=M_x⁻¹(Gmat−κM_s)=R/t+L₀+L₁t+…` by power-series arithmetic (no DFT). **Residue reconciled: the indicial exponent μ = 1−2κ EXACTLY** (−4.621105 for κ=2.81; −1.0 for κ=1; +0.286 for κ=0.357) — the −4.62 "expectation" was right; the DFT/polyfit were artifacts. All **3 analytic modes** (ker R, clean rank-1) now recovered; they solve the **direct operator ODE** to **9.8e-14 at t_m=−0.02**, 2.9e-10 at t_m=−0.03. **The Frobenius GATE the continuation prompt asked for PASSES.**

This fixes the exact inaccuracy the prior session named as the root blocker. It is banked (`nr_sonic`, `nr_laurent`, `nr_frob`) and reusable.

### B · With the exact machinery, HKA's OWN eigenvalue method finds ONLY the gauge mode — **κ=2.81 is not there**

`rflanl.tex` §V (read verbatim, 5.216–5.463) gives HKA's exact BVP + numerical method: construct the **unique** analytic-at-sonic solution (a₀∈ker R, + the identity eq:alg-PP, + gauge **N̄_p(0)=0**, normalized by Ā_p(0)=1 → 1 param κ), integrate **sonic→center**, and require the single expanding center mode `(0,0,0,1)e^{−2x}` (pure V_p) to vanish (eq:PPasmp1). +∞ is automatically bounded for Re κ>0 (not an extra condition). Implemented faithfully in **`nr_shoot.py`** with the exact Frobenius sonic data and the **verified identity** (`idc·Ψ_gauge ≈ 1e-16` at all x, both κ — the identity is correct).

**Result (six independent methods agree — β NOT measured, none faked):**

| method | finds gauge κ=1 | finds physical κ=2.81 |
|---|---|---|
| `nr_shoot` sonic→center, 4-D (accurate Frobenius) | ✅ dip (|y|≈9) | ❌ none |
| `nr_shoot` reduced-3D (identity built in, no drift) | ✅ **log\|y\|=0.0** (razor-clean) | ❌ none |
| `hka_beta_solve` sonic→center (old polyfit Frobenius) | — | ❌ monotonic, no min |
| `hka_beta_match` reduced-3D 2×2 determinant | ✅ sign-change | ❌ det explodes monotonically |
| `nr_evans2` two-sided Evans (gauge-removed) | (n/a) | ❌ flat \|E\|≈0.7 |
| complex-κ scan (Re∈[2.6,3.0], Im∈[0,2.5]) | — | ❌ flat log\|y\|≈11.5 |

Comprehensive real scan **κ∈[−1.5, 12]** (HKA's own search range Re κ≥−1.5): the **only** eigenvalue is the gauge mode at **κ=1** (= −N̄_x(x_s) = −(−1); the gauge-fixed gauge mode). The growing-mode coefficient is smooth and sign-consistent through 2.81 with **no zero crossing**; no scaled/complex sibling of 2.81 exists in range.

**Independent structural checks on the operator all PASS** yet the physical mode is absent: gauge-mode gate (1e-10), **center indicial = {−1,0,0,−2}** (matches eq:PPasmp1 exactly), **sonic indicial = 1−2κ** (exact). So the gauge gate is **necessary but NOT sufficient**.

### C · Conclusion + redirect (honest, RAYFORMER-style)

**The bottleneck is not the sonic-point numerics** (now exact, gated) — it is **upstream**: the verified-by-gauge-gate operator `hka_pert_hka99` does **not possess the known physical relevant eigenmode** via HKA's own method. The most likely cause is a **residual error in an operator coefficient** that annihilates the gauge mode and preserves the leading indicial exponents (hence invisible to every check run so far) but shifts/removes the physical spectrum — a class of error the M_x≠I fix (D-029) plausibly did not fully clear. (A subtle background/identity inconsistency is less likely: Stage-A is machine-precision and `idc·Ψ_g≈1e-16`.)

**Next effort (well-scoped):** verify the operator *beyond the gauge mode*.
1. Find a **second independent operator test** (e.g., derive eq:alg-PP as a first integral **from** `L` and check self-consistency of every row; or verify each `M_s`/`E`/`F` coefficient against a fresh linearization of the primary KHA95 EOM keeping ∂_s).
2. **Coefficient sensitivity**: the gauge gate fixes only the combinations `E2=As·ω̄_x+Bs·V_x`, `F2=Cs·ω̄_x+Ds·V_x`; the *individual* `As,Bs,Cs,Ds` (and `E3,F3` etc., which vanish on the gauge mode's structure) are under-constrained by the gauge gate. Re-transcribe/re-derive these specifically and re-run `nr_shoot` — the κ=2.81 dip should appear once the operator is right.
3. Only after κ=2.81 lands: port to `fluidcss_nexus.cpp` Stage-B, freeze the β golden.

**The scaffolds are all in place and validated** (`nr_shoot.py` is the clean eigenvalue reader; it correctly resolves the gauge mode to log|y|=0). The moment the operator is corrected, β falls out with no further numerical work.

```
κ, β: STILL NOT MEASURED, none faked. Sonic Frobenius: EXACT (μ=1−2κ, gate 9.8e-14).
Finding: operator lacks the physical κ=2.81 eigenmode (6 methods); wall moved numerics → operator.
```

---

## 2026-07-13 · SESSION 3 (cont.) ADDENDUM — D-031 CORRECTS §C above: the operator is PROVEN correct; the wall is the EXTRACTION method

**§C's conclusion ("most likely a residual error in an operator coefficient") is FALSIFIED by measurement (D-031).** The "next effort" §C scoped was run and came back the other way:

- `nr_rederive.py` symbolically linearizes the primary EOM (`rflanl.tex` 4.1/4.2, keeping ∂ₛ→κ) and re-derives ALL 16 fluid operator coefficients FROM SCRATCH: **every one matches `hka_pert_hka99` (=§V) EXACTLY (diff=0)** — including the gauge-under-constrained `As,Bs,Cs,Ds,E3,F3` that §C(2) suspected. Metric rows (G,H) verified by hand from 4.3/4.4.
- The background is the true EC solution (`N∝e^{−x}`, `A=1+(2−γ)ω` hold EXACTLY; two independent methods), and the identity is verified (`idc·Ψ_g≈1e-16`).
- The 0.35699 gauge-mode footnote is a **red herring**: with the EC sonic values the EOM (4.4) forces N̄'(0)=−1 ⇒ gauge κ=1 (which the reader finds); 0.35699 would require ω₀<0 (unphysical).

**So operator + background + identity are all PROVEN; κ=2.81 is absent from the sonic→center SHOOT (6 methods) ⇒ the wall is the eigenvalue-EXTRACTION method** (the known-hard numerics of critical collapse), not the operator. HKA's own *second* method sidesteps the sonic point entirely: the **Lyapunov time-evolution** (§V.G, "very regular" — M_s is non-singular, so the s-evolution has NO sonic singularity; power-iterate, dominant growth rate = the relevant κ). **NEXT = implement it** (sonic-aware upwinding + center-regular/gauge BCs). Full plan + ranked steps: `docs/CONTINUATION_2026-07-13_session3.md` §5.

**QC re-verification (2026-07-19, all cold):** gauge gate `|res|/|Ψ|`≈1e-10 ✓ · `nr_rederive` all 16 diff=0 ✓ · `nr_laurent` μ=1−2κ exact + Frobenius gate 9.8e-14 @ t=−0.02 ✓ · `nr_sonic` matches integrated EC to 1.5e-10 @ t=−0.02 ✓ · `nr_shoot` gauge positive control κ=0.99998705 (≈1 to 1.3e-5; the recorded 0.99999924 differs in the last digits — env/tolerance-sensitive, same mode) ✓.

```
κ, β: STILL NOT MEASURED, none faked. Operator: PROVEN (D-031). Wall: the extraction method.
Next: HKA Lyapunov §V.G power iteration → κ → β → C++ port → fluidcss_nexus v1.0.0.
```

---

## 2026-07-19 · SESSION 4 (autonomous) — THE WALL FALLS: the background was FRIEDMANN; on the true EC background β LANDS

**D-032.** The three-session wall (D-029/030/031) is closed. Chain of discovery, each step measured:

1. **Lyapunov §V.G implemented** (`nr_lyap.py`): ξ=e^x uniform grid (bounds the center advection speed — HKA's "coordinate transformations"), metric slaved by vectorized integrating-factor/cumtrapz, MOL-RK4 + KO, free sonic BC, center penalty BCs. **Built-in analytic control: the gauge mode grows at 1.000047 with frozen shape.** Formulation B adds the **fifth equation** (`rflanl.tex` l.5090 — the linearized momentum constraint, an evolution equation for Ā; its eigen-form + row 0 derives eq:alg-PP: cN=c₁, com=G₃+c₁, cV=G₄+c₂ verified). One-step-map spectrum: **on the old background the discrete dynamics genuinely lack κ=2.81** — {gauge 1.0006, then decaying}.
2. **The banked background is the collapsing flat FRIEDMANN solution, not Evans–Coleman.** Fingerprints: A=1+(2−γ)ω is the FRW radiation identity 1/(1−H²r²); oi\*=3/8=(3/2)H²t²; 2m/r=1/3=V₀²; measured **V≡−√(1−1/A) to 1.9e-10** along it. Cause: `hka_ec` used "V=−c_s" as the sonic criterion, but §IV parametrizes sonic points by FREE V₀ (4.7–4.9) and **V₀=−1/√3 IS the Friedmann point**. The gauge gate is background-blind — that is the true bottom of "necessary but not sufficient" (D-030).
3. **The 0.35699 footnote was the fingerprint of the true background** (corrects D-031's red-herring ruling): the true EC has N̄'(sonic)=−0.355699 (the paper's 0.35699 = dropped-digit typo).
4. **True EC built** (`nr_ec2.py`, HKA §sec-ss.practice verbatim): bisection between case1 (A<1) / case2 (second sonic) failure modes → **V₀\* = 0.112439401388**; one zero of V ✓; center-relaunch closes onto the sonic values to 6 digits ✓.
5. **β lands.** Spectrum on EC: κ_rel = 2.8105501 / 2.8105526 / 2.8105528 (N=200/300/400; edge-standoff-insensitive 1e-2..1e-3), gauge → 1 (1.0027→1.0009), next modes Re κ ≤ −1.4.

```
kappa_0 = 2.8105526(3)   (reference 2.8105525488)
beta    = 1/kappa = 0.3558019   (reference 0.35580192; G-ANCHOR |dbeta| ~ 1e-8 << 4e-3)
G-CONVERGE PASS · G-UNIQUE PASS · G-ANCHOR PASS
```

**NEXT:** `nr_shoot` on the EC background (redundant recovery; expect the sonic-gauge mode at 0.3557 as a second control) → C++ Stage-A re-freeze (supersedes `b4f4e463` — Friedmann-as-EC — with an operator-signed note) + Stage-B shoot port → β golden → `fluidcss_nexus` v1.0.0.
