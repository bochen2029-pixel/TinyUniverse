# Stage B UPDATE 2 — robust eigenvalue solvers built; corrected operator has NO physical mode at κ=2.81 (2026-07-12)

Follow-up to `STAGE_B_UPDATE.md` (v1) and `RESULTS_hka_beta.md`. v1 landed the CORRECTED, gauge-exact
perturbation operator (`hka_pert_derive.py`, gauge residuals ~1e-9). This pass built the robust
eigenvalue machinery the operator was "armed and waiting" for — and reached a **precisely-diagnosed
honest wall**: five independent solvers agree the corrected operator has **only the gauge mode (κ=1)** and
**no physical eigenvalue anywhere near κ=2.81**. No β is faked; the diagnosis is localized to the operator's
one physically-undetermined coupling.

## What was built (all runnable cold)
- **`hka_beta_evans.py` — COMPOUND-MATRIX EVANS FUNCTION (primary, definitive).** The gold-standard
  method for exactly this problem: propagate the WEDGE (2-form in Λ²ℝ⁴, 6-dim) of the admissible
  solutions by the induced compound flow L⁽²⁾, and form the Evans scalar E(κ)=U(x_m)∧W(x_m) (Λ⁴≅ℝ).
  Because the compound flow grows at SUMS of exponents, it is **immune to the exponential mode
  amplification** (the t^{1−2κ}≈t^{−4.6} sonic mode + strong non-normality, ‖[L,Lᴴ]‖≈4‖L‖) that defeats
  naive two-sided shooting. Zeros of E = eigenvalues.
- **`hka_beta_fdbvp.py` — global finite-difference BVP generalized eigenproblem** (2nd-order centered
  differences; `A0 Ψ = κ B Ψ`; center-kill + sonic-analyticity + gauge-fix BC rows; self-consistent w).
- **`hka_beta_spectral.py` — global Chebyshev collocation** GEP (Approach A of the brief) + σ_min sweep.

## The decisive result — the Evans function (validated, delta-robust)
`python hka_beta_evans.py` (xc=−14, δ=1e-3, x_m midpoint) scanning |E(κ)| over κ∈[0.5,4.0]:
```
 kappa      |E|     |E|/median
  1.00   6.94e-02     0.229   <-- DIP = the GAUGE mode (κ=1)   [correctness check: PASS]
  ...     3.02e-01    ~1.000  (flat baseline everywhere else)
  2.81    3.03e-01     1.001  (NO dip)
```
- **The method is validated:** it cleanly finds the gauge mode at **κ=1** (|E| drops 4.4×). In this EC
  background `(lnN)′≡−1` **identically**, so HKA's "sonic gauge" and "origin gauge" **coincide at κ=1**
  (the md's fn15 "κ≈0.35699" is a garbled secondary-transcription artifact; (4.1b) gives (lnN)′|₀=−1 ⇒
  N̄_p(0)=0 ⇔ κ̄=1 exactly). So the gauge cross-check is present and correct.
- **κ=2.81 is NOT an eigenvalue:** |E(2.81)| sits on the baseline and is **delta-ROBUST** — identical to
  4 s.f. across δ∈{4e-3, 1e-3, 2e-4} (an eigenvalue would drive |E|→0 as δ→0). Flat over κ∈[2.0,3.5]
  even at deep x_c=−18.

## Cross-method agreement (five independent solvers, same conclusion)
| method | finds κ=1 gauge | finds κ≈2.81 physical |
|---|---|---|
| compound-matrix **Evans** (`hka_beta_evans`) | **yes** (|E| dip 4.4×) | **NO** (flat, δ-robust) |
| Chebyshev GEP (`hka_beta_spectral.gep_eigs`) | yes | no (drifts to κ=1 as N↑) |
| σ_min collocation sweep (`hka_beta_spectral.find_eigs`) | yes (σ~4e-6) | no |
| principal-angle two-sided shoot | yes (cos₁≡1) | no (2nd angle never →1) |
| FD-BVP GEP (`hka_beta_fdbvp`) | yes | **ghost** at ~2.79 that is NON-convergent (κ↓ with M: 2.77→2.76→2.73, increments GROW) and strongly δ-dependent (δ:3e-3→1e-3 shifts κ 2.63→2.74) ⇒ a discretization artifact, not a true eigenvalue |

The FD-BVP ghost near 2.79 is seductive (close to the target) but fails every refinement test — a textbook
example of why the honest gate is refinement-stability, not proximity to the published number.

## The diagnosis — the operator is gauge-exact but its physical spectrum is wrong
The corrected operator (`hka_pert_derive.py`) satisfies **all three structural necessary conditions**:
- gauge-mode exact (residual ~1e-9, invariant in κ̄),
- center indicial eig(L)={−2,−1,0,0} (κ-independent),
- sonic residue eig(R)={0,0,0,1−2κ} exactly.

Yet it has no physical eigenvalue. **Why:** the gauge mode fixes only 6 of the 9 s-coupling (Q) entries;
`Q·e_N=0` and `Q g₀ + P e_N = 0`. After ALSO imposing the correct center exponents ({−2,−1,0,0}) and sonic
exponents ({0,0,0,1−2κ}), the residual freedom is a **1-parameter family** (call it t2; the second putative
dof t3 is FORCED to 0 by the center exponent −2). The derive picked t2=0. **No member of this
structurally-valid family (any t2, t3=0) produces a physical eigenvalue** (Evans |E| stays at baseline for
all tested t2∈[−3,8]) — so the missing physics is NOT recoverable from the gauge mode + indicial structure
alone. It is the genuine ∂_s content of HKA's original time-dependent fluid PDEs (KHA eq. 18), which the
available transcriptions do not pin down reliably:
- `hka_pert_derive.py` assumes ∂_s attaches to the self-similar-reduced (4.1d,e) via HKA's (5.5) matrix →
  gauge-exact + correct exponents, but **no physical mode** (this pass).
- The transcribed HKA (5.5)-block Q, or its sign flip, gives the WRONG sonic exponent (1+2κ, analytic) —
  broken.
- A **first-principles** re-derivation from ∇_aT^{ab}=0 keeping ∂_s (`derive_linear.py`→`fluid_sderiv.pkl`,
  whose ∂_s→0 limit reproduces the verified background) yields `hka_pert_firstprinc.py`: the Ā,N̄ rows and
  **center** exponents come out correct, but the sonic residue is rank-0 ({0,0,0,0}) — the raw covariant
  T^{aτ},T^{ar} basis hides the sonic singularity that lives in HKA's specific (4.4) combination. Needs the
  same 2×2 recombination HKA uses; not closed here.
- The equations.md eq(18) transcription (flagged `[VERIFY-BRACES]`) gives correct **sonic** exponents but
  **wrong center** ({−1,−0.5,0,0}) and fails the gauge gate — the flagged brace ambiguity in the 2nd fluid
  equation is real and blocks a clean derivation.

## Honest verdict
```
κ (physical, radiation-fluid CSS):  NOT REPRODUCED convergently.
   Evans function (amplification-immune, validated on the κ=1 gauge mode) shows the corrected operator
   has ONLY the gauge mode κ=1 and NO eigenvalue at κ=2.81 (|E(2.81)| δ-robust, on baseline).
β = 1/Re κ:  NOT EMITTED — no convergent physical κ to invert. (Target 0.35580192 NOT claimed; NOT tuned.)
Gauge-mode cross-check:  PASS — κ=1 origin/sonic gauge (|E| dips 4.4×); confirms the solvers are correct.
Stage A background:  LANDED (unchanged) — oi*=3/8, sonic=(3/2,2/√3,3/4,−1/√3) EXACT, (lnN)′≡−1, A=1+2ω/3.
```
This is a **precisely-diagnosed NEAR-MISS/FAIL**: the eigenvalue *method* is now robust and validated (Evans
function finds the gauge mode and is refinement/δ-stable); the *operator* is gauge-exact with correct
indicial structure but is missing one physically-undetermined coupling (the true ∂_s content of the KHA
fluid PDEs), so it does not carry the physical mode. A number that would be built on it — or extracted from
the non-convergent FD-BVP ghost — is not emitted.

## To close Stage B (next action, unchanged in spirit, sharper in target)
Obtain the correct s-coupling t2* from the primary source: transcribe **KHA gr-qc/9503007 eq. (18)** (the
two fluid PDEs with ∂_s) directly from the PDF, resolving the flagged brace in the 2nd equation, OR complete
the ∇_aT^{ab}=0 first-principles route (`hka_pert_firstprinc.py`) by applying HKA's (4.4) 2×2 recombination
so the sonic singularity appears. Then re-run `hka_beta_evans.py`: the physical κ is the non-gauge zero of
E(κ); β=1/Re κ; freeze the golden.

## Files (this pass)
- `hka_beta_evans.py` — compound-matrix Evans function (PRIMARY). `python hka_beta_evans.py`.
- `hka_beta_fdbvp.py` — FD-BVP generalized eigenproblem (shows the non-convergent ghost, honestly).
- `hka_beta_spectral.py` — Chebyshev collocation GEP + σ_min sweep (Approach A).
- `hka_pert_firstprinc.py` — first-principles operator from ∇_aT^{ab}=0 (Ā,N̄ + center correct; sonic
  residue not closed — documents the recombination gap).
- `derive_linear.py` (pre-existing) → `fluid_sderiv.pkl`; `arbiter_divT.py` (pre-existing) — ground-truth
  conservation checks.
```
κ, β: NOT MEASURED convergently (operator missing one physical ∂_s coupling; Evans-validated: only κ=1 gauge).
Stage A: oi*=3/8, sonic=(3/2,2/√3,3/4,−1/√3) EXACT. Gauge cross-check PASS (κ=1). No β faked or tuned.
```
