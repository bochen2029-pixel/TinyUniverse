# TINY UNIVERSE — CONTINUATION PROMPT (handoff from the 2026-07-13 session)

*Triple-purpose, self-contained: (a) a memory dump of the session that cracked the fluid-β operator, (b) a sitrep that situates a fresh session, (c) ranked next-steps that steer the dedicated β-eigenvalue numerics. Paste this whole file into a new session. Read §0 first; trust the files over this prose; verify against git.*

---

## 0 · RECONSTITUTION POINTER — read in this order, then VERIFY vs git

1. **This file, whole.**
2. `C:\Users\user\.claude\projects\C--TinyUniverse\memory\MEMORY.md` → linked `tiny-universe-rehydration.md` (deep memory) + `local-disk-search-tools.md` + `autonomous-loop-no-polling.md`.
3. In `C:\TinyUniverse`: **`RUN_STATE.md`** (position + OPERATOR WORK QUEUE + the CURRENT DIRECTIVE block) → **`tournament/gamma/phase4/RESULTS_hka_beta.md`** (the β-track saga — **the single most important file for the active thread**; read its 2026-07-13 update sections in full) → **`DECISIONS.md`** (read **D-029, D-028, D-027, D-026, D-025, D-024, D-023, D-022, D-021, D-020**) → `ROADMAP.md` → `CLAUDE.md`.
4. **`git log --oneline -20`** = GROUND TRUTH. HEAD should be `80c4770` (or later), tree clean apart from untracked `runs/*.png|gif` (non-declared). If docs disagree with git, **git wins**; reconcile.
5. **Do NOT re-run the golden harness to "check"** — receipts are two-passed. GPU shared (preflight exits 3 if <2 GB free; CPU oracles run regardless).

**Disk beats memory. Never resume blind. The β work is CPU-only Python — no GPU needed.**

---

## 1 · Tooling (operator flagged — do NOT use grep / `Get-ChildItem -Recurse`)

- files by name/date → `python C:\everything\search.py "QUERY"` (`ext: path: dm:today size: "phrase"`).
- content in files → `C:\everywhere\build\Release\everywhere.exe [-i -n -e PAT] [-g GLOB] PATH` (the rg replacement).
- fetch papers/weights → `python C:\fetcher\fetch.py url "URL" --dest DIR` (used this session to pull the arXiv source tarball; **never `hf download`**).
- The β-track Python needs numpy/sympy/scipy (present: numpy 2.4.2, sympy 1.14.0, scipy 1.17.1, Python 3.13). **Run all β scripts from `C:\TinyUniverse\tournament\gamma\phase4`** (they use relative imports + relative `goldens/` paths).

---

## 2 · THE HEADLINE (CORE — never drop)

**The fluid-β Stage-B perturbation OPERATOR is CRACKED and verified to machine precision — this was the wall that blocked the overnight autonomous run AND this session's first pass.** The remaining task, **extracting the eigenvalue κ₀=2.81055255 → β=0.35580192**, is walled at the *regular-singular sonic point* — the known-hard core of critical-collapse numerics. **Nine distinct extraction methods were tried on the correct operator; all recover the gauge mode but wall on the physical mode.** β is **withheld, none faked** (D-016/D-021). The next session's job is the *dedicated* eigenvalue numerics on the already-verified operator — **not** another brute-force sweep (that space is exhausted).

- **The crack (the scientific crux, DONE):** fetched the authoritative long paper **Hara–Koike–Adachi `gr-qc/9607010` = PRD 59 104008 §V** (source on disk: `tournament/gamma/phase4/HKA99_src/rflanl.tex`). The transcribed coefficients (5.5–5.10) were **correct all along**; the bug was the **ASSEMBLY**. The paper writes the eigenmode system as `M_x·Ψ' = (Gmat − κ M_s)·Ψ` with **M_x ≠ identity** (its fluid block is `[[Ax,Bx],[Cx,Dx]]`). So the correct operator is **`L = M_x⁻¹(Gmat − κ M_s)`**. Tool: **`hka_pert_hka99.py`** — passes the gauge-mode exactness gate at **`|res|/|Ψ_g| ≈ 1e-10` for every κ̄** (was O(1)). Committed `fe86cb7`. Decision: **D-029**.
- **Single next action:** rebuild the sonic-point analytic-expansion machinery accurately for the correct operator (`hka_pert_sonic.py` is inaccurate for it), then read κ off the existing Evans/match-determinant machinery (`hka_beta4.Delta`) with a compound-matrix center leg to kill mode collapse. See §5.
- **Hard constraints:** contract-first · golden-gated · **honesty-over-progress (never fake β)** · two-pass citable claims · don't-grind-a-wall · don't-poll-your-own-tasks. Exit codes 0/1/2 never conflated.
- **Source of truth on disk:** `RESULTS_hka_beta.md` (β saga, full), `hka_pert_hka99.py` (the operator), the `hka_beta_*.py` scaffolds (the 9 methods).

---

## 3 · SITREP — where the whole project stands (Ring 3 background; terse, git-verifiable)

**TINY UNIVERSE** (`C:\TinyUniverse`): a scale-compressed physics universe (native Windows C/C++/CUDA), the *inhabitable* sibling of ORRERY, under the same doctrine (contract-first, golden-gated, deterministic-or-it-doesn't-ship, every claim gated against an oracle or honestly walled). Three axes → one binary: **A** substrate physics (N0→N7), **B** renderer (R0→R5), **C** rigor crowns (fluid-β, scalar-γ).

- **v1 (M0–M7): COMPLETE** — the regime ladder + world (21 goldens): beauty → gravity → arrow → relativity → black holes → quantum → cosmos.
- **v2 substrate ladder: N0→N3 CLOSED** — N0 `substrate_nexus` (spherical EKG oracle; Choptuik *transition*, γ deferred, D-021), N1 `field` (SP weld 6/6), N2 `lapse` (redshift 2/2), N3 `curve` (light-bending + precession + Shapiro 3/3). **All four classical tests of GR pass on the substrate.** Plus 2.5PN inspiral (D-025) + Q-006 resolved (D-026).
- **Axis C crowns:** fluid-β **Stage-A BANKED** (golden `fluidcss_stageA` `b4f4e463`, Evans–Coleman background, oi*=3/8, sonic point exact; D-027). Stage-B = **this session's work** (operator cracked, β walled — see §4). Scalar-γ=0.374 = not started (the DSS/Floquet crown, follows the fluid warm-up).
- **The γ-tournament ruling (D-028, reconciled this session):** the Choptuik exponent route is **CSS×spectral eigenvalue (γ=1/Re λ₀)**, NOT AMR-in-N0 (graveyarded on the ≤5-min budget → **reinstates at N4-GPU**). D-021 **amended**: the uniform-grid wall is a *gauge×discretization* artifact, not a grid law (2nd-order FD needs ~9722× the DOF of Chebyshev; `experiments/spectral/cheb_convergence`). The fluid-β warm-up validates the CSS×spectral machine on a *known* exponent before the scalar-γ claim — which is exactly why closing Stage-B matters.
- **Axis B (renderer):** least-advanced; R0 Vulkan⇄CUDA interop parked (Vulkan SDK not installed, D-002). The biggest gap to the "inhabitable" north star, but orthogonal to the β work.
- **git:** master, HEAD `80c4770`. This session added `3d1e400`→`80c4770` (7 commits: reconciliation, the crack, the 9-method machinery, D-029, the diagnoses).

---

## 4 · THE β EIGENVALUE — full technical state (Ring 1: the active thread)

### 4a · The verified operator (DONE — `hka_pert_hka99.py`)

Radiation fluid, γ=4/3. Similarity vars `s=−ln(−t)`, `x=ln(−r/t)`; reduced fields `N=α/(a eˣ)`, `A=a²`, `ω=4πr²a²ρ`, `V`=3-velocity. State vector **Ψ=(Ā_p, N̄_p, ω̄_p, V_p)** = (δlnA, δlnN, δlnω, δV). Eigenmode `h_var = h_p(x) e^{κs}`. The system (rflanl.tex §V `eq:EOM-eigenmodes`):

```
[ M_x d/dx − Gmat ] Ψ = −κ M_s Ψ         =>   L(κ) = M_x⁻¹ (Gmat − κ M_s)
M_x  = [[1,0,0,0],[0,1,0,0],[0,0,Ax,Bx],[0,0,Cx,Dx]]          (NOT the identity — the fix)
Gmat = [[G1,0,G3,G4],[H1,0,H3,0],[E1,E2,E3,E4],[F1,F2,F3,F4]]
M_s  = [[0,0,0,0],[0,0,0,0],[0,0,As,Bs],[0,0,Cs,Ds]]
```
Coefficients (all in `hka_pert_hka99.coeffs`, verbatim from §V; g=γ=4/3, oV2=1−V²):
`As=1, Bs=gV/oV2, Cs=(g−1)V, Ds=g/oV2`; `Ax=1+NV, Bx=g(N+V)/oV2, Cx=(g−1)(N+V), Dx=g(1+NV)/oV2`;
`E1=−((g+2)/2)ANV`, `E2=((6−3g)/2)NV−((2+g)/2)ANV+(2−g)ωNV−NV·ω̄_x−gN·V_x/oV2`, `E3=(2−g)ωNV`, `E4=((6−3g)/2)N−((2+g)/2)AN+(2−g)ωN−N·ω̄_x−g(1+2NV+V²)V_x/oV2²`;
`F1=((2−3g)/2)AN`, `F2=(2−g)(g−1)ωN+((7g−6)/2)N+((2−3g)/2)AN−(g−1)N·ω̄_x−gNV·V_x/oV2`, `F3=(2−g)(g−1)ωN`, `F4=−(g−1)ω̄_x−g(N+2V+NV²)V_x/oV2²`;
`G1=−A, G3=2{1+(g−1)V²}ω/oV2, G4=4gωV/oV2²`; `H1=A, H3=(g−2)ω`.
(ω̄_x=ω'/ω, V_x=V' are BACKGROUND derivatives from `bg_state(x)` = `hka_beta4.bg_state`.)

**Verification (the gate a correct operator must pass):** the pure-gauge mode `Ψ_g=(Ā_x, N̄_x+κ̄, ω̄_x, V_x)` is an exact solution for EVERY κ̄. `hka_pert_hka99.py __main__` → `|res|/|Ψ_g| ≈ 1e-10` at all x, all κ̄∈{0.357,1,2.81}; the κ¹ gauge condition `E2=As·ω̄_x+Bs·V_x`, `F2=Cs·ω̄_x+Ds·V_x` holds to ~1e-15 (it *is* the background Eq 3/4). **This is the proof the operator is correct.**

### 4b · The eigenvalue problem (what β needs)

Target (paper Table I): **γ=4/3 → κ=2.81055255 → β_BH=1/κ=0.35580192**, the **unique** relevant (Re κ>0) mode. Conditions: (i) analytic through the sonic point; (ii) regular at the center (x→−∞); gauge `N̄_p(x_s)=0`. The algebraic **identity** (`eq:alg-PP`, first integral): `(κ−A)Ā_p + cN·N̄_p + com·ω̄_p + cV·V_p = 0` with `cN=2gNVω/oV2, com=2ω{1+(g−1)V²+gNV}/oV2, cV=2gω{N(1+V²)+2V}/oV2²` (in `hka_beta_solve.identity_coeffs`). Center asymptotic modes (`eq:PPasmp1`): {(0,1,0,0), (0,0,1,0)} bounded; (0,0,0,1)e^{−2x} blows up (killed by κ); (1,−1,…,0)e^{−x} excluded by the identity. At +∞: `1, e^{−κx}, e^{−κx}` — bounded iff Re κ>0.

### 4c · The nine methods tried + why each walled (measured; all in `phase4/`)

| # | method | file | result / wall |
|---|---|---|---|
| 1 | Frobenius match det[v_center\|s1,s2,s3] | `hka_beta4.Delta` | noisy; deepest min at κ≈1.90 (spurious), not 2.81 |
| 2 | rejection (non-gauge v_center vs sonic span) | inline | floors at ~2e-3, no dip at 2.81 |
| 3 | growth shoot, 4-D (sonic→center) | `hka_beta_solve.growth` | monotonic ~1e10, no min (4-row ODE doesn't preserve the identity) |
| 4 | `solve_bvp` collocation, 4-D | `hka_beta_bvp` | converges to −1.67 / 0.98 / 1.24 (gauge/spurious), never 2.81 |
| 5 | Chebyshev spectral, 4-D | `hka_beta_spec` (orig) | only κ≈0.994 (gauge mode) |
| 6 | Chebyshev spectral, identity-augmented (id row replaces Ā ODE) | `hka_beta_spec` (edited) | degenerate κ=1 (singular B-matrix) |
| 7 | reduced-3D growth shoot | inline | 1e57 blow-up (backward integration excites the sonic +7 mode) |
| 8 | reduced-3D forward 2×2 match determinant | `hka_beta_match` | clean gauge sign-change at κ≈1.0; det EXPLODES to −1e5 in the physical band |
| 9 | Lyapunov s-evolution generator, `eig(Q)` | `hka_beta_lyap` | spectrum swamped by O(N²) spurious modes; only one O(1) real (κ≈−13.4) |

**Two structural findings that MUST carry forward:**
1. **The identity is NOT preserved by the 4-row ODE** (measured: an identity-exact vector drifts to 24% off-identity by x≈−2). ⇒ ANY 4-row shoot/collocation admits identity-violating modes and misses the physical one. **Always use the identity-REDUCED 3-D system** (Ā_p eliminated via the identity → `hka_beta_match.L3`).
2. **This background's gauge puts the gauge mode at κ≈1, not 0.357.** (The paper fixes the sonic at x=0; our `hka_ec` background has the sonic at x_s≈−0.144 with `N̄'(x_s)≈−1`, so `N̄_p(x_s)=0 ⇒ κ̄=−N̄'(x_s)≈1`.) Every method finds this gauge mode; **discard it** — the physical mode is 2.81.

### 4d · The precise blocker (the sonic point) + the broken Frobenius machinery

The sonic point (x_s≈−0.144) is a **regular singular point**: `M_x` fluid block `det(AxDx−BxCx)→0`, so `L~R/(x−x_s)`. Consequences that defeat direct shooting/collocation: (a) a **violent +7 eigenvalue** of the reduced operator near x_s (backward integration blows up; forward match determinant explodes); (b) the **Frobenius series diverges for |t|≳0.3** (ratio test = 1.0; order-35 coeffs overflow); (c) **the existing sonic machinery `hka_pert_sonic.py` is INACCURATE for the new operator** — `bg_series_near_sonic` fits the background near x_s with a **low-order polyfit** of sampled points; the residue/indicial is not even pinned: DFT `L_laurent`→μ≈−4.48, an analytic `R=(1/d′)·adj(fluid)·(Gmat−κM_s)_fluid`→μ≈+10.5, old-code "expectation" 1−2κ=−4.62 — **three different values.** The physical Frobenius series has O(1)–O(100) ODE-residual at every |t|∈[0.05,0.20]. So the sonic analytic modes (needed by every match/Evans method) are currently unreliable.

### 4e · Files / scaffolds (all committed, all in `tournament/gamma/phase4/`)

- **`hka_pert_hka99.py`** — the VERIFIED operator (`coeffs`, `Lnum`, `_bg_derivs`, the gauge gate). **This is the anchor; everything else feeds it or reads it.**
- `HKA99_src/rflanl.tex` — the authoritative long paper (§V perturbation operator at `\subsection{Perturbation equations…}` ~line 4973; `eq:EOM-var`, `eq:EOM-eigenmodes`, `eq:alg-PP`, `eq:PPasmp1`). `KHA95_src/9503007.tex` = the PRL (background EOM only).
- `hka_beta_match.py` — reduced-3D system `L3` + forward 2×2 match determinant (**the best-designed scaffold**; `L3`, `lift`, `v_sonic` 3×3 physical-a₀ construction). `hka_beta_solve.py` — `identity_coeffs`, `v_sonic`, `Lx`. `hka_beta_lyap.py` — the Lyapunov generator `spectrum(N)`. `hka_beta_spec.py`, `hka_beta_bvp.py`, `hka_beta_validate.py` (G-ANCHOR/CONVERGE gates), `hka_beta_eigfn.py`.
- **Existing (overnight) machinery to reuse:** `hka_beta4.py` (`bg`, `bg_path`, `bg_state`, `v_center`, `Delta` = the det match), `hka_frobenius.py` (`analytic_modes` = ker(R) + Taylor recursion), `hka_pert_sonic.py` (`bg_series_near_sonic`, `L_laurent` — **needs the rebuild**), `hka_ec.py` (the verified EC background: `bg_state`, `A_of`, `_fluid_slopes`, `center_state3`, `shoot_to_sonic`). Monkeypatch idiom: `import hka_pert_core as PC; PC.Lnum = hka_pert_hka99.Lnum` routes ALL the old machinery through the correct operator.
- Contract: `contracts/fluidcss_nexus.contract.md` v0.9.1 (Stage-B is scoped there: β=1/Re κ₀, gates G-ANCHOR |β−0.3558|<4e-3 / G-CONVERGE / G-UNIQUE — **do not change the contract**; complete the module). Tool `substrate/fluidcss_nexus.cpp` (Stage-A C++; β golden pending).

---

## 5 · RECOMMENDATIONS + STEERING — the dedicated β numerics (do this, in order)

**Framing:** the operator is right; β is bounded numerics on it. Do NOT re-run the 9-method sweep. Pick ONE rigorous route and engineer it properly. **Primary = fix the sonic Frobenius, then Evans/compound-matrix.** Alt = worked time-evolution Lyapunov. Both are CPU-only Python, ≤ a focused session each.

### 5a · PRIMARY PATH — rebuild the sonic analytic expansion, then compound-matrix Evans (recommended)

1. **Rebuild `hka_pert_sonic` accurately for the new operator.** The current polyfit background is the root inaccuracy. Compute the **background sonic Taylor series analytically**: at x_s the *desingularized* flow (`hka_desing.py`) is regular, so `(A,N,ω,V)(t)` has an exact Taylor series in `t=x−x_s` — get its coefficients from the desingularized RHS recursion (or a very-high-accuracy `solve_ivp` from x_s on the desingularized flow + Chebyshev fit, tol 1e-13). Then the Laurent `L(t)=R/t+L0+L1 t+…` follows. **Reconcile the residue**: the analytic `R = (1/d′)·[[Dx,−Bx],[−Cx,Ax]]·(Gmat−κM_s)[fluid rows]` (embed as rows 2,3; `d′=d/dx(AxDx−BxCx)|_{x_s}`) vs the DFT `L_laurent` — they must agree, and the indicial exponent (the non-zero eig of R) must be a *clean, κ-consistent* value. **Gate: the rebuilt Frobenius series must have ODE-residual < 1e-8 at |t|=0.1** (the current fails at O(1)). This unblocks every match/Evans method.
2. **Compound-matrix (exterior-algebra) Evans function on the reduced-3D system** to kill the mode-collapse that broke the shooting. For the 3-D reduced ODE (identity built in, `hka_beta_match.L3`): the center-regular subspace is 2-D → evolve its **2nd compound** (the Λ² minors, a 3-vector for a 3×3 system) from the center forward — the compound evolves **stably** (no collapse) under the induced operator `L3^{(2)}`. The Evans function `E(κ)` = the wedge of the forward center-regular compound with the sonic-analytic compound (from the rebuilt Frobenius) at a match point; its zeros (away from the gauge κ≈1) are the eigenvalues. This is the textbook cure for exactly the stiffness that walled methods #3/#7/#8. *(Note: the overnight run reportedly built a compound-matrix Evans solver — search `everywhere.exe -e "compound|evans|wedge" phase4` for `hka_beta_evans.py` or similar; with the correct operator it may largely work.)*
3. **Read κ, gate, freeze.** With an accurate Frobenius, try `hka_beta4.Delta` (monkeypatched to `PC.Lnum=hka_pert_hka99.Lnum`) FIRST — it may now resolve κ=2.81 cleanly (its noise in method #1 was the inaccurate Frobenius). Bisect the Evans/Delta zero near 2.81, discard the gauge zero near 1.0. Validate: G-UNIQUE (only one relevant Re κ>0), G-CONVERGE (κ stationary under grid/order refinement), redundant recovery. Then **port to `fluidcss_nexus.cpp` Stage-B, freeze the β golden, bump the contract to v1.0.0**, write D-030, update RUN_STATE/MODULE/ARCHITECTURE §11.

### 5b · ALT PATH — worked time-evolution Lyapunov (HKA §V "second method", "very regular")

The `eig(Q)` shortcut (method #9) fails on O(N²) spurious modes, but the **actual time evolution** (power iteration) is robust — the physical growing mode dominates and high-k spurious modes are dissipated. Evolve `∂_s(ω̄,V) = Msinv[(E·Ψ,F·Ψ) − Mx_fl ∂_x(ω̄,V)]` in s (RK4 or implicit) from random data, with the metric constraints solved each step (`∂_xĀ=G·Ψ` BC Ā(x_c)=0; `∂_xN̄=H·Ψ` BC N̄(x_s)=0). Requirements the naive version lacks: **(i) sonic-aware upwinding** of the ∂_x term (characteristic-split on the sign of the fluid characteristic speeds, which change sign at x_s); **(ii) proper BCs** (center-regular; +∞ outflow — the modes are `e^{−κx}` there); **(iii)** measure Re κ from `d/ds ln‖Ψ‖` after the transient → β=1/Re κ. Generator matrices already assembled in `hka_beta_lyap.spectrum` (reuse the constraint-elimination + `Q`). This avoids the sonic ODE singularity entirely and is the paper's own method.

### 5c · Honest fallback (pre-registered, D-028/hybrid.md)

If both routes prove out-of-scope in a focused session: the tournament's **`reframe`/boson-star floor** ships as the crown and scalar-γ (and the fluid-β number) defer — *nothing wasted, the operator crack stands as a real result*. **Never fake β to close it** (D-016/D-021). The scalar-DSS γ=0.374 crown (`similarity_nexus`, Gundlach gr-qc/9604019) follows a *validated* warm-up — do not start it before β lands.

**Steering, one line:** rebuild the sonic Frobenius accurately (the true blocker), then let the compound-matrix Evans / `hka_beta4.Delta` read κ off the correct operator — that is the shortest honest path to β.

---

## 6 · Discipline + hard-won lessons (carry these)

- **The RAYFORMER lesson, twice-proven:** measure/read against what the object *actually is*, not its nominal label. It flipped Q-006 (D-026); and this session it flipped the β wall — two operators "failed the gate" and it *looked* like a coefficient/derivation bug for two sessions, but it was a **matrix-structure reading** of the primary source (M_x≠I). When stuck, fetch and read the primary source verbatim.
- **Honesty over progress (D-016/D-021):** β is NOT measured and NOT faked. A faked eigenvalue poisons the oracle farm. The wall is named with its exact resume state.
- **Don't grind a wall:** this session tried 9 methods — that was warranted under "go all out," but the space is now *exhausted*; the next move is a *different engineered method* (§5), not a 10th sweep. (Operator flagged grinding/polling twice — see `autonomous-loop-no-polling.md`.)
- **Contract-first, golden-gated, deterministic.** Each `_nexus` tool is a single self-contained file (blake2b + envelope + golden). CPU oracles build with `cl /std:c++17`; GPU with `nvcc -arch=sm_89` (no fast-math). Registration per milestone = ARCHITECTURE §11 + TASKLIST + DECISIONS + RUN_STATE + MODULE, then atomic commit.

---

## 7 · Cadence

Operator works directive-by-directive ("go", "keep going", "proceed") and this session said **"go all out"** — expects deep, thorough pushes, wants to *see* results (SendUserFile for visuals; the β work is numeric so there was none), and quota can cut mid-flight → commit early, keep RUN_STATE current, atomic-commit each close. Suspend and ask plainly if the operator returns confused/frustrated. This session's operator explicitly recommended **banking the operator crack and giving β the dedicated careful treatment** — honor that: open the next session by *scoping* the §5a rebuild, not by re-sweeping.

---

## 8 · VERBATIM TAIL (most recent last; if replaying, recall backwards)

> **Operator (with a screenshot showing context 36%, 5-hour limit 0%):** "you have plenty of context window, see screenshot, go for it, go all out:"
> **Me → (worked):** fetched `gr-qc/9607010`, found M_x≠I, built `hka_pert_hka99.py`, gauge gate → 1e-10. Committed `fe86cb7`. Then attacked the eigenvalue.
> **Me (mid-push report):** "🎯 The wall is broken … the corrected operator passes the gauge gate at |res|/|Ψ| ≈ 1e-10 … Closing on β: switched from shooting to collocation … the sweep now hunts the unique relevant mode (Re κ>0), κ=2.81055255 → β=0.35580192."
> **[then 9 methods, all walling on the sonic point; the identity-non-preservation + reduced-3D findings; the Frobenius-machinery inaccuracy for the new operator (residue μ: −4.48 vs +10.5 vs −4.62); the Lyapunov generator swamped by O(N²) spurious modes.]**
> **Me (final report):** "The wall is broken; the number is not (and I won't fake it). … nine distinct methods … all extract the gauge mode but wall on the physical κ=2.81 … the known-hard core … Choptuik/Maison/HKA/Gundlach each solved this with dedicated engineered schemes … The hard scientific crux is done. … I'd genuinely recommend banking this here … the β number deserves the careful dedicated treatment the literature gives it rather than another round of the sweep I've now exhausted. If you want, I can scope that dedicated effort (the specific compound-matrix Evans or worked-Lyapunov implementation) as the next session's plan."
> **Operator (this trigger):** "yes give me a continuation prompt for a new session that triples as a current memory dump/save of your state and also as a sitrep … and also your recommendations of the next steps and guiding and steering it towards a certain path … very rigorous, detailed and comprehensive and self-contained full context prompt. go"
> **Me:** wrote this file (`docs/CONTINUATION_2026-07-13.md`), updated the memory pointer, committed.

*Build one rung right. Freeze its golden. Prove its number — or name the wall exactly and hand off clean. The operator is cracked; go land β with the engineered method, not the sweep.*
