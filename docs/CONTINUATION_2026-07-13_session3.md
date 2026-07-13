# TINY UNIVERSE — CONTINUATION / SITREP (session 3 handoff, 2026-07-13)

*Written 2026-07-13 08:46 CST (Monday). Triple-purpose, self-contained: (a) a memory dump of the
session that **PROVED the fluid-β perturbation operator correct from first principles** and isolated
the β wall to the eigenvalue-extraction method, (b) a sitrep situating a fresh session, (c) ranked
next-steps + roadmap. Paste this whole file into a new session. Read §0 first; trust files over prose;
verify against git. Supersedes `docs/CONTINUATION_2026-07-13.md` for the β thread.*

---

## 0 · RECONSTITUTION POINTER — read in this order, then VERIFY vs git

1. **This file, whole.**
2. `C:\Users\user\.claude\projects\C--TinyUniverse\memory\MEMORY.md` → linked `tiny-universe-rehydration.md` + `local-disk-search-tools.md` + `autonomous-loop-no-polling.md`.
3. In `C:\TinyUniverse`: **`RUN_STATE.md`** (position + OPERATOR WORK QUEUE + CURRENT DIRECTIVE) → **`tournament/gamma/phase4/RESULTS_hka_beta.md`** (the β saga — **the single most important file for the active thread**; read the 2026-07-13 *session 3* sections in full) → **`DECISIONS.md`** (read **D-031, D-030, D-029, D-028, D-027**) → `ROADMAP.md` → `CLAUDE.md`.
4. **`git log --oneline -12`** = GROUND TRUTH. HEAD should be `73eb99f` (or later); tree clean apart from untracked `runs/*` (non-declared visuals), `tournament/gamma/phase4/HKA99_src/*` (the arXiv paper source), and `*.pkl` (gitignored caches). If docs disagree with git, **git wins**; reconcile.
5. **The β work is CPU-only Python** (numpy/sympy/scipy; Python 3.13). No GPU needed. **Run all β scripts from `C:\TinyUniverse\tournament\gamma\phase4`** (relative imports + relative `goldens/`).

**Disk beats memory. Never resume blind.**

---

## 1 · Tooling (operator-flagged — do NOT use grep / `Get-ChildItem -Recurse`)

- files by name/date → `python C:\everything\search.py "QUERY"` (`ext: path: dm:today size: "phrase"`).
- content in files → `C:\everywhere\build\Release\everywhere.exe [-i -n -e PAT] [-g GLOB] PATH` (the rg replacement).
- fetch papers/weights → `python C:\fetcher\fetch.py url "URL" --dest DIR` (**never `hf download`**).
- Anchor recency with the system clock: `Get-Date -Format "yyyy-MM-dd HH:mm:ss K (dddd)"`.

---

## 2 · THE HEADLINE (CORE — never drop)

**The fluid-β Stage-B perturbation OPERATOR is now PROVEN CORRECT FROM FIRST PRINCIPLES — this corrects
the prior session's (and my own D-030's) suspicion of an operator error.** The β number
(κ₀=2.81055255 → β=0.35580192) is **still not landed, none faked** (D-016/D-021), but the wall is now
**precisely isolated to the eigenvalue-EXTRACTION numerics** (the known-hard core of critical-collapse),
NOT the operator, background, or identity — all three are verified. The next move is HKA's *second*
method, the **Lyapunov time-evolution** (§V.G), which sidesteps the sonic singularity entirely.

### What is PROVEN (banked, this session)
- **The operator is faithful to first principles.** `nr_rederive.py` symbolically linearizes the primary
  self-similar EOM (`rflanl.tex` 4.1/4.2, keeping ∂ₛ→κ) and re-derives the fluid operator coefficients
  FROM SCRATCH: **all 16 (As,Bs,Cs,Ds, Ax,Bx,Cx,Dx, E1–E4, F1–F4) match `hka_pert_hka99` (=§V) EXACTLY,
  diff=0.** Metric rows (G,H) verified by hand from 4.3/4.4. ⇒ §V *is* the linearized EOM; my operator
  *is* §V. **D-030's "residual operator error" inference is FALSIFIED (D-031).**
- **The sonic Frobenius machinery is exact + gated** (the prior session's named blocker, now fixed):
  - `nr_sonic.py` — analytic background sonic Taylor series. **Key structural finding:** the order-1 fluid
    solvability is **quadratic** in the null-coefficient α₁ (two analytic branches through the sonic
    point); the Evans–Coleman branch is picked via the desingularized-flow eigenvector (`hka_desing`).
    Validated reference-free: ODE-residual ~1e-15; matches the integrated background to **3.7e-11** at
    t=−0.02. Radius ≈0.12 (match inside |t|≤0.05).
  - `nr_laurent.py` — EXACT operator Laurent by power-series arithmetic (no DFT). **Residue μ = 1−2κ
    EXACTLY** (−4.621105 for κ=2.81), reconciling the prior three-value disagreement (−4.48/+10.5/−4.62).
    All **3 analytic modes** solve the DIRECT operator ODE to **9.8e-14** at t_m=−0.02. Frobenius gate PASSES.
- **The background is the true Evans–Coleman solution.** `N ∝ e^{−x}` (N̄'≡−1) and `A = 1+(2−γ)ω` hold
  EXACTLY along it — confirmed by BOTH `hka_ec` integration AND the independent `nr_sonic` branch
  selection (A_k=(2−γ)ω_k, N_k=N₀(−1)ᵏ/k! to ~1e-13). Sonic values (3/2, 2/√3, 3/4, −1/√3), center
  (NV=−2/(3γ)=−1/2), and EOM (4.3/4.4) all match §IV verbatim.
- **A validated eigenvalue reader.** `nr_shoot.py` implements HKA's *own* method (`rflanl.tex` §V,
  5.452–5.463): build the unique analytic-at-sonic solution (ker R + identity eq:alg-PP + gauge
  N̄_p(0)=0, norm Ā_p(0)=1), integrate **sonic→center**, kill the lone `e^{−2x}` expanding mode (pure
  V_p). Positive control: it resolves the **gauge mode to κ=0.99999924** (= −N̄'(x_s) = 1).

### What is WALLED (the active problem)
- **κ=2.81055255 is ABSENT from the sonic→center shoot.** Across κ∈[−1.5, 12] real + the complex
  neighborhood, SIX independent methods (`nr_shoot` 4-D & reduced-3D, `hka_beta_solve`, `hka_beta_match`,
  `nr_evans2`, complex scan) surface ONLY the gauge mode (κ=1). Since the operator + background + identity
  are all PROVEN correct and the sonic-data conditioning is fine at 2.81 (cond≈10, Ā-content≈0.80), this
  is a genuine **known-hard-numerics** failure of the *shoot*, not an upstream error. **I could not
  explain it by inspection.** HKA warn this method "is almost impossible to search every eigenvalue."

### The red herring (do NOT chase)
- **The 0.35699 gauge-mode footnote** (`rflanl.tex` §V.B footnote) states the sonic-gauge gauge mode is
  at κ=−N̄'(0)≃0.35699. This is IMPOSSIBLE with the EC sonic values: the EOM (4.4) forces
  N̄'(0)=−2+A₀−(2−γ)ω₀ = −2+1.5−0.5 = **−1.0** (⇒ gauge κ=1, which the reader finds); 0.35699 would
  require ω₀<0 (unphysical, since ω=4πr²a²ρ≥0). I burned real effort chasing it as a coordinate/variable
  scaling; it is a dead end (likely a paper typo or a different-solution/convention quantity).

---

## 3 · SITREP — where the whole project stands (Ring 3; terse, git-verifiable)

**TINY UNIVERSE** (`C:\TinyUniverse`): a scale-compressed physics universe (native Windows C/C++/CUDA),
the inhabitable sibling of ORRERY, under the same doctrine (contract-first, golden-gated,
deterministic-or-it-doesn't-ship, every claim gated against an oracle or honestly walled). Three axes:
**A** substrate physics (N0→N7), **B** renderer (R0→R5), **C** rigor crowns (fluid-β, scalar-γ).

- **v1 (M0–M7): COMPLETE** — the regime ladder + world (21 goldens): beauty → gravity → arrow →
  relativity → black holes → quantum → cosmos.
- **v2 substrate ladder: N0→N3 CLOSED** — N0 `substrate_nexus` (spherical EKG oracle; Choptuik
  *transition*, γ deferred, D-021), N1 `field` (SP weld 6/6), N2 `lapse` (redshift 2/2), N3 `curve`
  (light-bending + precession + Shapiro 3/3). **All four classical tests of GR pass on the substrate.**
  Plus 2.5PN inspiral (D-025) + Q-006 resolved (D-026).
- **Axis C crowns:** fluid-β **Stage-A BANKED** (golden `fluidcss_stageA` `b4f4e463`, Evans–Coleman
  background, oi*=3/8, sonic point exact; D-027). Stage-B = **this thread** (operator PROVEN correct;
  β walled at extraction — §4). Scalar-γ=0.374 = not started (DSS/Floquet crown; `similarity_nexus`,
  Gundlach gr-qc/9604019) — follows the fluid warm-up.
- **The γ-tournament ruling (D-028):** the Choptuik-exponent route is **CSS×spectral eigenvalue
  (γ=1/Re λ₀)**, NOT AMR-in-N0 (graveyarded on the ≤5-min budget → reinstates at N4-GPU). D-021 amended:
  the uniform-grid wall is a gauge×discretization artifact. The fluid-β warm-up validates the machine on
  a KNOWN exponent before the scalar-γ claim — which is why closing Stage-B matters.
- **Axis B (renderer):** least-advanced; R0 Vulkan⇄CUDA interop parked (Vulkan SDK not installed, D-002).
  The biggest gap to the "inhabitable" north star, but orthogonal to β.
- **git:** master, HEAD `73eb99f`. This session added `2c68cc5`→`73eb99f` (7 commits): exact background
  series, exact Laurent, the eigenvalue reader + κ=2.81-absent finding, D-030, first-principles operator
  proof, D-031.

**OPERATOR WORK QUEUE (from RUN_STATE):** #3 v1-polish (2.5PN ✅ + Q-006 ✅ + galaxy-orbit ✅;
GARGANTUA Kerr art pass remains) → **#4 crowns (ACTIVE — fluid-β Stage-B: operator ✅ proven, β extraction
walled; scalar-γ follows)** → #2 renderer (Vulkan SDK gap) → #1 AMR contract → N4-GPU.

---

## 4 · THE β EIGENVALUE — full technical state (Ring 1: the active thread)

### 4a · The verified operator (`hka_pert_hka99.py` — PROVEN, do not re-derive)
Radiation fluid γ=4/3. Similarity vars `s=−ln(−t)`, `x=ln(−r/t)` (β_scale=1). Log fields
`Ā=lnA, N̄=lnN, ω̄=lnω`, `V` linear. State `Ψ=(Ā_p, N̄_p, ω̄_p, V_p)`. Eigenmode `h_var=h_p(x)e^{κs}`.
System (`rflanl.tex` eq:EOM-eigenmodes): `M_x Ψ' = (Gmat − κ M_s)Ψ` ⇒ **`L(κ)=M_x⁻¹(Gmat − κ M_s)`**,
with `M_x = diag-block(1,1,[[Ax,Bx],[Cx,Dx]])`, `M_s = block(0,0,[[As,Bs],[Cs,Ds]])`. All coefficients in
`hka_pert_hka99.coeffs`, **verified diff=0 vs the first-principles linearization** (`nr_rederive.py`).
Gauge mode `Ψ_g(κ̄)=(Ā_x, N̄_x+κ̄, ω̄_x, V_x)` is an exact solution ∀κ̄ (gauge gate `|res|/|Ψ|`≈1e-10).
Identity (eq:alg-PP, verified `idc·Ψ_g≈1e-16`): `(κ−A)Ā_p + cN N̄_p + com ω̄_p + cV V_p = 0`
(`hka_beta_solve.identity_coeffs`; `cN=2γNVω/(1−V²)`, `com=2ω[1+(γ−1)V²+γNV]/(1−V²)`,
`cV=2γω[N(1+V²)+2V]/(1−V²)²`).

### 4b · The eigenvalue BVP (HKA §V, the "first method" = the shoot)
Conditions (5.216–5.340): (i) analytic ∀x∈ℝ; (ii) regular at center (Ā_p=0 at x→−∞). At the sonic point,
analyticity (a₀∈ker R) + identity + gauge N̄_p(0)=0 make the analytic solution a **unique 1-param(κ)**
family (norm Ā_p(0)=1). Integrate sonic→−∞; the only growing center mode is `(0,0,0,1)e^{−2x}` (pure V_p;
the `e^{−x}` mode is identity-excluded); regularity ⇒ **kill it** ⇒ discrete κ. At +∞ modes are
`1, e^{−κx}, e^{−κx}` — bounded ∀ Re κ>0 (not an extra discrete condition). Target: **κ=2.81055255,
β_BH=1/κ=0.35580192, the unique relevant (Re κ>0) mode.**

### 4c · The wall (measured, this session)
`nr_shoot.py` implements 4b faithfully on the exact Frobenius + verified identity, sign-consistent
(norm Ā_p(0)=1). **Positive control passes** (gauge κ=0.99999924, log|y|→0). But the growing-mode
coefficient is **smooth and monotonic through κ=2.81 with NO zero-crossing**; a wide fine scan
(κ∈[−1.5,12], step≤0.02 near 2.81) + a complex-κ scan (Re∈[2.6,3.0], Im∈[0,2.5]) show **only the gauge
dip at κ=1**. `hka_beta_match` (reduced-3D determinant) and `hka_beta_solve` (independent shoot)
corroborate. The amplification e^{−2x_c}≈e^{26}≈1e11 makes eigenvalue dips razor-sharp, but the gauge
dip *is* caught, and there is genuinely no feature near 2.81. **Why the shoot misses 2.81 despite the
proven-correct operator is unexplained — this is the open question.**

### 4d · Files / scaffolds (all committed, `tournament/gamma/phase4/`)
**This session (the `nr_` family — the new clean tools):**
- `nr_sonic.py` — EXACT background sonic Taylor series (`bg_series(order)` → A,N,om,V coeff arrays).
- `nr_laurent.py` — EXACT operator Laurent (`laurent_exact`, `analytic_modes`, `gate`; residue μ=1−2κ).
- `nr_frob.py` — Frobenius helper (bgser dict, `residue_report`, `gate`).
- `nr_rederive.py` — **THE operator proof** (symbolic linearization of 4.1/4.2 vs §V; all diff=0).
- `nr_shoot.py` — **THE validated eigenvalue reader** (HKA's shoot; `psi_sonic`, `E`, `refine`; nails gauge to 1e-7).
- `nr_evans.py`/`nr_evans2.py`/`nr_evans_diag.py` — two-sided Evans experiments (gauge-robust; flat, superseded by the shoot).
- `nr_match.py`, `nr_scan.py`, `nr_bgcheck.py`, `nr_check2.py`, `nr_diag.py` — scans + validation drivers.
**Prior machinery still valid:** `hka_pert_hka99.py` (operator), `hka_ec.py`/`hka_desing.py` (background),
`hka_beta_solve.py` (identity_coeffs + shoot scaffold), `hka_beta_match.py` (reduced-3D L3 + lift),
`hka_beta4.py` (bg_state/bg_path, v_center), `hka_frobenius.py`/`hka_pert_sonic.py` (OLD DFT/polyfit path — superseded by nr_laurent/nr_sonic; keep for reference only).
**Primary source on disk:** `HKA99_src/rflanl.tex` (HKA gr-qc/9607010 = PRD 59 104008). Key line anchors:
background EOM 4.3/4.4 (≈4341–4351), sonic values 4.7–4.9 (≈4719–4731), perturbation operator eq:EOM-var
(≈4991), coefficients (≈5039–5088), identity eq:alg-PP (≈5192), center asymptotics eq:PPasmp1 (≈5300),
shoot method (≈5452), gauge-mode footnote/0.35699 (≈5502), Lyapunov method §V.G (≈5651).
**Contract:** `contracts/fluidcss_nexus.contract.md` v0.9.1 (Stage-B scoped: β=1/Re κ₀, gates
G-ANCHOR |β−0.3558|<4e-3 / G-CONVERGE / G-UNIQUE — **do not change**; complete the module). Tool
`substrate/fluidcss_nexus.cpp` (Stage-A C++; β golden pending).

---

## 5 · RECOMMENDATIONS + ROADMAP — the dedicated β numerics (do this, in order)

**Framing:** the operator/background/identity are PROVEN correct; the shoot won't surface κ=2.81. The
literature-standard cure is HKA's **own second method**, which they call *"very regular."* This is the
recommended primary path.

### 5a · PRIMARY — implement the Lyapunov time-evolution (HKA §V.G) — recommended
The insight that makes it robust: **M_s is NON-singular** (det M_s,fl = γ[1−(γ−1)V²]/(1−V²) ≠ 0 anywhere),
so the s-evolution has **NO sonic singularity** (unlike the shoot). Structure:
- **Fluid (ω̄,V):** `∂ₛ(ω̄,V) = M_s,fl⁻¹ [ Gmat_fl·Ψ − M_x,fl ∂ₓ(ω̄,V) ]`. M_x,fl is finite (singular only
  in the sense det→0 at the sonic point, but it MULTIPLIES ∂ₓ here — it is not inverted).
- **Metric (Ā,N̄):** solved from the SPATIAL constraints each s-step — `∂ₓĀ = Gmat[row0]·Ψ`
  (BC Ā(x_c)=0, center-regular), `∂ₓN̄ = Gmat[row1]·Ψ` (BC N̄=0, gauge). (These are eq:EOM-var rows 1–2
  with M_s=0, M_x=I on the metric block.)
- **Evolve** from generic (ω̄,V)(x) data on a grid x∈[x_center, x_sonic(+a bit)]; the **dominant growth
  rate `d/ds ln‖Ψ‖` after the transient = the relevant Re κ = 2.81** (it dominates because it is the
  largest-Re eigenvalue; gauge κ=1 is subdominant). β=1/Re κ.
- **Requirements the naive `eig(Q)` attempt (`hka_beta_lyap.py`, the overnight "method #9") lacked, and
  which are why it drowned in O(N²) spurious modes:** (i) **sonic-aware upwinding** of the ∂ₓ term
  (characteristic-split on the sign of the fluid characteristic speeds, which flip at the sonic point);
  (ii) proper BCs (center-regular + gauge); (iii) **power iteration / actual time-stepping** (RK4 or
  implicit) — NOT the raw generator eigendecomposition — so high-k spurious modes dissipate and the
  physical growing mode dominates. Reuse `hka_beta_lyap.spectrum(N)` for the constraint-elimination
  scaffold but replace `eig(Q)` with time-evolution.
- **Gate + freeze:** measure Re κ, check G-UNIQUE (one relevant mode), G-CONVERGE (κ stationary under
  grid/order refinement), redundant recovery (β vs the C++ Stage-A). Then **port to
  `fluidcss_nexus.cpp` Stage-B, freeze the β golden, bump the contract to v1.0.0**, write the decision,
  update RUN_STATE/MODULE/ARCHITECTURE §11.

### 5b · SECONDARY — crack the shoot (diagnostic, may reveal a subtle method bug)
Since the operator is PROVEN correct, κ=2.81 *must* be a bounded sonic→center solution — so `nr_shoot`
has a subtle bug I did not find. Worth a fresh look: (1) construct the reduced-3D sonic data DIRECTLY
from the reduced residue R3 = R4[1:4]·T(0) (rather than reducing the 4-D `psi_sonic`), and re-shoot;
(2) verify the reduced L3 is self-consistent (row-0 Ā' from L4 vs d/dx of the identity-slaved Ā);
(3) reconsider whether the "kill e^{−2x}" detection needs the e^{−x} identity-drift subtracted in the
4-D case. If found, β falls out of the *already-built* `nr_shoot` with no new machinery.

### 5c · Honest fallback (pre-registered, D-028)
If both prove out-of-scope in a focused session: the tournament's `reframe`/boson-star floor ships as the
crown and scalar-γ (and the fluid-β number) defer — nothing wasted; the **operator proof + exact
Frobenius stand as real results**. **Never fake β** (D-016/D-021).

### 5d · Broader roadmap (after β lands)
1. **Scalar-DSS γ=0.374 crown** — new `similarity_nexus` contract (Gundlach gr-qc/9604019, period-Δ
   Floquet spectrum), built ONLY after the fluid-β warm-up proves the CSS×spectral machine on a known
   number. Redundant recovery: eigenvalue γ vs echo Δ=3.4453 agree <0.03.
2. **Axis B renderer** — R0 Vulkan⇄CUDA interop (needs Vulkan SDK install; D-002) → the biggest gap to
   the inhabitable north star.
3. **N4 `horizon` (GPU)** — the Choptuik crown / AMR contract reinstatement (D-028): 2^L Berger–Oliger
   subcycling affordable on GPU; the fully-specified recipe is `tournament/gamma/phase1/amr-berger-oliger.md`.
4. **v1 polish tail** — GARGANTUA Kerr geodesic art pass (the last #3 item).

**Steering, one line:** the hard scientific crux (operator + background + Frobenius) is DONE and PROVEN;
implement HKA's Lyapunov power-iteration (their "very regular" second method) on the verified operator —
that is the shortest honest path to β. Don't chase the 0.35699 footnote; don't re-run the shoot sweep.

---

## 6 · Discipline + hard-won lessons (carry these)

- **RAYFORMER, twice this session.** (1) D-030 *inferred* "operator error" from the shoot's failure;
  *measuring* (the first-principles re-derivation, `nr_rederive`) FALSIFIED that — the operator is right.
  (2) The 0.35699 footnote *looked* like a coordinate-scaling clue for an hour; *computing* N̄'(sonic)
  from the EOM (=−1.0, forced) killed it. **Measure against what the object IS; keep skepticism on the
  method, not the (verified) object.**
- **Honesty over progress (D-016/D-021):** β is NOT measured and NOT faked. A faked eigenvalue poisons
  the oracle farm. The wall is named with its exact resume state.
- **Don't grind the same wall:** the shoot sweep is exhausted (6 methods, wide+complex scans). The next
  move is HKA's *different* method (Lyapunov), not a 7th sweep. (Operator flagged grinding/polling —
  see `autonomous-loop-no-polling.md`.)
- **Contract-first, golden-gated, deterministic.** Each `_nexus` tool is a single self-contained file
  (blake2b + envelope + golden). CPU oracles `cl /std:c++17`; GPU `nvcc -arch=sm_89` (no fast-math).
  Registration per milestone = ARCHITECTURE §11 + TASKLIST + DECISIONS + RUN_STATE + MODULE, atomic commit.

---

## 7 · Cadence

Operator works directive-by-directive ("go", "keep going", "proceed", "go all out"), wants to *see*
results (SendUserFile for visuals; the β work is numeric so there is none), and quota can cut mid-flight
→ commit early, keep RUN_STATE current, atomic-commit each close. Suspend and ask plainly if the operator
returns confused/frustrated. This session's operator said **"yes do that and proceed"** (pursue the
operator verification + land β); I verified the operator (a real result) but did NOT land β, and
checkpointed at the verified milestone because the finding changed the strategy (operator ✓ → wall =
extraction) and the Lyapunov method is a substantial fresh implementation.

---

## 8 · VERBATIM TAIL (most recent last; if replaying, recall backwards)

> **Operator (input to this session):** pasted `docs/CONTINUATION_2026-07-13.md` (the prior handoff:
> "operator cracked via M_x≠I, gauge gate 1e-10; β walled at the sonic numerics after 9 methods; next =
> rebuild the sonic Frobenius + compound-matrix Evans").
> **Me (session 3):** rebuilt the sonic Frobenius EXACTLY (`nr_sonic`, `nr_laurent`; μ=1−2κ, gate
> 9.8e-14). Built the eigenvalue reader (`nr_shoot`, HKA's method). Found κ=2.81 ABSENT (6 methods);
> inferred a residual operator error (D-030).
> **Operator:** "yes do that and proceed" (verify the operator beyond the gauge mode, land β).
> **Me:** re-derived the operator from first principles (`nr_rederive`) → **all 16 fluid coeffs match §V
> exactly, diff=0** → the operator is PROVEN correct; D-030's "operator error" FALSIFIED (D-031). Proved
> the background is the true EC solution (N∝e^{−x}, A=1+(2−γ)ω exact). Showed the 0.35699 footnote is a
> red herring (needs ω0<0). Isolated the wall to the eigenvalue-EXTRACTION method. Recommended HKA's
> Lyapunov "second method" as next. Committed `2c68cc5`→`73eb99f`.
> **Operator (this trigger):** "check date and time, write a sitrep and comprehensive status update
> including memory dump and your own best recommendation of what is next and roadmap pathways etc, as
> markdown md to disk."
> **Me:** wrote this file (`docs/CONTINUATION_2026-07-13_session3.md`), updated the memory pointer + git.

*Build one rung right. Freeze its golden. Prove its number — or name the wall exactly and hand off clean.
The operator is PROVEN; go land β with HKA's Lyapunov method, not the exhausted shoot.*
