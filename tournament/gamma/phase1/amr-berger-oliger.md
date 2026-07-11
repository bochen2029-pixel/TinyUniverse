# APPROACH — Berger–Oliger AMR: resolve γ the way γ was discovered

**Persona:** the AMR advocate. **Fork:** F2 (resolution strategy), with a load-bearing claim on F1 (formulation) and F3 (observable). **Stance assigned (phase 0):** *AMR is the proven route — Choptuik found γ and the echoing WITH Berger–Oliger adaptive mesh refinement; it is the reference answer, not a gamble.* **My hard job (CHARTER determinism lens):** show a Berger–Oliger scheme can be built single-file, made byte-deterministic and golden-freezable, and run inside ~5 min CPU — or concede honestly where it cannot. I do not hand-wave the determinism; it is my exposed flank and I attack it myself first.

---

## 0 · Thesis in three sentences

Critical collapse is self-similar: as `p → p*` the solution echoes down through scale factors of `e^Δ ≈ 31×` per echo, and a *fixed* grid caps the smallest feature it can hold — this is exactly why the incumbent's uniform grid saturated at `≈(φ/dr)²` and its bisected p\* wandered `0.40→0.356` under refinement (D-021, measured). Berger–Oliger AMR does not add resolution *everywhere* (unaffordable); it recursively adds fine subgrids **only where a deterministic refinement criterion fires**, buying `≥ 2^L` effective resolution at the echo core for a cost that grows like the *echo's* shrinking support, not the domain — which is precisely how Choptuik (1993) first saw γ ≈ 0.374 and Δ ≈ 3.44. My claim: a *restricted, deterministic* Berger–Oliger scheme (fixed refinement schedule, curvature-threshold flagging with fixed buffers, deterministic subcycling, no floating tie-breaks) is buildable single-file and golden-freezable, and it is the highest-probability route to a *real* γ receipt in N0's spirit — but its honest home is the deepest-strained corner of the doctrine, and it may cost the budget.

---

## 1 · The AMR scheme, fully specified (no hand-waving)

Everything below is a *concrete recipe* an implementer can turn into `substrate_amr.cpp`. It reuses the incumbent's physics verbatim — same polar-areal EKG, same constrained metric solve, same RK4, same KO dissipation, same parity/Sommerfeld boundaries (`substrate/substrate_nexus.cpp`) — and wraps it in a Berger–Oliger mesh hierarchy. **What changes is the grid, not the equations.** That keeps the oracle relationship intact: the base level *is* the incumbent, so an AMR run with zero refinement must reproduce the incumbent golden bit-for-bit (my reversibility lemma, §3.4).

### 1.1 The hierarchy
- **Base grid `L0`:** uniform polar-areal, `N0 = 400` cells on `r ∈ [0, r_max=24]`, `dr0 = 0.06`. (Half the incumbent's 800 — AMR earns its resolution locally, not globally.)
- **Levels `L1…Lmax`, `Lmax = 6`:** each level `ℓ` is a set of **patches** (contiguous cell ranges) refined from its parent by a **fixed integer refinement ratio `ρ_r = 2`** in space. Effective finest resolution `dr_fine = dr0 / 2^Lmax = 0.06/64 ≈ 9.4e-4` — a **64×** resolution gain at the echo core, reaching scales the uniform grid provably cannot (§2). Level count `Lmax = 6` is a *declared, golden-pinned* dial, not adaptive.
- **1-D simplification (the doctrine gift):** spherical symmetry ⇒ every patch is an *interval* `[r_a, r_b]`, never a 2-D box. Clustering (Berger–Rigoutsos in the general case) collapses to **interval merging** — deterministic, tie-free, `O(N)`. This is why AMR is *tractable single-file here* where it would not be in 3-D.

### 1.2 Refinement criterion (deterministic, no floating tie-breaks)
A cell `j` on level `ℓ` is **flagged for refinement** iff a **fixed dimensionless curvature proxy** exceeds a **fixed threshold**:

```
χ_j  =  dr_ℓ² · ( Π_j² + Φ_j² )           # local field-energy density × cell-area  (the Garfinkle–Duncan
                                          # subcritical scaling variable, made scale-relative)
flag(j)  ≡  χ_j  >  χ_thr                  # χ_thr a DECLARED, golden-pinned constant (e.g. 0.02)
```

Why this proxy: `Π²+Φ²` is exactly the source that drives both constraints and *is* the quantity the incumbent measured saturating at the grid ceiling (`substrate_nexus.cpp` `zmax`, line ~306). Multiplying by `dr_ℓ²` makes the criterion **scale-relative** — a feature that is "sharp for this level's resolution" flags regardless of absolute magnitude, which is the self-similar-correct trigger (each echo looks the same at its own scale). **Determinism guards:**
- `χ_thr` is a compile-time/declared constant; **no percentile, no max-relative, no "top-k"** criterion (those introduce data-dependent tie-breaks). A strict `>` with a fixed constant has **no ties to break**.
- Comparison is `fp64 χ_j > fp64 χ_thr`. Equality (`==`) is treated as **not flagged** (deterministic, documented) — the measure-zero boundary is resolved by a fixed rule, never by iteration order.

### 1.3 Buffer zones & regridding (deterministic clustering)
- **Buffer `B = 4` cells** (declared): every flagged cell dilates its flag to `±B` neighbours before clustering, so a feature cannot propagate off the fine patch between regrids (the classic Berger–Oliger safety margin; `B` chosen so `B·dr_ℓ >` the fastest characteristic can cross in one regrid interval — a CFL-tied, checkable inequality).
- **Clustering = interval union:** sort flagged-cell indices (they are already ordered), merge any two runs separated by `< 2B` gap into one patch, enforce a **minimum patch width `W_min = 8` cells** and **parent-alignment** (patch edges land on even parent-cell boundaries). All integer arithmetic on a sorted index list ⇒ **bit-identical patch layout for identical field data, every time.**
- **Regrid schedule is FIXED, not adaptive-in-time:** regrid every `K = 8` base-level steps (declared). Between regrids the hierarchy is frozen. This is the single most important determinism decision: *the regrid times are a function of the step counter alone*, never of a convergence test or a wall-clock trigger. Same `(params, seed)` ⇒ same regrid ticks ⇒ same hierarchy history ⇒ same declared hash.

### 1.4 Time subcycling (Berger–Oliger, deterministic order)
Berger–Oliger integrates finer levels with **proportionally smaller timesteps**: level `ℓ` takes `2^ℓ` substeps of `dt_ℓ = dt0 / 2^ℓ` per one base step (CFL `dt_ℓ = λ·dr_ℓ` preserved at every level — the whole point of subcycling is to keep CFL tight without globally shrinking `dt`). The recursion (Berger & Oliger 1984, Algorithm; Berger & Colella 1989 for the conservative form) is a **fixed-order tree walk**:

```
advance(ℓ, dt_ℓ):
    for substep s in 0 .. ρ_t-1:              # ρ_t = 2, fixed
        set fine ghost cells by PROLONGATION from level ℓ-1 (time-interpolated)   # deterministic
        take one RK4 step of the EKG fields on all patches of level ℓ
        if ℓ < Lmax: advance(ℓ+1, dt_ℓ/2)     # recurse, fixed pre-order
    RESTRICT level ℓ+1 solution onto the overlapped coarse cells of level ℓ        # deterministic
```

**Determinism guards:** the tree walk is **pre-order, patch-index-ordered, fp64 sequential, single-thread** — no parallel reduction, no atomics, no OpenMP. The *only* place order could bite is summation across patches; there is none (patches are disjoint intervals, updated independently). This is strictly simpler than the incumbent's RK4-with-per-substage-metric-solve, which already ships byte-deterministic.

### 1.5 Grid-transfer operators (the Berger–Oliger machinery)
- **Prolongation (coarse→fine ghost fill):** cubic Lagrange interpolation in space, linear in time (standard BO). Coefficients are **fixed rationals** (e.g. `[-1, 9, 9, -1]/16` for the mid-cell); fp64 evaluation is deterministic. Cubic (not linear) matters: linear prolongation injects a `O(dr²)` truncation kink at every refinement boundary that *looks like* a spurious feature and re-triggers flagging — a determinism-safe but accuracy-poison choice. Flag if the implementer must drop to linear for a single-file budget: `[CITE-UNVERIFIED]` on whether cubic is strictly required for γ vs a convenience.
- **Restriction (fine→coarse):** straight injection (copy the fine value at the collocated point) for the *field* variables `φ, Φ, Π`; this is the constraint-consistent choice for polar-areal EKG because the metric `a, α` is **re-solved by the outward radial integral on the composite grid**, not evolved — so there is no flux to conserve at coarse-fine boundaries, sidestepping the hardest part of Berger–Colella conservative AMR. **This is a genuine simplification the constrained scheme grants us** and a central plank of feasibility (§3.3).
- **Coarse-fine boundary flux fix (the Berger–Colella refluxing step):** *not needed for the constrained metric* (no conservative flux across the boundary — the metric is elliptic-constrained each substep, not hyperbolic-fluxed). The scalar field's characteristic flux is handled by the ghost-cell prolongation + KO dissipation at the boundary. **This is my most attackable physics claim** — I flag it (§6, R3): a naive drop of refluxing can leak a small constraint violation at coarse-fine interfaces; the honest test is whether the Hamiltonian-constraint residual at the boundary stays `< tol` under refinement (a metamorphic gate the deciding experiment must include).

### 1.6 The composite metric solve (the elegance that saves us)
Because `a(t,r), α(t,r)` are solved by **outward integration of an ODE in `r`** (not evolved in `t`), the hierarchy needs only a **single ordered radial pass over the composite grid** each substep: walk `r` from 0 to `r_max`, stepping with `dr_ℓ` wherever a fine patch is active and `dr0` elsewhere. The Heun/trapezoidal integrator (incumbent lines ~157–173) is `O(N_composite)` and **inherently sequential and deterministic** — it is the *same* integrator, just over a variable-spacing 1-D mesh. No linear solve, no iteration, no matrix. This is the single biggest reason AMR is feasible *single-file* for this problem and would not be for a 3-D BSSN evolution.

---

## 2 · Why AMR resolves γ where the uniform grid (D-021) cannot — mechanism, not assertion

D-021's corpse reason, quoted exactly: *"a fixed grid caps the self-similar central curvature (≈(φ/dr)²), the grid ceiling; refining does not converge, it makes near-critical chaotic … p\* moved 0.40→0.356 from N=800→1600."* Two failures, two AMR answers:

1. **The curvature ceiling.** Self-similarity means the interesting structure at echo `n` lives on scale `L_n ∝ e^{-nΔ} ≈ 31^{-n}`. A uniform `dr` resolves echo `n` only while `dr < L_n`; beyond that the feature falls *between grid points* and the measured `max_t(Π²+Φ²)` **saturates at `≈1/dr²`** — not physics, arithmetic (this is the incumbent's `zmax` ceiling, `substrate_nexus.cpp` ~line 306, and the reason the Garfinkle–Duncan subcritical scaling variable saturates rather than diverging as `(p*−p)^{-2γ}`). AMR **follows the echo down**: each time the criterion (§1.2) fires because a feature sharpened past the current level's resolution, a new finer level is born *at that feature*, so `dr_effective` shrinks in lockstep with `L_n`. With `Lmax=6` we hold `≥ 64×` finer scale — roughly `log(64)/Δ ≈ 1.2` extra echoes past the base grid, and adaptively *more* near p\* where they concentrate. That is the difference between seeing `~1` echo (uniform, γ invisible) and seeing enough log-periods to fit the exponent.

2. **The "refining goes chaotic" pathology.** The incumbent's uniform refinement moved p\* because a *globally* finer grid resolves *new* short-wavelength dynamics **everywhere at once**, including far-field radiation and origin-regularity error, faster than the bisection's fixed collapse criterion could track — the near-critical regime became stiff and the floating collapse flag chased noise. AMR is **surgical**: resolution is added *only at the flagged echo core*, the far field stays coarse and quiet, and — critically for the tournament — the search is handed to `autotune`'s **pre-registered level-crossing with a fixed metric** (§4, F5-steal), so p\* cannot wander with the grid. The instability D-021 measured was *partly the global refinement and partly the floating search*; AMR fixes the first, ORRERY-armed autotune fixes the second. **This is the crux of my case: D-021 killed "uniform + hand-bisection," not "adaptive + pre-registered locate."**

**Reference-answer weight (the stance).** This is not a hopeful analogy. Choptuik (1993) *discovered* γ ≈ 0.374 and the echo Δ ≈ 3.44 **using Berger–Oliger AMR on exactly this system** (massless real scalar, spherical, polar-areal-family). Gundlach & Martín-García (2007) review a decade of AMR-and-similarity confirmations converging on γ = 0.3737. AMR is the *only* method on the table with a published, reproduced γ for this exact PDE. The burden of proof for AMR is "can we fit it in the box," **not** "does it work" — that is settled physics.

---

## 3 · Guarantee-audit against the fixed constraints (the honest determinism/golden story)

The CHARTER's four constraints, each answered plainly — with the strains marked, not buried.

### 3.1 Determinism-or-it-doesn't-ship → **PASS, with a named discipline**
Every source of nondeterminism in general AMR is eliminated *by construction* here:
| general-AMR nondeterminism | this scheme's fix | residual risk |
|---|---|---|
| data-dependent regrid *times* | **fixed schedule** (regrid every `K=8` base steps) | none — pure function of step counter |
| floating-point tie-breaks in flagging | strict `>` vs a **fixed constant** `χ_thr`; `==`→not-flagged | none — no ties exist |
| clustering heuristic (Berger–Rigoutsos) ambiguity | **1-D interval merge**, integer, sorted | none — 1-D collapses the 2-D ambiguity |
| parallel reduction / atomics ordering | **single-thread, disjoint patches, no cross-patch sum** | none — no reduction in the declared path |
| patch iteration order | **pre-order, patch-index-ordered** tree walk | none — fixed traversal |
| variable-spacing metric integral | **one ordered radial pass**, same Heun integrator | none — sequential by construction |

**The golden is freezable exactly as the incumbent's is:** serialize the composite-grid declared state (the same battery metrics: p\*, m_ratio, γ_fit, R², plus a hierarchy-summary hash) → blake2b → `goldens/substrate_amr/golden.hash`, and re-run in a fresh process to confirm `GOLDEN OK`. The in-process determinism probe (incumbent's S6) extends to "two full AMR runs → identical declared object." **I assert byte-determinism is achievable and I put it on the record as the thing to attack.** Where it strains: the discipline is *heavier* — every one of the six rows above is a place a careless implementer reintroduces nondeterminism, so the code needs a determinism-selftest that a stranger runs cold (KAT-style: a fixed field snapshot → fixed patch layout → fixed hash). That selftest is part of the tool contract (§4).

### 3.2 CPU ≤ ~5 min → **STRAINED — my exposed budget flank (honest)**
Cost model. Base grid `N0=400`. AMR adds work `∝ Σ_ℓ (cells on level ℓ) × (substeps 2^ℓ)`. Near p\* the flagged region is the echo core — empirically a small fraction of the domain, but it *deepens* (more levels) as p→p\*, and subcycling multiplies work by `2^ℓ` per level. Rough envelope: if the finest levels touch `~O(100)` cells each and we run `Lmax=6` with subcycling, the per-step cost is `~ N0 + Σ_{ℓ=1}^{6} 100·2^ℓ ≈ 400 + 100·126 ≈ 1.3e4` cell-updates vs the incumbent's `800`/step — call it **~16×** the incumbent's per-step cost, but on a **coarser base and shorter effective evolution** (the collapse resolves faster once tracked). The incumbent's full battery is ~20 s; a single well-tracked near-critical AMR evolution is plausibly **30–120 s**, and the deciding experiment needs a *scaling ladder* of ~10–15 such evolutions for the fit → **plausibly 8–20 min for the full γ golden.** **This may exceed the 5-min budget.** Honest positions:
- *Best case:* the fit needs only supercritical points (mass-scaling) at moderate `Lmax`, ~10 evolutions × ~20 s ≈ **3–4 min → PASS.**
- *Likely case:* subcritical curvature-scaling (Garfinkle–Duncan) needs deeper levels near p\* → **6–12 min → the operator prices budget-vs-crown (A4 in the ledger is a real trade, exactly as predicted).**
- *This is the single place I concede the doctrine may not stretch to fit me.* I do not pretend otherwise. My mitigation: `Lmax` is a dial; a `Lmax=4` "argument-grade" golden (~2 min) can ship as the *frozen* artifact while a `Lmax=6` "citable" run is a documented longer-cold-run — the two-pass discipline the project already uses for citable claims.

### 3.3 Single-file-where-possible → **PASS (feasible), with honest weight**
The scheme is single-`.cpp`-feasible because spherical symmetry collapses AMR's two hardest general-case components: **(a)** 2-D/3-D box clustering → 1-D interval merge (~40 lines), and **(b)** conservative flux refluxing → *nothing* (constrained metric, no conservative flux; §1.5). What remains — a patch list (a `std::vector<struct Patch{int lo,hi,level;}>`), the recursive `advance()`, cubic prolongation, injection restriction, the variable-spacing metric pass — is **~400–600 lines on top of the incumbent's ~650**, all stdlib, fp64, single-thread. It is *more* code than the incumbent and the recursion is genuinely trickier to get byte-stable, but it is not a "research-grade AMR library" (no PARAMESH, no Carpet, no AMReX) — it is a *restricted, hand-rolled, doctrine-shaped* BO scheme. **Strain, stated:** this is the most complex single file in the v2 tree, and its correctness surface (grid transfers, ghost fills, regrid) is where subtle bugs hide. That is an argument for a *ferocious* selftest, not against the approach.

### 3.4 Exit codes 0/1/2 + reversibility → **PASS**
Standard envelope (0 pass · 1 declared gate fired · 2 error), never conflated — inherited verbatim from the incumbent. **Reversibility lemma (my graveyard insurance):** *set `Lmax = 0` and the composite grid degenerates to the base `L0` grid with no patches; the AMR tool must then reproduce a `N=400` uniform-grid run bit-for-bit.* This is a built-in metamorphic anchor: it proves the AMR machinery adds *only* refinement and corrupts nothing, and it means the incumbent's honest Type-II transition result is a strict *subset* of what the AMR tool computes (`Lmax=0` ⊂ `Lmax=6`). **Pre-registered reinstatement trigger** (if AMR is later killed): *if a future double-null or CSS approach resolves γ to tol at < 5 min, AMR is retired to the graveyard with corpse-reason "correct but budget-dominated"; it reinstates the moment a GPU N3 `horizon` stage exists where its `2^L` subcycling cost is affordable* — AMR's natural home may simply be N3, not N0 (the ledger's F4/A4 tension, which I concede is real).

---

## 4 · Pre-registered falsifier + costed deciding experiment

### 4.1 The tool contract (the "AMR EKG gear" — does NOT exist yet; specified here)
`substrate_amr` · v0.1.0 (proposed) · single-file C++17, fp64, stdlib, no GPU. Universal envelope: `substrate_amr [--Lmax n] [--N0 n] [--chi-thr x] [--regridK n] [--collapse-flag] [--json|--golden|--selftest]` · exit 0/1/2 · `(params)→byte-identical declared JSON→blake2b`. Faces:
- `--dev p` — one AMR evolution at amplitude `p`; prints `collapsed` (bool), `m_bh`, `max_t(Π²+Φ²)` (the un-capped curvature), hierarchy depth reached. **`--collapse-flag` makes it emit a single `{"collapsed": 0|1}` metric** — the level-crossing observable `autotune` reads to locate p\* (F5-steal, §4.3).
- `--gamma` — the crown: locate p\*, run the scaling ladder, fit γ, emit `{gamma_fit, r2, pstar, echo_delta}`.
- `--selftest` — blake2b KAT + **the determinism KAT** (fixed field snapshot → asserted patch layout → asserted hash) + **the reversibility KAT** (`Lmax=0` reproduces the frozen uniform-grid hash).
- `--golden` — freeze/check `goldens/substrate_amr/golden.hash`.

### 4.2 PASS bar (pre-registered, mechanism-precise)
**DE-γ (AMR):** build `substrate_amr` → `autotune --locate crossing` locates p\* against pre-registered target → fit `M_BH ∝ (p−p*)^γ` over the supercritical ladder →
> **PASS ≡ |γ_fit − 0.374| < 0.05 AND R² > 0.99 AND declared golden reproduces on a fresh process AND wall-clock ≤ 5 min (`Lmax` set to the deepest that fits).**

**Redundancy gate (DE-redundant, F3-steal):** on the same substrate, also recover the **echo period** `Δ = 3.44 ± 0.3` (log-periodicity of `max_t(Π²+Φ²)` vs `ln(p*−p)`) *and* the **subcritical curvature exponent** (`max_t R ∝ (p*−p)^{-2γ}`, Garfinkle–Duncan). **Redundant-recovery PASS ≡ the mass-scaling γ and the curvature-scaling γ agree to < 0.03.** Two independent observables agreeing is the oracle-honesty firewall against "grid-noise dressed as γ."

**The falsifier (what kills AMR outright):** *if `substrate_amr` with `Lmax` up to the budget ceiling produces γ_fit that (a) does not converge toward 0.374 as `Lmax` increases `4→5→6`, OR (b) whose p\* still wanders `> 0.02` between `Lmax` levels (i.e. AMR did not cure the D-021 wander), OR (c) whose Hamiltonian-constraint residual at coarse-fine boundaries grows with refinement (the missing-refluxing risk, §1.5), THEN AMR is falsified for N0 and the corpse reason is "the refinement did not converge in the affordable depth" — a strictly stronger, more honest corpse than D-021's.* This is pre-registered *before* the build; convergence-with-`Lmax` is the single graph that decides it.

### 4.3 Costed plan
| step | tool | cost | gate |
|---|---|---|---|
| 1. build `substrate_amr` | `cl /std:c++17 /O2` | ~1 day dev; compile < 30 s | compiles clean, selftest PASS |
| 2. reversibility check | `--Lmax 0 --golden` | < 30 s | reproduces incumbent uniform-grid hash |
| 3. locate p\* | `autotune --tool substrate_amr --sweep p --metric collapsed --locate crossing --level 0.5 --target <p*_expected>` | ~10–15 evolutions × 20–60 s ≈ 3–15 min | p\* located, `on_target`, cite blake2b |
| 4. γ ladder + fit | `substrate_amr --gamma` | ~10 evolutions ≈ 3–10 min | PASS bar §4.2 |
| 5. redundancy | echo Δ + curvature exponent | reuses ladder | agree < 0.03; Δ = 3.44 ± 0.3 |
| 6. freeze golden | `--golden` | < 1 min | GOLDEN OK on fresh process |

Any γ figure this experiment ultimately reports is **`[ARGUMENT-GRADE]` until step 4 runs** — the real AMR gear does not exist yet; §5's ORRERY runs prove only the *search-and-gate mechanism*, not γ.

---

## 5 · ORRERY-armed evidence (run now, no GPU) — the search-and-gate mechanism, PROVEN

The CHARTER's RESOLVABILITY lens says: *drive ORRERY to run the critical-point search and show the number, rather than assert it.* I cannot run an AMR EKG solve today (gear doesn't exist). What I **can** make evidence-grade *now* is the load-bearing claim of §2(2) and §4.3: **that a pre-registered, deterministic level-crossing locate of a critical point against a fixed target, with a real gate, works and is byte-reproducible.** This is exactly `autotune`'s proven job (it found ratchet's ρ_c = 0.2581 vs analytic 0.25 blind — CHARTER). I drove it three ways (artifacts: `experiments/amr-berger-oliger/`, reproduce with `run_autotune_locate.ps1`):

**Command (verified, CWD `C:\ORRERY`):**
```
python tools\autotune\autotune.py --objective threshold --obj-center 0.5 --obj-width 0.05 \
  --lo 0 --hi 1 --locate crossing --level 0.5 --target 0.5 --tol 0.02 --seed 0 --json
```
This models the p\* locate: a sigmoid "collapse fraction" crossing 0.5 at the critical amplitude, located against a **pre-registered target 0.5**.

| run | what it proves (maps to my claim) | result | declared blake2b |
|---|---|---|---|
| **positive** | the locate hits the pre-registered critical point exactly (§4.3 step 3) | `x_located=0.5000`, `error=0.000`, `on_target=true`, **exit 0** | `db490a3111b770996fcec05b9bc59981f35395980be5f1ad617101dd13552c80` |
| **neg-control** | the gate is *real*, not decorative — a wrong target FIRES `G-OFF-TARGET` | target 0.60 → `error=0.100`, `on_target=false`, **exit 1** | `bf7023ee399c91946e010eab1fb0bd39c7aa686e40c05a38013dd31e37329c10` |
| **tight-tol** | the locate is *sharp*, not scraping the threshold | `tol=0.005` still PASS, `error=0.000`, **exit 0** | `6c6e2adac5c7fd549e1eda4d0ba73c2466165a21061240f0760b1be9307a5131` |

**Determinism receipt:** the positive run repeated in-process yields a **byte-identical** declared object (hash `db490a31…` both times — verified in the harness). `autotune --selftest` → SELFTEST PASS (blake2b KAT + determinism + off-target-fires, all green). ORRERY version 1.0.0; autotune golden `c79002f23cf2…` (frozen, `goldens/autotune/declared.hash`).

**What this does and does not establish.** It establishes — **evidence-grade** — that the F5 search-methodology fix (pre-registered crossing with a fixed metric, replacing the incumbent's floating hand-bisection that let p\* wander) is deterministic, gated, and hash-reproducible: the *second* of D-021's two failure modes is curable *today* with a tool that exists. It does **not** establish γ; the AMR resolution claim (the *first* failure mode) is `[ARGUMENT-GRADE]` pending `substrate_amr` (§4). I do not dress the sigmoid as physics — it is a mechanism proof, labelled as such.

---

## 6 · Honest risks

- **R1 · Budget (my exposed flank).** Subcycling's `2^ℓ` cost may push the citable γ golden past 5 min (§3.2). *This is the most likely reason AMR loses on doctrine rather than physics.* Mitigation: `Lmax` dial + two-pass (short frozen golden, long citable cold-run); honest fallback: AMR's home is N3-GPU, not N0 (§3.4 reinstatement trigger).
- **R2 · Determinism surface area.** Six distinct nondeterminism sources are each *fixable* (§3.1) but each is a place a careless refactor breaks byte-identity. Higher maintenance burden than any other approach. Mitigation: the determinism KAT in `--selftest` (§4.1) makes a break loud and cold-detectable.
- **R3 · Missing refluxing / constraint leak.** Dropping Berger–Colella conservative refluxing (§1.5) is justified for the *constrained* metric but can leak a Hamiltonian-constraint violation at coarse-fine boundaries. *This is my most attackable physics claim.* Mitigation: the falsifier §4.2(c) gates constraint-residual-vs-refinement explicitly — if it leaks, AMR fails honestly, loudly, early.
- **R4 · Prolongation truncation re-triggering flags.** Linear prolongation kinks re-flag boundaries (§1.5); cubic is safer but heavier and its necessity is `[CITE-UNVERIFIED]`. Mitigation: measure it (convergence-with-`Lmax`, §4.2a).
- **R5 · Complexity vs the single-file ethos.** ~500 extra lines of the trickiest code in the tree (§3.3). Mitigation: the reversibility lemma (§3.4) bounds the blast radius — anything that breaks `Lmax=0`≡incumbent is a bug, cold-detectable.
- **R6 · F1/A1 coupling.** AMR is specified on the incumbent's *polar-areal* gauge; if the double-null advocate is right that polar slicing crowds resolution into a coordinate singularity the grid can't follow, then AMR-on-polar-areal may need *even more* depth than budgeted, and AMR-on-double-null could be the real winner. *I concede F1 and F2 may not be independent* (ledger A1). This is a steal opportunity, not a fatal wound: BO subcycling composes with *either* gauge.

## 7 · STEAL list (harvest regardless of AMR's fate)

1. **Pre-registered `autotune` level-crossing for p\*** (F5) — *proven today* (§5), deterministic, gated, hash-cited. **Every** approach should replace hand-bisection with this; it cures D-021's *search* pathology independent of the grid. **Highest-value steal in the tournament.**
2. **The scale-relative refinement proxy `χ = dr² (Π²+Φ²)`** — the self-similar-correct trigger; useful to *any* adaptive or similarity-coordinate scheme as a "have I resolved this echo?" diagnostic.
3. **The reversibility lemma as a built-in metamorphic anchor** (`Lmax=0` ≡ incumbent) — a pattern any refinement/coordinate-change approach can adopt: *your fancy method, dialled to nothing, must reproduce the honest baseline bit-for-bit.*
4. **The constrained-metric "no refluxing needed" observation** (§1.5) — if CSS or double-null also uses constrained (not free) evolution, they inherit this simplification.
5. **The determinism KAT** (fixed snapshot → fixed hash) as a selftest pattern for *any* approach whose determinism is non-obvious.
6. **The convergence-with-refinement falsifier graph** (§4.2) — the single deciding plot; the reframe advocate should demand it of *every* γ claim before the crown is awarded.

---

## 8 · Cited literature (precise; unsure ones flagged)

- **Choptuik, M. W.**, "Universality and scaling in gravitational collapse of a massless scalar field," *Phys. Rev. Lett.* **70**, 9–12 (1993). — The origin: γ ≈ 0.37, discrete self-similarity (echoing), **found with Berger–Oliger AMR** on the spherical massless real scalar. The reference answer this tournament anchors to. *(Volume 70, page 9 — matches the ledger's citation.)*
- **Berger, M. J. & Oliger, J.**, "Adaptive mesh refinement for hyperbolic partial differential equations," *J. Comput. Phys.* **53**, 484–512 (1984). — The BO algorithm: recursive refinement, time subcycling, the regrid/prolong/restrict structure §1.4 is built on.
- **Berger, M. J. & Colella, P.**, "Local adaptive mesh refinement for shock hydrodynamics," *J. Comput. Phys.* **82**, 64–84 (1989). — Conservative AMR + the flux-refluxing step at coarse-fine boundaries that §1.5 argues is *not needed* for the constrained metric (and the clustering machinery that 1-D collapses).
- **Gundlach, C. & Martín-García, J. M.**, "Critical phenomena in gravitational collapse," *Living Rev. Relativity* **10**, 5 (2007). — The review: γ = 0.3737, Δ ≈ 3.44 for the massless scalar; a decade of AMR-and-similarity confirmations. The source of the ledger's precise Gundlach value.
- **Garfinkle, D. & Duncan, G. C.**, "Scaling of curvature in subcritical gravitational collapse," *Phys. Rev. D* **58**, 064024 (1998). — The subcritical curvature-scaling observable `max R ∝ (p*−p)^{-2γ}` used as the redundant second observable (§4.2, DE-redundant) and the source of the incumbent's `zmax` diagnostic. *(Ledger cites this exact volume/article for the subcritical variable.)*
- **Garfinkle, D.**, "Choptuik scaling in null coordinates," *Phys. Rev. D* **51**, 5558 (1995). — Reports γ from a *uniform double-null* grid; cited by the ledger (F1) as the evidence the resolvability failure is gauge-linked. Relevant to R6. `[CITE-UNVERIFIED]` on the exact page/exponent value I have not independently re-verified beyond the ledger's reference.

*All non-flagged citations match the assumption-ledger's own references (Choptuik PRL 70:9; Gundlach LRR 10:5; Garfinkle–Duncan PRD 58:064024; Berger–Oliger and Berger–Colella are standard, volume/page as given). The two Garfinkle-1995 and the cubic-prolongation-necessity claims are flagged `[CITE-UNVERIFIED]` — I did not fetch primary PDFs this session; do not treat those two specifics as verified.*

---

*The spec is the product; the golden is load-bearing. AMR is how γ was found — my job was never to prove it resolves the exponent (Choptuik did), but to prove it fits the box. It fits everywhere except, honestly, the clock — and that one strain is a trade the operator prices, not a fault I hide.*
