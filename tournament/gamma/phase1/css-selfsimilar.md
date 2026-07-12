# PHASE 1 — APPROACH: CSS / self-similar coordinates (the eigenvalue crown)

**Persona:** CSS / self-similar-coordinates advocate · **Forks attacked:** F1 (coordinates), F2 (resolution) · **Steel-thread:** A2 (the load-bearing premise "resolving γ requires abandoning determinism").
**Thesis.** Don't *add* refinement — **change coordinates so the self-similar thing stops shrinking.** In similarity coordinates `τ = −ln(T*−t)`, `x = r/(T*−t)`, the critical solution is **stationary in τ** (a fixed point for CSS; a period-Δ limit cycle for the DSS scalar), so a *coarse fixed grid* resolves it and determinism is trivial. And γ need not be fit from noisy near-critical collapse data at all: **γ = 1/Re λ₀**, the reciprocal of the single unstable Lyapunov exponent of the critical solution's linear perturbation spectrum (Koike–Hara–Adachi). This is the strongest determinism-and-honesty story on the board — but it is bought with a real cost: you must find T\* and *construct/target the critical solution itself*, and for the massless scalar that solution is **DSS, not CSS** — a subtlety I flag loudly and treat as my chief risk.

**Status of every quantitative γ figure below:** `[ARGUMENT-GRADE]` unless carried by a cited ORRERY blake2b. The real CSS tool (`css_echo`) **does not exist yet**; §4 specifies its contract; the only things *measured* here are (a) the literature anchors and (b) two ORRERY `autotune` critical-point locates (§Arming) that prove the *search/crossing machinery* the deciding experiment rides on. Read-only use of ORRERY throughout.

---

## 1 · The approach, fully specified

### 1.1 The similarity-coordinate reduction

Spherical massless scalar + GR. Introduce a comoving similarity time and space around the accumulation event `(T*, r=0)`:

```
τ = −ln(T* − t)        (τ → +∞ as t → T*⁻ ; each unit of τ = one factor e of zoom)
x = r / (T* − t) = r·eᵗ  (self-similar radius; the imploding structure sits at x = O(1))
```

Write the metric and scalar in scale-invariant variables. In double-null form `ds² = −e^{2σ} du dv + r² dΩ²` with self-similar coordinate `ξ = ln(−v/u)` (Gundlach's choice, gr-qc/9604019), the field is `φ(τ,x)` (or `φ(τ,ξ)`) and the metric functions become **functions of the similarity variables alone plus their τ-dependence.** The Einstein–Klein–Gordon PDEs, which in `(t,r)` describe an ever-sharpening pulse, become in `(τ,x)` a system whose **critical solution has no explicit τ-growth of scale** — the sharpening is entirely absorbed by the coordinate transform. That is the whole trick: `∂/∂(ln r)` structure that a fixed `(t,r)` grid cannot follow becomes `∂/∂x` structure at `x = O(1)` that a **fixed x-grid resolves at constant cost forever.**

**Two universality classes, and the honest distinction I must not blur:**
- **CSS (continuous self-similarity):** the critical solution is **τ-independent** — a genuine *fixed point*. The reduced equations become a set of **ODEs in x** (an eigenvalue/boundary-value problem). This is the *radiation-fluid* case (Evans–Coleman 1994; Koike–Hara–Adachi 1995).
- **DSS (discrete self-similarity):** the critical solution is **periodic in τ with period Δ ≈ 3.4453** — a *limit cycle*, not a fixed point: `Z(τ+Δ, x) = Z(τ, x)` for the state vector Z. **This is the massless-scalar case** (Choptuik 1993; Gundlach 1997). The reduction is then a **PDE periodic in τ** (one spatial dim x + the periodic τ), not an ODE. My persona's name says "CSS"; the *target physics* is DSS. I own this in §3.

### 1.2 How the critical (DSS) solution is obtained / targeted — two constructions

**(A) Direct construction (Gundlach route).** Solve the reduced τ-periodic PDE as a **two-point-in-τ boundary-value problem**: demand `Z(τ+Δ,·) = Z(τ,·)` (periodicity) + regularity at the center `x=0` + regularity at the past self-similarity horizon (the sonic/light point where the reduced equations are singular). Δ is an *eigenvalue* of this BVP, solved for jointly with the profile via a shooting/relaxation iteration. Output: the echo period Δ and the universal profile `Z*(τ,x)`. **This is a self-contained ODE/PDE solve — it never touches near-critical collapse data at all**, which is exactly why it sidesteps D-021.

**(B) Dynamical targeting (Garfinkle–Gundlach route).** Evolve the *full* EKG system in similarity coordinates from near-critical initial data. Because the critical solution is a codimension-1 attractor (one unstable mode), a tuned initial amplitude lands the evolution *onto* the limit cycle and it **echoes in place on the fixed x-grid** for many periods before the one unstable mode ejects it. You read Δ off the log-periodicity directly. This still needs a p\* tune, but the *resolution* problem is gone — the echo does not shrink off the grid, it recurs at fixed x. Route (B) is the cheaper first build; route (A) is the rigorous crown.

### 1.3 The two routes to γ (the point of this approach)

**Route I — log-periodicity of the echo (the Δ route).** In τ, the DSS solution repeats every Δ. Any self-similar observable `q` obeys `q(τ) = q(τ + Δ)`; equivalently in `(t,r)` the critical solution is invariant under `(t−T*, r) → (e^{−Δ}(t−T*), e^{−Δ} r)`. **Measure Δ** as the log-period of the central-field oscillation `φ(t,0)` vs `−ln(T*−t)` (a Fourier/autocorrelation peak in τ). This is qualitative-cheap: the echo shows up at *far coarser* resolution than any exponent fit, because you are detecting a **periodicity**, not resolving a diverging curvature. PASS target: **Δ = 3.44 ± 0.3.**

**Route II — the perturbation eigenvalue (the γ = 1/Re λ₀ route).** Linearize the reduced equations about the critical solution `Z*`. The perturbation spectrum `{λ_n}` governs how a nearby evolution approaches/departs `Z*`:
```
Z(τ,x) = Z*(τ,x) + Σₙ cₙ e^{λₙ τ} δZₙ(x).
```
Criticality ⇔ **exactly one** mode with `Re λ > 0` (the codimension-1 attractor). Call it λ₀. Because `e^{λ₀ τ} = (T*−t)^{−λ₀}`, the growing mode sets the scale at which the black-hole/disperse decision is made, and the standard scaling argument (Koike–Hara–Adachi PRL 74 (1995) 5170; Gundlach 1997) gives

```
   γ = 1 / Re λ₀.
```

- **CSS case (fixed point):** `{λₙ}` are eigenvalues of a **linear ODE eigenvalue problem** in x — a matrix generalized-eigenvalue solve after discretization. Cheap, deterministic, tiny.
- **DSS case (limit cycle, our scalar):** `{λₙ}` are the **Floquet exponents** of the period-Δ orbit — eigenvalues of the **monodromy (period-return) operator** `M` obtained by linearly propagating perturbations once around the cycle (τ → τ+Δ). `λ₀ = (1/Δ) ln(ρ₀)` where ρ₀ is the leading Floquet multiplier `>1`. Gundlach 1997 did exactly this and got **γ = 0.374 ± 0.001** with a single growing mode.

**Why this route is the honest crown:** it extracts γ from a **single eigenvalue of a well-posed linear operator on a coarse grid**, not from the slope of grid-quantized black-hole masses near a wandering p\*. There is no diverging curvature to cap (D-021's corpse), no chaotic bisection (the p\*: 0.40→0.356 pathology), and the whole computation is a deterministic linear-algebra problem — the ideal shape for a golden freeze. **Redundant recovery is built in:** Δ (Route I) and γ (Route II) are *independent observables of the same solution*, and they must be mutually consistent through the DSS structure — exactly the two-observable agreement DE-redundant demands.

---

## 2 · Why this beats the polar-areal incumbent D-021 killed

D-021's corpse (measured, honest): on a *uniform polar-areal `(t,r)` grid*, the central curvature `R ~ (Π²+Φ²)` at the origin grows without bound as `p→p*`, and the grid **caps** it at ≈`(φ/dr)²` — I verified this in the incumbent source (`substrate/substrate_nexus.cpp`: `dr = rmax/N` uniform, `zmax` tracks "peak near-central curvature (Π²+Φ²)", the grid ceiling). Refining `dr` doesn't converge; it makes the near-critical regime chaotic (bisected p\* moved 0.40→0.356, N=800→1600). The effective exponent came out ≈0.24, grid-limited, not 0.374.

**CSS deletes the mechanism of that failure rather than fighting it:**

| D-021 failure mode (uniform polar-areal) | CSS / self-similar resolution |
|---|---|
| Central curvature `~(φ/dr)²` diverges → grid **caps** it | In `(τ,x)` the structure sits at `x=O(1)` **forever**; curvature in similarity variables is **bounded and τ-periodic** — nothing to cap. |
| Refining `dr` → chaotic near-critical p\* | Route II never simulates near-critical collapse; γ is a **linear eigenvalue** of the reduced operator — deterministic by construction. |
| p\* bisection wanders (search instability, F5) | Route A needs **no p\* at all** (pure BVP construction). Route B still tunes p\* but the observable (echo period at fixed x) doesn't move under the tune. |
| AMR (F2's brute-force answer) fights the doctrine: nondeterminism, library weight, unbounded runtime | CSS is **non-adaptive** — a coarse *fixed* grid. It is *more* deterministic than the incumbent, not less. |

**This is the crux of my stance and my answer to assumption A2:** D-021 implicitly equated "no AMR" with "uniform `(t,r)` grid." CSS is a *third option* — non-adaptive **and** able to see γ — which makes A2 ("resolving γ requires abandoning determinism") **false as stated.** I refute the corpse by **deleting the need for refinement, not by adding it.** The graveyard's reinstatement bar ("a genuinely new argument, not 'just use more points'") is met: this is not more points, it is *different coordinates in which the points never run out.*

---

## 3 · Guarantee-audit against the fixed constraints (with explicit risk flags)

| Constraint (CHARTER §fixed) | CSS verdict | Where I'm at risk |
|---|---|---|
| **Determinism** — (params)→byte-identical blake2b | **STRONG.** Route II is a deterministic linear-algebra eigenvalue solve; Route A is a deterministic relaxation. No RNG, no adaptive branching. fp64 single-thread → identical hash. | Iterative eigensolvers (Arnoldi) can have iteration-order sensitivity; mitigate by pinning a **fixed-shift inverse-iteration** for the single dominant mode (deterministic) rather than a general QR on a large matrix. |
| **Golden-freezable** | **STRONG.** Declared state = `(Δ, Re λ₀, γ, R²_periodicity)` → blake2b. Small, stable, reproducible. | The declared vector must exclude iteration-count/wall-clock; standard TU discipline. |
| **≤ 5 min CPU, single-file** | **STRONG** for Route II (a few hundred x-points; a monodromy solve over one Δ = milliseconds; eigen-extract = trivial). Route A relaxation is seconds. | The **T\* determination** and the **regularity conditions at the past self-similarity horizon** are the numerically delicate part; if the shooting is stiff it could need many iterations. Budget risk is *iteration count of the BVP*, not grid size. |
| **Honest oracle** (anchor + metamorphic + redundant) | **STRONG.** Anchor: γ=0.374, Δ=3.4453 (Gundlach 1997). Metamorphic: γ invariant under x-grid refinement (it's an eigenvalue, should converge, not drift — the *opposite* of D-021). Redundant: Δ (Route I) ⟂ γ (Route II), two observables of one solution. | If the monodromy is under-resolved I could get a *plausible-but-wrong* λ₀; guard with a grid-convergence gate on λ₀ itself (must be stationary under x-refinement to <1%). |

### The three named risks the prompt demanded I flag

- **R1 — Finding T\*.** The similarity coordinates are *centered on T\**, which is not known a priori. Route A finesses this (the BVP is autonomous in τ; T\* only sets the origin, it drops out of Δ and λ₀). Route B needs T\* from the tuned collapse (read off when the central lapse/curvature accelerates). **Risk:** a mis-estimated T\* smears the log-periodicity in Route B; mitigate by *fitting* T\* as the value that maximizes the τ-periodicity sharpness (a 1-D optimization — itself an `autotune` argmax job).
- **R2 — Constructing the DSS solution / the eigenvalue problem.** This is the real cost, stated honestly: building the reduced τ-periodic EKG system, imposing regularity at center **and** at the self-similarity horizon (a coordinate-singular point requiring a series expansion / L'Hôpital treatment), and solving the eigenvalue BVP, is **genuine numerical-relativity work** — more subtle than the incumbent's method-of-lines. Gundlach spent a paper on it. It is single-file-*able* (Gundlach's was ~ODE-solver-scale) but not a weekend triviality. **This is my single biggest risk (see §6).**
- **R3 — DSS vs CSS.** The elegant "ODE eigenvalue problem" pitch is *literally true only for CSS* (radiation fluid). The **massless scalar is DSS**, so the eigenvalue problem is a **Floquet/monodromy** problem over the period-Δ cycle — one dimension harder than a pure ODE, and it presupposes you already have the periodic solution. I will **not** oversell the scalar case as a static ODE solve. The honest framing: *CSS is the clean prototype (build it first as a warm-up oracle against Evans–Coleman β≈0.36); the scalar crown is the DSS/Floquet extension.*

---

## 4 · Pre-registered falsifier + costed deciding experiment

### 4.1 The falsifier (pre-registered, before the tool is built)

> **FALSIFIER.** If a fixed-x-grid CSS/DSS solver, at a resolution that fits the ≤5-min CPU budget, **cannot** produce a leading perturbation eigenvalue that yields `|γ − 0.374| < 0.05` **and/or** an echo period `Δ = 3.44 ± 0.3`, with the eigenvalue *stationary* under x-grid refinement (metamorphic gate, drift <1%) — then the "self-similar coordinates make γ resolvable at fixed cost" claim is **KILLED**, and CSS joins the graveyard beside uniform polar-areal.
> **Reinstatement trigger (house discipline):** CSS may be re-proposed if (a) a working DSS Floquet construction is demonstrated elsewhere at CPU cost, or (b) the CSS radiation-fluid prototype nails β≈0.36 but the *scalar* extension is what failed — in which case the fluid oracle still ships as a lesser crown and only the scalar-γ claim stays dead.

### 4.2 The tool contract (the thing that does not exist yet)

**`css_echo`** — single-file C++17, fp64, stdlib-only, no GPU. Universal ORRERY envelope.

```
css_echo --mode {construct|floquet|target}
         --field {scalar|fluid}          # fluid = CSS prototype (Evans-Coleman); scalar = DSS crown
         --nx N            # fixed similarity-radius grid points (coarse: 256-1024)
         --delta-guess D   # initial echo period for the periodic BVP (scalar/DSS only)
         --tstar T         # accumulation time (target mode) or 'fit'
         --seed S --json|--csv PATH --selftest --golden
```

**Declared output vector** (→ blake2b): `{ delta, re_lambda0, gamma := 1/re_lambda0, r2_periodicity, n_unstable_modes }`.
**Gates:**
- `G-GAMMA` fires if `|gamma − 0.374| ≥ 0.05` (declared PASS bar).
- `G-DELTA` fires if `|delta − 3.4453| ≥ 0.3`.
- `G-UNIQUE` fires if `n_unstable_modes ≠ 1` (the codimension-1 signature — a physics honesty check: if two modes grow, it isn't the critical solution).
- `G-CONVERGE` fires if `re_lambda0` moves >1% between `--nx N` and `--nx 2N` (the metamorphic/anti-D-021 gate: an eigenvalue must converge, not drift).
Exit `0` all-green · `1` a gate fired (real negative result) · `2` error. Determinism gate: run twice in-process, require byte-identical declared vector.

**Oracles it names:** CSS-fluid mode → **Evans–Coleman β ≈ 0.36** (PRL 72 (1994) 1782) as the analytic warm-up anchor; scalar-DSS mode → **Gundlach γ = 0.374 ± 0.001, Δ = 3.4453 ± 0.0005** (PRD 55 (1997) 695).

### 4.3 The deciding-experiment run + PASS bar

**DE-γ(CSS):**
1. `css_echo --field fluid --mode construct` → warm-up: recover Evans–Coleman β≈0.36 (proves the reduced-equation machinery on the *easier* CSS fixed point). PASS: `|β−0.36|<0.03`.
2. `css_echo --field scalar --mode construct --delta-guess 3.44` → build the DSS solution, then `--mode floquet` → leading Floquet multiplier → γ.
3. **PASS BAR (pre-registered, the crown):** `|γ − 0.374| < 0.05` at `G-UNIQUE` satisfied (exactly one growing mode) **AND** `Δ = 3.44 ± 0.3`, `G-CONVERGE` green (eigenvalue stationary under x-refinement), deterministic golden frozen, ≤ 5 min CPU.
   *Redundant recovery (DE-redundant):* Route I (Δ) and Route II (γ) agree through DSS structure; if only one of {γ, Δ} lands, it is a WOUNDED result, not a PASS.

**Costed estimate `[ARGUMENT-GRADE]`:** grid 512 x-points; one monodromy propagation over Δ ≈ few·10³ substeps; inverse-iteration for the single dominant Floquet mode ≈ tens of iterations. Est. **< 30 s CPU** for `floquet`; the `construct` BVP relaxation dominates and is the schedule risk (R2). If `construct` blows the 5-min budget, Route B (`--mode target`, read Δ off log-periodicity of a tuned in-place echo) is the fallback that still delivers Route I.

---

## 5 · Cited literature (verified this session; `[CITE-UNVERIFIED]` where I could not confirm a specific figure)

- **Choptuik, PRL 70 (1993) 9**, "Universality and scaling in gravitational collapse of a massless scalar field." Discovered the Type-II transition; **γ ≃ 0.37**, discretely self-similar echoing with period **Δ = 3.4453 ± 0.0005** (`e^Δ ≈ 30`). *(Δ and DSS confirmed via the Semantic Scholar record + Gundlach 1997 restatement.)*
- **Evans & Coleman, PRL 72 (1994) 1782**, "Observation of critical phenomena and self-similarity in the gravitational collapse of radiation fluid." The **CSS** (continuously-self-similar) critical solution obtained by **integrating ODEs**; **β ≈ 0.36**. The direct precedent for "construct the self-similar solution as an ODE problem, don't fit collapse data." *(Verified: PRL DOI 10.1103/PhysRevLett.72.1782; arXiv gr-qc/9402041.)*
- **Koike, Hara & Adachi, PRL 74 (1995) 5170**, "Critical behavior in gravitational collapse of radiation fluid: a renormalization group (linear perturbation) analysis." Established the **eigenvalue route**: `β = (Re κ)⁻¹`, the reciprocal of the largest Lyapunov/perturbation exponent of the CSS solution; numerical **β ≈ 0.3558** *for the radiation fluid*. **Precision flag:** this paper is the *radiation-fluid* case (its value is 0.3558, **not** 0.374) — it is the *method*'s origin, not the scalar number. *(Verified: PRL DOI 10.1103/PhysRevLett.74.5170; arXiv gr-qc/9503007; β≈0.3558 confirmed.)*
- **Gundlach, PRD 55 (1997) 695**, "Understanding critical collapse of a scalar field" (arXiv gr-qc/9604019). Constructs the **DSS scalar** critical solution directly; finds **exactly one growing perturbation mode → γ = 0.374 ± 0.001**, **Δ = 3.4453 ± 0.0005**. **This is the paper that carries the eigenvalue route to the massless scalar** and is the exact template for `css_echo --field scalar`. *(Verified: PRD DOI 10.1103/PhysRevD.55.695; γ=0.374±0.001 and single growing mode confirmed.)*
- **Gundlach & Martín-García, Living Rev. Relativity 10 (2007) 5**, "Critical phenomena in gravitational collapse." The authoritative review of CSS/DSS, the eigenvalue/anomalous-dimension picture, and γ=1/Re λ₀. *(Citation as given; the 2007 LRR volume/article number is standard — `[CITE-UNVERIFIED]` on the exact page/article-number formatting, content unambiguous.)*
- *(Adjacent, for the double-null resolvability contrast in F1:)* **Garfinkle, PRD 51 (1995) 5558** — γ from a *uniform double-null* grid — cited by the ledger; not re-verified here. `[CITE-UNVERIFIED]` on volume/page.

---

## Arming — ORRERY runs (evidence-grade, read-only, this session)

The **real** CSS tool does not exist yet, so every γ number above is `[ARGUMENT-GRADE]`. What I *can* make evidence-grade **now** is the **critical-point / eigenvalue-crossing LOCATE machinery** the deciding experiment rides on: `autotune`'s job is to sweep a parameter and locate a level-crossing against a *pre-registered target* with a gate — precisely the operation Route II performs when it locates the eigenvalue-crossing (Re λ = 0 boundary → γ) or Route B performs when it locates p\*. I drove `autotune` (Python, no GPU, CWD `C:\ORRERY`) against pre-registered targets:

**Instrument trust:** `autotune --selftest` → **SELFTEST PASS** (6/6, incl. the negative control "wrong target → G-OFF-TARGET fires, exit 1"); the tool's own frozen golden is `c79002f23cf236baab5ecdb5753603a7a3853199f750b0139771d4e6cdd55bbe`.

**Run A — CSS critical-point locate** (pre-registered p\*/critical-value = 0.37, tol 0.02):
```
python tools/autotune/autotune.py --objective threshold --obj-center 0.37 --obj-width 0.05 \
       --lo 0 --hi 1 --locate crossing --level 0.5 --target 0.37 --tol 0.02 --seed 0 --json
```
→ `x_located = 0.369950`, `located_error = 5.0e-5`, `on_target = true`, `G-OFF-TARGET` **not** fired, **verdict pass, exit 0**.
**Declared blake2b-256 = `6dd1fb038bd12a69795732d7a142e317c454dad9196c93c1424e475d73344585`.**

**Run B — eigenvalue-route locate** (pre-registered γ = 0.374, tol 0.02):
```
python tools/autotune/autotune.py --objective threshold --obj-center 0.374 --obj-width 0.03 \
       --lo 0.2 --hi 0.6 --locate crossing --level 0.5 --target 0.374 --tol 0.02 --seed 0 --json
```
→ `x_located = 0.374004`, `located_error = 4e-6`, `on_target = true`, **verdict pass, exit 0**.
**Declared blake2b-256 = `0602eb2d6a02b03060826087c12239077afed315ce2f609a277ff5511f857acd`.**

**What this proves / does not prove.** It proves — evidence-grade — that a **deterministic, pre-registered level-crossing locate lands on a target to <5×10⁻⁵ and fires its off-target gate honestly** (the negative control in the selftest confirms it *would* have flagged a miss). That is the exact search primitive the deciding experiment uses to locate the eigenvalue-crossing / p\*. It does **not** prove γ=0.374 for the scalar — the objective here is `autotune`'s analytic sigmoid stand-in, not `css_echo` (which does not exist). Any γ figure remains `[ARGUMENT-GRADE]` until `css_echo` is built and its golden frozen. Artifacts: `experiments/css-selfsimilar/autotune_receipts.txt`.

---

## 6 · Honest risks + STEAL list

**Biggest risk (single):** **R2 — constructing the DSS scalar solution and its Floquet spectrum is real numerical-relativity work, not a quick script.** The persona's elegant pitch ("just an ODE eigenvalue problem, determinism is trivial") is *cleanly true only for the CSS radiation fluid*; the **massless scalar is DSS**, so the crown requires (i) solving a τ-periodic PDE BVP with regularity at both the center and the coordinate-singular past self-similarity horizon, then (ii) a monodromy/Floquet eigen-extraction. Gundlach 1997 shows it is *doable* and single-file-scale, but the schedule risk is the `construct` BVP's iteration count and horizon-regularity handling, which could blow the 5-min budget or need care I can't guarantee blind. **This is where CSS is most likely to die** — and I flag it rather than hide it (D-016/D-021 discipline).

**Other risks:** T\* estimation smearing Route B's log-periodicity (R1, mitigated by fitting T\* as a periodicity-sharpness argmax); a plausible-but-wrong λ₀ from an under-resolved monodromy (guarded by `G-CONVERGE`); and the honest possibility that the *fluid* prototype nails β≈0.36 but the *scalar* DSS extension is what fails (covered by the reinstatement trigger — the fluid oracle still ships as a lesser crown).

**STEAL list (harvest regardless of CSS's fate):**
1. **γ = 1/Re λ₀ as the crown observable** — the eigenvalue is *by far* the most golden-friendly path to γ (a single number from a linear operator, converging under refinement) versus any mass-scaling slope near a wandering p\*. **Any** winning approach should extract γ this way if it can, even a double-null or AMR one that first *finds* the critical solution then linearizes about it.
2. **`G-UNIQUE` (exactly-one-growing-mode) as a physics-honesty gate** — a cheap, decisive check that you are on the *critical* solution and not a nearby imposter; portable to every approach.
3. **`G-CONVERGE` (eigenvalue stationary under x-refinement)** — the exact anti-D-021 metamorphic gate: it *demands* the thing D-021's grid failed to do (converge), turning "does the grid see it?" into a pass/fail.
4. **Similarity coordinates `(τ,x)` as a resolution primitive** — even if CSS-the-crown loses, evolving *in* similarity coordinates (so the echo recurs at fixed x) is a resolution trick a double-null or spectral approach can adopt to stop the structure shrinking off the grid.
5. **Fit-T\*-by-periodicity-sharpness** — a clean, deterministic way to pin T\* (a 1-D argmax) that any echo-detecting approach can reuse.

---

*CSS advocacy, phase 1. The crown I offer is the eigenvalue: γ = 1/Re λ₀, a deterministic linear-algebra number on a coarse fixed grid — the most honest determinism story on the board. The cost I do not hide: the massless scalar is DSS, so that number lives behind a period-Δ Floquet construction that is genuine NR work. Build the CSS fluid prototype first (Evans–Coleman β≈0.36) as the warm-up oracle; extend to the scalar DSS crown; freeze both goldens. Refute D-021 by deleting refinement, not adding it.*
