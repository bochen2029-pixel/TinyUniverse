# APPROACH — REFRAME THE CROWN: N0's crown is a boson-star bound-state receipt, not γ

**Persona:** reframe-the-crown advocate (F4; the stance solo-me never surfaced). **Fork attacked:** F4 (meta) primarily; F3 (observable) and A3 (γ is not the only anchor) consequentially. **Phase:** 1 (approach doc; argument-grade until DE-reframe runs). **ORRERY-armed:** posit (parsimony), autotune (crossing-locate) — both run below, hashes cited.

---

## 0 · Thesis (two sentences)

N0's actual job is to be a **trustworthy fp64 gravitational oracle before any GPU code exists** — not to solve the single hardest problem in numerical relativity (a self-similar, discretely-self-similar critical attractor requiring adaptive mesh refinement) on a uniform CPU grid. A **boson-star mass–frequency relation M(ω)** — the self-gravitating ground state of a massive complex scalar, obtained by shooting the Einstein–Klein–Gordon eigenvalue problem on the *same* uniform grid N0 already runs — is a genuine strong-gravity, self-gravitating **bound-state** receipt that nails a known analytic/NR constant (the maximal mass **M_max ≈ 0.633 M_Pl²/m**) to **< 1 %** at the N0 budget, satisfying the CHARTER's oracle-honesty lens far better than a heroic, grid-noise-limited γ; precise γ belongs to the GPU **N3 `horizon`** stage where AMR is affordable.

---

## 1 · The reframed crown, fully specified

### 1.1 The object

A **mini boson star** (Kaup 1968; Ruffini & Bonazzola 1969): a self-gravitating equilibrium of a *massive complex* scalar field with **no self-interaction**, held up against its own gravity by the dispersive/quantum pressure of the Klein–Gordon field. It is the scalar-field analogue of a neutron star and a genuine strong-gravity object: the compactness `2M/R` reaches ~0.1–0.3 near the maximal-mass configuration — deep in the nonlinear GR regime, not a weak-field toy.

Ansatz (harmonic time dependence makes the stress-energy **static**, so the metric is time-independent):

```
  φ(t,r) = σ(r) · e^{ -i ω t }            (complex scalar, frequency ω, real profile σ(r))
```

Metric (static, polar-areal — the SAME gauge substrate_nexus.cpp already integrates):

```
  ds² = − α(r)² dt² + a(r)² dr² + r² dΩ²
```

### 1.2 The equations (the EKG eigenvalue/shooting system)

With potential `V = μ²|φ|²` (mass μ; set μ = 1 to fix units, so lengths are in 1/μ and masses in M_Pl²/μ), the coupled first-order ODE system for the ground state is:

```
  (ln a)_{,r} = (1 − a²)/(2r) + 4π r [ (ω²/α² + μ²) a² σ² + σ'² ]           (Hamiltonian constraint)
  (ln α)_{,r} = (a² − 1)/(2r) + 4π r [ (ω²/α² − μ²) a² σ² + σ'² ]           (slicing)
  σ''          = −( 2/r + (ln α)_{,r} − (ln a)_{,r} ) σ' + a² ( μ² − ω²/α² ) σ   (KG, static)
```

Boundary conditions:

```
  r = 0 :  σ(0) = σ_c (the shooting family parameter),  σ'(0) = 0,  a(0) = 1,  α(0) = α_c
  r → ∞ :  σ → 0  (bound state), and  σ has NO interior nodes  (the GROUND state)
```

The two constraint equations replace the incumbent's exact machinery (`Fa`, `Fal`, `solveMetric` in `substrate_nexus.cpp`, lines 141–173) — the `4π r a² V` and `2π r(Π²+Φ²)` source terms become the static `4π r(ω²/α² ± μ²)a²σ² + σ'²` terms. **This is not a new solver; it is the incumbent's constraint integrator with a static scalar source.**

### 1.3 The shooting/eigenvalue method (uniform grid, deterministic)

For a chosen central amplitude `σ_c`, `ω` is an **eigenvalue**: only a discrete set of ω admits a solution that decays (rather than diverging to ±∞) at large r. The ground state is the *lowest* such ω (nodeless σ).

```
  1.  Fix σ_c.  Guess ω.  Set a(0)=1, α(0)=1 (rescale α at the end so α·a → 1 at r_max).
  2.  Integrate the 3 ODEs OUTWARD from r=0 to r=r_max with RK4 on the uniform grid
      (dr fixed — the incumbent's grid; no adaptivity, no refinement).
  3.  Read the residual: for the eigen-ω, σ(r_max) → 0 from above without a node; for ω too
      low σ diverges to +∞, for ω too high σ crosses zero (a node) and diverges to −∞.
      Define R(ω) = sign of the large-r divergence (a monotone step across the eigenvalue).
  4.  BISECT / level-cross ω on R(ω) until |R| < tol  →  the ground-state ω(σ_c).
  5.  The ADM mass is M = m(r_max) = (r_max/2)(1 − 1/a(r_max)²);  R_99 = radius enclosing
      99% of M.  Record the point (ω, M, R_99).
  6.  Sweep σ_c over a range  →  the curve M(ω) (and M(R_99)).
```

Step 4 is **exactly** the incumbent's `bisect_pstar` pattern (lines 333–343) with the collapse flag replaced by the divergence sign — and exactly what ORRERY's `autotune --locate crossing` does mechanically (demonstrated in §4.3).

### 1.4 The observable + its analytic/NR anchor

The curve `M(ω)` (equivalently `M(σ_c)`) rises from the low-amplitude Newtonian limit, reaches a **maximum**, then turns over — the maximum marking the boundary between stable (dω branch) and unstable configurations. The anchor is that maximum:

> **M_max ≈ 0.633 M_Pl² / μ** for the ground-state, non-self-interacting mini boson star, at central amplitude σ_c ≈ 0.27 (in the μ = 1, √(4π)-absorbed convention) and frequency ω/μ ≈ 0.85.

This constant is the boson-star field's most-cited number (Kaup 1968; Ruffini & Bonazzola 1969; tabulated in Liebling & Palenzuela's Living Review). Its precise value in a given convention depends on the normalization of the scalar-gravity coupling; **the *robust, convention-independent* receipt is the existence and location of the maximum and M_max to the stated ~0.63 M_Pl²/μ, matched to < 1 %.** (The exact numeric — 0.633 vs 0.6330 — is fixed by the built tool against a stated convention; do not over-precision it here. `[VALUE: 0.633 — standard, but reproduce in-tool against the declared convention]`.)

**PASS bar (pre-registered):** the shot ground-state M(ω) curve exhibits a single maximum, and **M_max matches the analytic constant to < 1 %** on the uniform grid at the N0 budget, deterministic, ≤ 5 min. A second, independent receipt for redundancy: the maximal-mass config also satisfies **dM/dω = 0 AND dN/dω = 0 coincident** (mass and particle-number extrema coincide — the standard stability turning-point theorem), giving a metamorphic cross-check that does not depend on the mass normalization.

---

## 2 · Why this is a BETTER N0 crown than γ (the oracle-honesty lens)

The CHARTER's **oracle-honesty** lens (constraint 2; D-016/D-021 discipline) asks: is the reported number **anchored**, **metamorphic-stable**, and **redundantly recovered** — or is it grid-noise dressed as a result?

| criterion | precise γ on a uniform CPU grid | boson-star M(ω) |
|---|---|---|
| **anchored** | γ = 0.374 exists, but the *uniform-grid measurement* yields ≈ 0.24 (D-021, measured) — the anchor is missed by 36 % | M_max ≈ 0.633 M_Pl²/μ; equilibria are smooth (no self-similar blow-up), so a uniform grid resolves them to < 1 % |
| **metamorphic-stable** | **NO** — the measured p\* itself wandered 0.40 → 0.356 as N: 800 → 1600 (D-021); refining makes it *worse* (chaotic) | **YES** — refine dr and M_max converges monotonically (2nd-order); the curve is a fixed function, not a moving target |
| **redundantly recovered** | the two observables (mass-scaling, curvature-scaling) both saturate at the grid ceiling ≈ (φ/dr)² — they agree on *noise* | **YES** — dM/dω = 0 and dN/dω = 0 coincide at M_max independently of normalization; two receipts, one config |
| **honest verdict** | must ship as a *labelled diagnostic* `gamma_eff_gridlimited`, never gated (D-021) — i.e. N0 currently has **no gated crown number**, only a gated *phenomenon* | a **gated crown number** (M_max < 1 %) a stranger runs cold from BUILD.md in one pass |

The decisive point: **γ on a uniform grid cannot pass its own oracle-honesty test — the incumbent already conceded this (D-021) and ships it ungated.** A foundation-stone oracle whose crown number is un-gatable is a weak foundation stone. The reframe gives N0 a crown number that *is* gated, *is* metamorphic-stable, and *is* redundantly recovered — precisely the D-016 discipline (the 7π withdrawal: ship what you can gate, defer what you cannot).

### 2.1 Why γ → N3 is correct scoping (not a dodge of the hard problem)

The v2 proposal's **own §6** prices this verbatim: *"fixed-grid horizon resolution bounds the smallest honest BH; declare it; **AMR is a later contract**."* The proposal §5 milestone table further assigns **N3 `horizon`** the gate *"N0 Choptuik reproduced on GPU at declared tolerance"* — i.e. the architecture *already intends γ to be re-attempted where refinement is affordable.* Chasing precise γ *on the CPU in N0* therefore duplicates, at the worst possible resolution budget, a receipt the ladder already schedules for the stage built to afford it. The scoping is not my invention; it is the proposal's, and D-021 is the measured confirmation that N0 is the wrong stage for it.

A boson star, by contrast, is *on the N1–N4 critical path*: the proposal's "dream sentence" (§S3/§S4) is **"a ball of field-hydrogen contracts → pressure supports it → fuel depletes → collapse resumes."** The pressure-supported equilibrium the star **relaxes toward is a boson star** (for the massive-scalar substrate) / an oscillaton. N0 shipping the boson-star ground state is N0 shipping the *first receipt N4 will cross-check against* — a foundation stone that bears load downstream, whereas a critical exponent (an *unstable* attractor no star sits at) is a beautiful cul-de-sac at N0.

---

## 3 · Guarantee-audit against the fixed constraints (this is my strength — and one honest weakness)

**Constraint 1 — Determinism-or-it-doesn't-ship.** A static ODE shoot on a fixed grid has **no time evolution, no chaotic near-critical regime, no floating collapse criterion** — the three sources of the incumbent's p\* wander. `(σ_c, ω-bracket, grid) → byte-identical (ω, M, R_99)`. Bisection to a fixed residual tol is deterministic. **Passes comfortably** — arguably a *stronger* determinism story than the incumbent's evolution battery.

**Constraint 2 — Golden-gated, honest.** Anchor = M_max (analytic). Metamorphic = convergence under dr-refinement (must *improve*, unlike γ) + node-count invariance (ground state stays nodeless as σ_c varies below max). Redundant recovery = the coincident dM/dω = 0, dN/dω = 0 turning point. **All three legs present** — the D-016 triad, satisfied. Golden freezes `(σ_c-grid, dr, r_max) → blake2b` of the M(ω) table.

**Constraint 3 — CPU ≤ ~5 min, single file.** One shoot = one outward RK4 integration to r_max (N ~ 10³ points). One eigenvalue = ~20–30 shoots (bisection). One M(ω) curve = ~20–40 eigenvalues = ~600–1200 integrations of a 3-ODE system on 10³ points — **seconds, not minutes** (the incumbent's *evolution* battery, far heavier, runs in ~20 s). Single-file C++17, fp64, stdlib-only — same shape as `substrate_nexus.cpp`. **Passes with orders of magnitude of headroom.**

**Constraint 4 — Exit codes.** 0 = M(ω) shows a maximum and M_max within 1 %; 1 = declared gate fired (no maximum found, or M_max off by > 1 % — a real negative result); 2 = error (integration overflow, no eigenvalue bracketed). Clean.

**The honest weakness (stated plainly, per the persona brief): does it duck the challenge?**
The tournament question is literally *"how should N0 resolve γ."* Reframing answers *"it shouldn't — here is a different, cleaner crown."* That is a legitimate move **only if** the reframe is *at least as deep* as γ and the deferral is *principled, not evasive*. Two honest concessions:
1. **A boson star is a static eigenvalue problem; γ is a dynamical, self-similar, symmetry-breaking phenomenon.** γ is, by a defensible measure, *deeper physics* — a universal critical exponent versus a bound-state mass. The reframe trades **depth-of-phenomenon** for **honesty-of-measurement**. A judge who weights "hardest correct physics" over "cleanest honest receipt" should score γ-if-resolvable above this.
2. **It leaves the tournament's headline question partially unanswered.** If a *sibling* approach (double-null, CSS) can resolve γ deterministically in budget (their claim), then the reframe was unnecessary and the honest verdict is "resolve γ, don't reframe." **My claim is conditional:** reframe *iff* no buildable N0-budget instance lands γ within tol. The reframe is the correct floor, not necessarily the ceiling. I hold this doubt in the register (§4.1) and do not pretend the deferral is free.

This weakness is real and I will not paper over it. But note its shape: it is a weakness of *ambition*, not of *correctness or feasibility* — the opposite failure mode from the incumbent's γ, whose weakness is that it cannot be measured honestly at all. Between "less ambitious but certainly right" and "maximally ambitious but currently un-gatable," a *foundation stone* should be the former.

---

## 4 · Pre-registered falsifier + costed deciding experiment

### 4.1 Reversibility lemma (CHARTER house discipline — the reinstatement trigger)

**If a sibling approach demonstrates, via a built single-file CPU tool at the N0 budget, `|γ_fit − 0.374| < 0.05` at `R² > 0.99` with a deterministic golden, THEN the reframe is demoted from "N0's crown" to "an N0 companion receipt," and γ is reinstated as the crown.** The reframe deletes *nothing* (it adds a boson-star tool); it only *reprioritizes* the crown. Reinstating γ costs no rework — the boson-star tool ships regardless as an N4 cross-check. This is the cheapest possible reversibility: the graveyard entry (if any) is "γ-as-N0-crown," and its reinstatement trigger is a sibling's green γ golden.

### 4.2 The pre-registered falsifier (what kills THIS approach)

The reframe is **KILLED** if any of:
- **(F-anchor)** the built shooting tool cannot reproduce M_max to < 1 % on the uniform grid at ≤ 5 min — i.e. boson stars are *also* grid-limited at the N0 budget (falsifies the "clean crown exists" claim). *Prediction: it will pass — equilibria are smooth; this is the falsifiable bet.*
- **(F-metamorphic)** M_max does **not** converge monotonically under dr-refinement (would mean the receipt is as grid-noisy as γ — falsifies the oracle-honesty advantage).
- **(F-depth)** the operator/judges rule that a bound-state mass is *not deep enough* to be N0's foundation-stone crown regardless of its honesty (falsifies the *sufficiency* of the reframe — the §3 weakness realized).

### 4.3 The costed deciding experiment — DE-reframe

**Tool contract (does not exist yet — specified here; label any figure `[ARGUMENT-GRADE]`):**

```
  bosonstar_nexus  ·  single-file C++17, fp64, stdlib-only, NO GPU  ·  units μ = 1 (G = c = 1)
  ── the incumbent's blake2b + declared-JSON/golden idiom, lifted verbatim (D-020 shape) ──

  face:  bosonstar_nexus [--json] [--selftest] [--golden] [--sigmac S] [--N n] [--rmax R]
  method: outward RK4 of the 3-ODE static-EKG system; ω bisected on the large-r divergence sign;
          ADM mass from the mass aspect at r_max; sweep σ_c → M(ω) table.

  battery:
    B1  weak-field limit  — σ_c → 0 gives ω/μ → 1⁻ and M → 0 (Newtonian boson-star limit)   [oracle: KG dispersion]
    B2  ground-state node — the located eigen-σ is nodeless; the next eigenvalue has exactly 1 node   [oracle: Sturm–Liouville]
    B3  MASS CURVE (crown) — M(ω) has a single maximum; M_max within 1 % of ≈0.633 M_Pl²/μ    [oracle: Kaup / Ruffini–Bonazzola]
    B4  turning-point     — dM/dω = 0 and dN/dω = 0 coincide at M_max (redundant recovery)      [oracle: Friedberg–Lee–Pang stability theorem]  [CITE-UNVERIFIED: FLP exact ref]
    B5  determinism       — battery re-run in-process → byte-identical declared JSON

  PASS bar (crown, pre-registered):  B3 |M_max − M_max^analytic| / M_max^analytic < 0.01,
                                     with B1,B2,B4,B5 green, deterministic golden, wall-clock ≤ 5 min.
  exit: 0 all green · 1 a declared gate fired (e.g. M_max off > 1 %, or no maximum) · 2 error.
```

**Cost:** ~1 engineer-day to adapt `substrate_nexus.cpp`'s constraint integrator + bisection into the static shoot (the solver, bisection, blake2b, JSON, golden idiom all already exist and are lifted). Runtime: **seconds**. This is the *cheapest* deciding experiment of any approach in the tournament — because it reuses the incumbent almost entirely and adds no time-evolution.

---

## 5 · ORRERY runs (evidence-grade, not argument-grade)

Read-only ORRERY, CWD `C:\ORRERY`. Both envelopes re-run 2× byte-identical. Full ledger: `experiments/reframe-crown/RESULTS.md`.

### 5.1 posit — the Occam argument, made numeric

**Baseline golden:** `python tools/posit/posit.py --golden --json` → `GOLDEN OK`, blake2b `7a22dd229a42ce46a6c102f0545f83022b975dc39d5f1794cd6019e6f5a20e44`, exit 0.

**The audit** (`posit_case_crown.json`, pinned `case_blake2b = bfd12420…`): PATCHWORK = *"chase precise γ on the CPU"*; UNIFIED = *"reframe to a boson-star receipt."* Targets: {strong-gravity self-gravitating bound state, analytic anchor matched < 1 %, deterministic single-file CPU ≤ 5 min, generalizes to N1–N4}. Result:

```
  patchwork physics budget = 5.0   (2 posits + 3 UNBACKED bridges: AMR-is-deterministic,
                                     AMR-resolves-γ-in-budget, γ-generalizes)
  unified   physics budget = 1.4   (1 new posit [the stationary complex-scalar ansatz]
                                     + 2 imports [already-owned constraint solver, shooting
                                     root-find] + 4 free derivations)
  delta_physics = +3.6   same_reach = true   floating = []   parsimony = "WIN"   exit 0
```

Declared envelope blake2b: `13e2248608a14e089002969d7fd0bb2d5cab429b15a7e7ff712b55a096529e99`.

**The Occam claim, now evidence-grade:** at **equal reach** and **no floating**, adding a boson-star receipt spends a physics-layer posit budget **3.6 units below** chasing γ on the CPU. The reframe's *entire* new epistemic cost is **one brute posit** — the stationary ansatz `φ = σ(r)e^{−iωt}` (Kaup/RB); everything else derives from the constraint solver and shooting root-find the incumbent *already owns*. Chasing γ must instead posit a **new** AMR/CSS mechanism (the shipped uniform grid provably can't, per D-021) and carry **three unbacked bridges**. *(Honesty, per the ledger: a stricter equal-*items* framing returns `parsimony=reject`/exit 1 because chasing-γ cannot even COVER the "generalizes to N1–N4" target — an even sharper result. The cited win-variant restores equal reach by pricing γ's generalization as an explicit bridge, so +3.6 is a licensed parsimony number, not a reach artifact. Reject-variant retained: `posit_case_crown.reject_variant.json`, exit 1, blake2b in RESULTS.md.)*

### 5.2 autotune — the ground-state finder is a level-crossing

**Baseline golden:** `python tools/autotune/autotune.py --golden --json` → `GOLDEN OK`, blake2b `c79002f23cf236baab5ecdb5753603a7a3853199f750b0139771d4e6cdd55bbe`, exit 0.

**The demonstration** (`--objective threshold --obj-center 0.8267 --locate crossing --level 0.5 --target 0.8267 --tol 0.01`): a level-crossing locate recovers a **pre-registered eigenfrequency target ω/μ = 0.8267** to `located_error = 2e-6`, `on_target = true`, `G-OFF-TARGET` not fired, exit 0. Declared envelope blake2b: `11b87fb84cb8f7d9bb884073a293febf00f835c1b09ac766f7e7c92f7b1ea9ea`.

`[ARGUMENT-GRADE]`: the built-in logistic is a **stand-in** for the real EKG shooting residual — the point is that autotune's `--locate crossing` against a pre-registered target **is** the boson-star ground-state finder (tune ω until the nodeless-decay residual crosses zero), mechanically and with a gate. autotune already located ratchet's ρ_c = 0.2581 blind against a *real* tool (MODULE.md), so the primitive is proven; here it hits a pre-registered eigenvalue to 2e-6. **0.8267 is a plausible placeholder ω/μ, NOT a measured boson-star eigenvalue** — the real number is DE-reframe's B3 output. No fabricated physics numbers.

---

## 6 · Honest risks + STEAL list

### Risks (owned, not hidden)

1. **Depth-sufficiency (the big one).** A judge weighting "hardest correct physics" may rule a bound-state mass is not a worthy *foundation-stone* crown versus a universal critical exponent. **This is the reframe's real exposure — it may be read as ducking the tournament's headline question.** Mitigation: the reversibility lemma (§4.1) makes reframe the *floor*, reinstatable the instant a sibling lands γ in budget; and a boson star is on the N1–N4 critical path where γ is not — a foundation stone that bears downstream load.
2. **The maximal-mass constant is convention-sensitive.** 0.633 M_Pl²/μ assumes a specific normalization of the scalar–gravity coupling; a sloppy in-tool convention could "miss" a correct answer. Mitigation: B3 states its convention and cross-checks with the convention-*independent* turning-point receipt (B4).
3. **`[CITE-UNVERIFIED]` on the FLP stability turning-point exact reference** (Friedberg–Lee–Pang, ~1987) — the *theorem* (dM/dω and dN/dω coincide at the stability boundary) is standard and in Liebling & Palenzuela's review, but I have not pinned the primary citation; verify before the golden.
4. **Ground-state node-selection robustness.** The bisection must reliably separate the nodeless branch from the 1-node branch near the maximum; a coarse grid could misidentify. Mitigation: B2 gates node count explicitly; refine dr if B2 is marginal (cheap — it's a static shoot).

### STEAL list (harvest regardless of the reframe's fate)

- **From the echo-period advocate (F3):** if γ is reinstated, the *echo Δ ≈ 3.44* is a far cheaper log-periodic receipt than the exponent — the reframe and echo compose (ship the honest boson-star crown *now* AND the echo when N3 resolves it).
- **From the search-methodology attacker (F5):** the incumbent's p\* wander was *partly the search* (floating collapse criterion). The reframe's `autotune --locate crossing` with a **pre-registered target + G-OFF-TARGET** is the disciplined search primitive — steal it for *whichever* observable wins (γ or M(ω)), so the located value cites its own hash.
- **From the CSS advocate (F2/A2):** the insight that "change coordinates so the thing you resolve stops shrinking" — the boson-star ansatz does exactly this in the *time* direction (harmonic dependence makes the source static, killing the evolution entirely). Same spirit, applied to N0's crown.
- **For N4 regardless of this tournament's verdict:** `bosonstar_nexus` is the ground-state oracle the "hydrogen ball → star" sentence relaxes toward — build it whether or not it wins N0's crown; it is a load-bearing N4 receipt either way.

---

## 7 · Cited literature (precise; unsure ones flagged)

- **Kaup, D. J.** "Klein–Gordon Geon." *Phys. Rev.* **172**, 1331–1342 (1968). — first self-gravitating massive-scalar equilibria (the original mini-boson-star; the maximal-mass scale).
- **Ruffini, R. & Bonazzola, S.** "Systems of Self-Gravitating Particles in General Relativity and the Concept of an Equation of State." *Phys. Rev.* **187**, 1767–1783 (1969). — the ground-state boson star as a second-quantized scalar; the canonical M(ω) / maximal-mass ≈ 0.633 M_Pl²/m result.
- **Seidel, E. & Suen, W.-M.** "Oscillating soliton stars." *Phys. Rev. Lett.* **66**, 1659–1662 (1991). — oscillatons: the *real* massive-scalar analogue (fundamental frequency as the alternative receipt).
- **Seidel, E. & Suen, W.-M.** "Formation of solitonic stars by gravitational cooling." *Phys. Rev. Lett.* **72**, 2516–2519 (1994). — oscillaton formation/stability (dynamical robustness of the receipt).
- **Liebling, S. L. & Palenzuela, C.** "Dynamical Boson Stars." *Living Reviews in Relativity* (2012; updated eds.). — the modern review; tabulated M_max, the M(ω) curve, the turning-point stability criterion, conventions. `[CITE-UNVERIFIED: exact volume/article number & latest revision year — verify before golden]`.
- **Friedberg, R., Lee, T. D. & Pang, Y.** "Mini-soliton stars" / scalar-field stability turning-point theorem, *Phys. Rev. D* ~**35** (1987). `[CITE-UNVERIFIED: exact volume/page — the dM/dω=dN/dω=0 coincidence is standard (in L&P) but the primary ref is unpinned]`.
- **Choptuik, M. W.** "Universality and scaling in gravitational collapse of a massless scalar field." *Phys. Rev. Lett.* **70**, 9–12 (1993). — the γ ≈ 0.374 crown being *deferred to N3*, cited for the scoping argument (not the reframe's receipt).

*Numbers used: M_max ≈ 0.633 M_Pl²/μ (Kaup/RB — reproduce in-tool against the stated convention), ω/μ ≈ 0.85 at max (order-of-magnitude, from the review), γ = 0.374 (Choptuik/Gundlach, incumbent's anchor). No numbers fabricated; the autotune 0.8267 is an explicit placeholder, not a physics claim.*

---

*Reframe the crown: give N0 a bound-state receipt it can nail honestly (M_max < 1 %, gated, metamorphic-stable, redundantly recovered), and let the unstable critical attractor go to N3 where AMR is affordable — the proposal's own §6. The one honest cost: it is less ambitious than resolving γ, and it holds that doubt in the reversibility lemma rather than hiding it.*
