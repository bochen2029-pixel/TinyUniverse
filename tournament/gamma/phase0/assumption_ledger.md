# ASSUMPTION LEDGER — phase 0 (the γ fork)

**Scope.** Mined from the N0 build session (commit `8060bd0`, `substrate/substrate_nexus.cpp` v1.0.0, contract `contracts/substrate_nexus.contract.md`, decision `DECISIONS.md` D-021). Part I: the forks I settled **by argument, inside one epistemic household (me), with no head-to-head experiment** — each with the argument used, the pre-registered deciding experiment that would settle it, and an assigned attack stance. Part II: assumptions the whole N0 build rests on that nothing challenged.

**Standing fact that colors every entry (borrowed from the_brain phase 0).** Every verdict below was reached by a single builder iterating against reality alone. By the tournament's own epistemics, that makes each one **argument-grade, not evidence-grade** — until its deciding experiment runs through a built tool / ORRERY. D-021 (the crown deferral) is honest, but it is *a conclusion I reached solo*; this ledger exists to let five perspectives attack it and an experiment settle it.

**Ground truth (the anchor every deciding experiment cites).** Massless real scalar, spherical, 4D: **γ ≈ 0.374** (Choptuik, PRL 70 (1993) 9; Gundlach γ = 0.3737); discrete self-similarity echo period **Δ ≈ 3.44**. These are the numbers a winning approach must recover.

---

## PART I — FORKS SETTLED BY ARGUMENT

### F1 · Coordinate formulation
- **Fork:** polar-areal (Schwarzschild-like, `ds² = −α²dt² + a²dr² + r²dΩ²`, constrained evolution) — the incumbent — vs **double-null / characteristic** (`ds² = −α²du dv + r²dΩ²`; the grid follows the collapse) vs **CSS / self-similar coordinates** (`τ = −ln(T*−T)`; discrete self-similarity becomes periodic-in-τ, resolvable on a fixed τ-grid).
- **Settled:** polar-areal, by default — it is Choptuik's original and the most-documented constrained scheme; I never evaluated the alternatives.
- **Argument (solo):** constraint integration is a clean radial ODE; log-metric integration guarantees positivity; I had it working for the transition.
- **Deciding evidence:** build each formulation's minimal EKG tool; drive `autotune` (level-crossing) to locate p\*; fit the scaling exponent; PASS = |γ_fit − 0.374| < 0.05 at R² > 0.99, deterministic golden, ≤ 5 min. Garfinkle (PRD 51 (1995) 5558) reports γ from a **uniform double-null grid** — so the resolvability claim against polar-areal is directly testable.
- **Attack stance:** **double-null advocate.** Line: the uniform grid failed *for polar-areal* because polar slicing crowds the resolution into a coordinate singularity the grid can't follow; double-null coordinates advect the grid with the imploding pulse, so the self-similar structure stays resolved at fixed cost — the D-021 corpse is a *gauge* artifact, not a grid law.

### F2 · Resolution strategy
- **Fork:** uniform grid (**graveyard, D-021**) vs **Berger–Oliger AMR** (recursive refinement where the field sharpens — Choptuik's original route to γ) vs **spectral / pseudospectral** vs **mesh adapted in the similarity coordinate** (rides on F1's CSS).
- **Settled:** uniform, then declared unable to resolve γ (D-021) and the exponent deferred.
- **Argument (solo):** AMR is heavy, hard to make deterministic/single-file/golden-freezable; a uniform grid at the N0 budget cannot see the echoing; so defer γ to the GPU stage.
- **Deciding evidence:** for the AMR path — can a *deterministic* Berger–Oliger refinement (fixed refinement criterion, fixed subcycling) be built single-file and golden-frozen, and does it land γ? For CSS — does a fixed τ-grid resolve γ with **no** adaptivity at all (the strongest determinism story)?
- **Attack stance:** **CSS/self-similar advocate.** Line: AMR is the *brute-force* answer and it fights the doctrine (nondeterminism risk, library weight); the *elegant* answer changes coordinates so the thing you're resolving stops shrinking — in similarity coordinates the discretely-self-similar solution is periodic, so a coarse fixed grid suffices and determinism is trivial. Refuting D-021 by deleting the need for refinement, not by adding it.

### F3 · The crown observable
- **Fork:** mass-scaling `M_BH ∝ (p−p*)^γ` (supercritical) vs **subcritical curvature scaling** `max_t R ∝ (p*−p)^{−2γ}` (Garfinkle–Duncan, PRD 58 (1998) 064024) vs **echo period** Δ ≈ 3.44 (a discrete-self-similarity signature, arguably easier than the exponent) vs **REFRAME the crown** entirely.
- **Settled:** I tried mass-scaling (grid-quantized, noisy) then subcritical curvature (grid-capped); both failed on the uniform grid; reported an effective exponent as a labelled diagnostic.
- **Deciding evidence:** on the winning F1/F2 substrate, measure ≥ 2 independent observables (mass-scaling AND curvature-scaling AND/OR echo Δ) and require **redundant recovery** — they must agree on γ. Echo Δ is its own PASS: detect log-periodicity with period 3.44 ± tol.
- **Attack stance:** **echo-period advocate.** Line: the *exponent* is the hardest thing to extract on any finite grid, but the *discrete self-similarity* (the echoing, period Δ≈3.44) is a qualitative log-periodic signature that shows up at far coarser resolution; N0's crown should be "we see the echo" — a receipt as deep as γ and far cheaper — with the number left to the fold where it's cheap.

### F4 · What N0's crown *should be* (the meta-fork)
- **Fork:** N0's crown must be **precise γ** (the v2 proposal §5 named it) vs N0's crown should be a **deepest-GR receipt a CPU-modest grid nails** — e.g. **boson-star mass–radius relation** or **oscillaton fundamental frequency** — with γ **deferred to the GPU N3 `horizon` stage** where resolution is cheap.
- **Settled:** I deferred γ (D-021) but kept the *transition* as N0's crown and reported γ_eff as a diagnostic — a middle path nobody argued for on purpose.
- **Deciding evidence:** does a boson-star M(R) curve / oscillaton frequency come out to <1% of the analytic value on a uniform grid at the N0 budget (making it a *clean, honest* crown), and is it "deep GR" enough to be the foundation-stone receipt? Versus: does any F1/F2 approach make precise γ cheap enough that reframing is unnecessary?
- **Attack stance:** **reframe advocate (the persona solo-me never surfaced).** Line: the proposal *assumed* γ is N0's crown; but N0's job is to be a trustworthy fp64 oracle before GPU code — and a boson-star mass–radius relation (a genuine strong-gravity, self-gravitating bound-state receipt) is exactly that, nails to <1% on a uniform grid, and doesn't require solving the single hardest problem in numerical relativity on a CPU. γ belongs to N3, where AMR-on-GPU is affordable. Chasing γ in N0 is scope error.

### F5 · Initial data & the critical search
- **Fork:** time-symmetric Gaussian (incumbent) vs ingoing pulse (Garfinkle's choice) vs a scale-parameter family; and the p\* search: hand bisection (incumbent) vs **`autotune` level-crossing against a pre-registered target** vs `mcts`.
- **Settled:** time-symmetric Gaussian + hand-rolled bisection (20 iters).
- **Deciding evidence:** replace the bisection with `autotune --locate crossing --target <p*_expected>` driving the EKG tool's collapse flag; does it locate p\* deterministically and cheaply? (This is the direct ORRERY-arming test — and autotune's proven job.)
- **Attack stance:** **search-methodology attacker.** Line: hand bisection with a floating collapse criterion is exactly why p\* wandered; a pre-registered `autotune` crossing with a fixed metric is deterministic and cites its own hash — the instability was partly the *search*, not only the *grid*.

---

## PART II — UNCHALLENGED ASSUMPTIONS

- **A1 · The formulation choice is free of the observable choice.** N0 assumed you pick a formulation then pick an observable; in reality CSS coordinates *make* the echo/period the natural observable, and double-null *makes* the mass-scaling clean — F1 and F3 may not be independent. A refuter should test whether the winning formulation forces the observable.
- **A2 · A uniform grid is the only "deterministic, single-file" option.** D-021 implicitly equated "no AMR" with "uniform." CSS coordinates and characteristic grids are *also* non-adaptive and deterministic — the assumption that resolving γ requires abandoning determinism is the load-bearing, untested premise of the whole deferral.
- **A3 · γ ≈ 0.374 is the only anchor.** The oscillaton frequency, boson-star M(R), and echo Δ≈3.44 are *also* anchors of the same physics; an approach that nails one of them is not worthless just because it doesn't nail γ. The crown-reframe fork (F4) lives here.
- **A4 · The N0 CPU budget (~5 min) is fixed.** If a CSS or double-null approach resolves γ in 30 s, the budget was never the constraint; if only AMR does and it needs 20 min, the budget vs the crown is a real trade the operator should price, not the builder.

---

## Pre-registered deciding experiments (the ORRERY-armed settlement)

| DE | settles | mechanism | PASS bar |
|---|---|---|---|
| **DE-γ** | F1, F2 | build formulation-X EKG tool → `autotune` locates p\* → fit mass-scaling exponent | \|γ_fit − 0.374\| < 0.05, R² > 0.99, deterministic golden, ≤ 5 min |
| **DE-redundant** | F3 | on the winning substrate, measure mass-scaling AND curvature-scaling (and/or echo Δ) | the two observables agree on γ to < 0.03; echo Δ = 3.44 ± 0.3 |
| **DE-reframe** | F4 | build boson-star M(R) / oscillaton-frequency tool on a uniform grid | match analytic to < 1% at the N0 budget (⇒ a clean crown exists without γ) |
| **DE-search** | F5 | `autotune --locate crossing` drives the collapse flag vs a pre-registered p\* | locates p\* deterministically; cite the declared blake2b |

*The graveyard already holds one corpse (uniform-ultrafine polar-areal, D-021). Phase 1 opens the five stances against these forks; the operator reviews this ledger before it does.*
