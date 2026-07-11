# PHASE 1 — the SPECTRAL / PSEUDOSPECTRAL approach to γ ≈ 0.374

**Persona:** spectral / pseudospectral advocate (attacks F2 with the spectral fork; leans on F1's CSS option).
**Thesis (2 sentences).** The uniform grid failed (D-021) not because γ is unreachable on a CPU but because a
*fixed algebraic-order* discretization spends its degrees of freedom uniformly while the physics concentrates all
its information in a shrinking central region — a Chebyshev pseudospectral scheme spends DOF where the spectrum
lives and converges *geometrically*, buying the near-critical resolution for a few dozen–hundred points instead of
the ~10⁴–10⁶ a finite-difference grid needs. The one hard part — the self-similar feature does keep shrinking, so a
*static* spectral grid eventually Gibbs-rings just as a static FD grid saturates — is dissolved not by adaptive mesh
refinement but by a **deterministic similarity/compactified coordinate map** in which the discretely-self-similar
solution becomes (log-)periodic and therefore fixed-resolution-resolvable, which is the strongest determinism story
of any fork in this tournament.

**Standing honesty banner.** Every γ figure and every search figure in this doc is `[ARGUMENT-GRADE]`. The real
spectral Einstein–Klein–Gordon instrument does not exist yet; §4 specifies its contract. What *is* evidence-grade
here: (a) the convergence-law demonstration (`experiments/spectral/cheb_convergence.py`, declared blake2b
`b6ad2eba…`), which shows exponential-vs-O(N⁻²) on a proxy field; (b) three ORRERY `autotune` runs that show a
deterministic critical-point locate against a pre-registered target, each citing its declared blake2b. Neither
touches the Einstein equations — they make the *method claim* and the *search claim* evidence-grade, not γ.

---

## 1 · The scheme, fully specified

### 1.1 Formulation (couples F1 and F2 — and A1 predicted it would)

Spectral method and coordinate choice are not independent (ledger A1): a global Chebyshev basis wants a smooth,
bounded solution on a fixed interval, and the *only* way the near-critical EKG solution is smooth and bounded on a
fixed domain is in **similarity coordinates**, where discrete self-similarity (DSS) turns the blow-up into a
log-periodic bounded orbit. So the spectral fork *forces* the CSS formulation. I adopt it and own it.

- **Spatial coordinate.** Compactified radial similarity variable. Start from the areal radius `r` and the central
  proper time to accumulation `T*`. Define log-time `τ = −ln(T* − t)` and a similarity radius `ξ = r / (T* − t)`
  (Gundlach's `x`). Under exact self-similarity the fields `Z(τ, ξ)` are `τ`-independent (CSS) or `Δ`-periodic in
  `τ` (DSS, the scalar-field case). Map the semi-infinite `ξ ∈ [0, ∞)` to the Chebyshev interval by a rational
  (algebraic) map `ξ = L (1+y)/(1−y)`, `y ∈ [−1, 1]` (Boyd 2001, ch.17, rational Chebyshev `TB` functions), or a
  finite truncation `ξ ∈ [0, ξ_max]` with `y = 2ξ/ξ_max − 1` if the outer boundary is causally disconnected from
  the self-similar core over the evolution window. The **map parameter `L` (or `ξ_max`) is a declared, golden-pinned
  scalar** — it is the map's only freedom, chosen once, frozen; there is no runtime adaptation, so determinism is
  trivial (see §3).
- **Basis.** Chebyshev–Gauss–Lobatto (CGL) collocation, nodes `y_j = cos(πj/N)`, `j = 0..N`. Spatial derivatives via
  the dense first-derivative differentiation matrix `D` (Trefethen construction with the negative-sum-trick diagonal;
  Boyd 2001 §5). Second derivatives as `D²` (or a separately-formed `D2` for better conditioning near the ends).
- **Field variables.** The same first-order EKG variables the incumbent uses — `φ`, `Φ ≡ ∂_ξ φ`, `Π` — but expressed
  in `(τ, ξ)`; the areal metric functions `a, α` remain **constraint-solved by radial integration**, which in the
  spectral setting is a spectral integration (or a small triangular solve against the constraint ODE) rather than a
  time-marched field. Spherical symmetry keeps the constraints as radial ODEs — the incumbent's one genuine elegance,
  preserved verbatim.

### 1.2 Time integration (in log-time τ)

- Method of lines: `∂_τ Z = RHS(Z)` with the CSS/DSS source terms (the `τ`-derivative picks up the self-similar
  rescaling terms `+ (something)·ξ ∂_ξ Z + weight·Z` relative to the areal form — the standard CSS reduction).
- Integrator: **explicit RK4** as the default (matches the incumbent, keeps the golden idiom identical, no implicit
  solve → no nondeterministic linear-algebra dependency). Time step `dτ` fixed and declared.
- **Conditioning caveat, stated up front (this is a real flank).** The CGL clustering makes `‖D‖ ~ O(N²)` and
  `‖D²‖ ~ O(N⁴)`; the explicit-RK stability limit is therefore `dτ ≲ C·N⁻²` for first-order and `≲ C·N⁻⁴` if a
  second-derivative operator sets the stiffest mode. For a *hyperbolic* first-order EKG system the binding operator is
  `D` (`N⁻²`), which is affordable at the modest N spectral needs (N ~ 40–120 per domain). If a formulation exposes
  `D²` stiffness, the fallback is (i) an operator split with an implicit treatment of the stiff linear part
  (`(I − dτ·A)⁻¹` is a *fixed, deterministic* dense solve — LU with partial pivoting is bit-reproducible
  single-threaded), or (ii) domain decomposition to cap N per patch. Both are deterministic; both are costed in §3.

### 1.3 Boundary treatment

- **Origin `ξ = 0`.** Parity/regularity exactly as the incumbent: `φ, Π` even, `Φ` odd. On a CGL grid regularity is
  imposed either by a parity basis (even Chebyshev `T_{2k}` for even fields) or by a boundary-bordering condition
  `Φ(0) = 0`, `∂_ξ φ(0) = 0` enforced in the collocation system. The parity-basis route is cleaner and halves the
  DOF; I take it and declare it.
- **Outer boundary.** With the rational compactified map the outer point *is* `ξ → ∞`; the boundary condition is the
  self-similar outgoing/decay condition (fields → their asymptotic self-similar tail). With a finite `ξ_max`, an
  outgoing (Bondi-type / Sommerfeld-in-`τ`) condition, penalty-enforced (a SAT/penalty term — deterministic, standard
  in spectral hyperbolic solvers).
- **Filtering.** An exponential spectral filter (Boyd 2001, ch.9) applied each step to the highest ~⅓ of modes —
  the spectral analogue of Kreiss–Oliger dissipation the incumbent already uses. **The filter is a fixed, declared
  linear operator** (a diagonal multiply in coefficient space) → it does not break determinism and is part of the
  frozen config.

### 1.4 The map that keeps the shrinking feature resolved — and stays golden-freezable

The load-bearing move. Two layers, in order of preference:

1. **Similarity coordinates alone (no adaptivity).** If `T*` is known well enough, DSS makes the solution
   `Δ`-periodic in `τ`; a *fixed* CGL grid in `ξ` resolves it for all `τ` because the feature no longer shrinks in
   `ξ` — it *echoes in place*. This is the dream: zero adaptivity, a static frozen grid, determinism for free. It is
   exactly Garfinkle's observation that null/similarity coordinates "beautifully avoid the need for mesh refinement"
   with "perfect numerical tuning in double-precision arithmetic in spherical symmetry" (Garfinkle, gr-qc/9412008;
   see §5).
2. **If `T*` drifts (the honest complication).** `T*` is not known a priori — it is *part of what the critical search
   finds*. A mis-set `T*` makes `ξ` slowly mis-scale and the feature drifts across the grid. The fix that keeps
   determinism is a **fixed feedback law**, not an AMR library: rescale `ξ` (equivalently, shift the map's `L`) by a
   *closed-form function of a measured invariant* — e.g. hold the location of the maximum of the compactified Ricci
   scalar at a fixed collocation index by an analytic update `L_{n+1} = L_n · g(ξ_peak^n)`. This is one scalar,
   updated by a pinned formula, from a declared observable — **byte-reproducible and golden-freezable**, unlike a
   Berger–Oliger tree whose refinement *topology* is data-dependent. This is the spectral fork's answer to F2's
   determinism objection: adaptivity, if needed at all, is a scalar map parameter under a frozen control law, not a
   mesh hierarchy. **I flag this as my primary strain point** (§3, §6): layer-1 is clean; layer-2 reintroduces an
   adaptivity story I must keep to a single deterministic scalar or the DETERMINISM lens kills it.

---

## 2 · Why exponential convergence beats the uniform grid (the D-021 rebuttal — evidence-grade for the *law*)

D-021's corpse reason: a fixed grid caps the resolvable central curvature at ≈(φ/dr)², and refining does not
converge — it goes chaotic (p\* wandered 0.40→0.356, N=800→1600). The spectral claim is that the *ceiling* is an
artifact of algebraic-order accuracy, not of "fixed grids" as a class.

**Demonstrated, not asserted.** `experiments/spectral/cheb_convergence.py`
(declared blake2b `b6ad2eba341d427ffa8c825058cb6a383c4500660fabbf3f576b9920d63dcaaf`, deterministic, stdlib, runs in
<1 s) differentiates a *sharp but smooth* Gaussian pulse — the analytic caricature of a near-critical, steep-but-not-
yet-singular field — with (a) Chebyshev collocation and (b) the incumbent's centered 2nd-order stencil on a uniform
grid of the same node count:

| N (nodes) | Chebyshev max-err | centered-FD2 max-err |
|---:|---:|---:|
| 24 | 2.06e-1 | 9.91e-1 |
| 32 | 1.73e-2 | 5.37e-1 |
| 48 | 1.59e-5 | 2.68e-1 |
| 64 | 1.64e-9 | 1.51e-1 |
| 96 | 6.66e-15 (fp floor) | 6.85e-2 |

- Chebyshev: geometric collapse to the fp64 floor by N≈96 (measured converging-window semilog slope ≈ −0.17
  decade/node — a rate no algebraic scheme can match).
- FD2: measured log–log order **−1.977 ≈ O(N⁻²)** — precisely the incumbent's scaling, and precisely the ceiling
  D-021 measured.
- **Headline number for the tribunal:** to reach Chebyshev's N=64 accuracy (1.6e-9) the O(N⁻²) grid would need
  **≈622,000 points — about 9,700× the degrees of freedom.** *That* is why the CPU budget was never the real
  constraint; the discretization order was. A 1D radial spectral EKG at N ~ 10²–10³ is milliseconds per step.

The map (§1.4) is what carries this smooth-field result into the near-critical regime: it keeps the evolved field
*smooth on the fixed grid* (self-similar and bounded), so the spectrum stays exponentially convergent instead of
developing the ever-steeper gradient that both breaks FD and, on a *static* spectral grid, would Gibbs-ring — which
brings us to the flank.

---

## 3 · Guarantee-audit against the fixed constraints (where I strain, stated plainly)

| Constraint (Charter §"fixed constraints") | Verdict | Where it strains |
|---|---|---|
| **Determinism** — (params,seed)→byte-identical | **PASS for layer-1** (static map: RK4 + dense `D` + diagonal filter are all fixed deterministic linear algebra, single-threaded → bit-reproducible, exactly the incumbent's idiom). **CONDITIONAL for layer-2** (a feedback-rescaled map is deterministic *iff* the control law is one closed-form scalar update from a declared observable; a data-dependent multi-domain *re-split* would not be). | **Strain #1.** I must prove layer-2 stays a single scalar under a frozen law. If γ needs true multi-domain adaptivity whose *topology* responds to data, the DETERMINISM lens wounds the claim and I fall back to layer-1 + fixed domain decomposition. |
| **Golden-freezable** `(params)→blake2b` | **PASS.** The declared state is the spectral coefficient vector (or nodal values) at final `τ` → same blake2b idiom as `substrate_nexus`. The map parameter `L`/`ξ_max`, N, `dτ`, filter strength are all in the frozen config. | The filter and map params must be in the golden or the hash is meaningless — noted, mechanical. |
| **CPU, ≤ ~5 min, single-file** | **PASS with margin, expected.** N ~ 10²–10³ dense `D` apply is O(N²) ~ 10⁴–10⁶ flops/step; even 10⁵ steps is sub-minute. Domain decomposition adds a fixed block-tridiagonal solve (deterministic). Single-file C++17/fp64, stdlib-only, no GPU — identical build story to the incumbent. | If `D²` stiffness forces `dτ ~ N⁻⁴`, step count balloons; mitigated by keeping the system first-order (binding operator `D`, `dτ ~ N⁻²`) and/or implicit-splitting the stiff linear part. Costed, not hand-waved. |
| **Honest oracle** (anchor + metamorphic + redundant) | **PASS by construction.** Anchor: γ=0.374, Δ=3.44. Metamorphic: γ invariant under N, `L`, `dτ`, domain count (a spectral scheme that is *converged* must show γ flat as N rises — the incumbent's fatal test, which spectral is built to pass). Redundant: mass-scaling AND subcritical-curvature-scaling AND the DSS echo period Δ, three observables, must agree. | **Strain #2 (the real one).** Near a *genuine* singularity Chebyshev loses exponential convergence and Gibbs-rings — `cheb_convergence.py` measures exactly this on a C¹ cusp (algebraic order **−0.68**, no geometric regime, blake2b as above). So the *un-mapped* spectral EKG would fail near p\* just like FD. The entire bet is that the similarity map keeps the evolved field smooth enough (self-similar, `Δ`-periodic) that the spectrum never sees the bare singularity. **If the map cannot hold smoothness — if residual non-self-similar structure keeps sharpening on the grid — the ORACLE-HONESTY lens kills the γ claim.** This is my pre-registered failure mode (§4). |

**Net:** determinism and budget are strong; the crown risk is entirely whether the coordinate map holds the field
smooth on a fixed grid through the critical regime. That is a *measurable* question, which is the point of §4.

---

## 4 · Pre-registered falsifier + costed deciding experiment

### 4.1 The falsifier (registered before building)

> **SPECTRAL IS FALSIFIED IF:** with the similarity/compactified map in place, the Chebyshev EKG solver's recovered γ
> does **not** converge as N rises — specifically, if `|γ(N=2N₀) − γ(N=N₀)|` fails to shrink toward zero (the same
> non-convergence D-021 measured for FD), OR if the coefficient spectrum fails to decay geometrically at the
> resolutions where the echo is present (Gibbs plateau), OR if achieving `|γ−0.374|<0.05` at `R²>0.99` requires a
> map control law that cannot be reduced to a single deterministic scalar update (determinism kill). Any one of these
> ⇒ spectral is not the crown; harvest the map idea (§steal) and yield to double-null or the reframe.

This is a *strong* pre-registration: it stakes the claim on N-convergence of γ (the exact thing the uniform grid
could not do) and on spectrum decay (a spectral-specific, directly observable smoothness receipt).

### 4.2 The deciding experiment DE-γ-spectral (the tool contract)

**Tool:** `ekg_spectral` (single-file C++17, fp64, stdlib, no GPU; universal envelope
`[--params] --seed N [--json|--csv] [--selftest] [--golden]`; exit 0 pass / 1 declared-gate-fired / 2 error).

**Declared formulation:** §1 — CGL collocation in compactified similarity `ξ`, RK4 in `τ`, parity origin, exponential
filter, constraint-solved `a,α`, static map (layer-1) with declared `L`; layer-2 feedback available behind a declared
flag and only if layer-1 under-resolves.

**Declared state → golden:** final-`τ` Chebyshev coefficient vector (+ scalar diagnostics) → blake2b.

**Battery:**
- **X1 flat/smooth** — spectral derivative & evolution reproduce a known self-similar test profile to fp-floor
  (the `cheb_convergence` law, now inside the real tool as a selftest). *This is the anchor that a stranger runs cold.*
- **X2 p\* locate** — drive the search with **ORRERY `autotune --locate crossing`** against the collapse flag
  (see 4.3), replacing hand bisection (ledger F5). The metric is a fixed observable (e.g. `max_ξ 2m/r` crossing 1).
- **X3 γ fit** — supercritical sequence in `(p−p*)`; least-squares `ln M_BH` vs `ln(p−p*)`.
- **X4 metamorphic** — repeat X3 at N ∈ {40, 60, 90, 135} and `L ∈ {L₀, 2L₀}`; **require γ flat** (the D-021-killer).
- **X5 redundant** — subcritical `max_τ,ξ R ∝ (p*−p)^{−2γ}` and DSS echo period Δ; require agreement with X3.

**PASS bar (pre-registered, quantitative):**
`|γ_fit − 0.374| < 0.05` at `R² > 0.99`, **AND** γ stable across X4 (`spread < 0.03` over the N/L set — this is the
teeth), **AND** Δ = 3.44 ± 0.3, **AND** deterministic golden reproduced twice, **AND** wall-clock ≤ 5 min CPU.
Miss any ⇒ exit 1 (a real negative result, honestly logged), never faked to 0.374.

**Costed.** Build: ~1 file, ~600–900 lines (the incumbent is comparable); the dense `D`, RK4, filter, and
constraint solve are all textbook. Runtime: N≤135, O(N²) steps, ≤10⁵ steps ⇒ seconds–minutes, comfortably inside
budget. Risk-adjusted effort: **medium** — the numerics are standard; the *physics tuning* (getting the CSS source
terms and `T*` handling right) is the real cost, shared with the double-null and CSS forks.

### 4.3 ORRERY runs already executed (evidence-grade search, `[ARGUMENT-GRADE]` numbers)

The p\* search *mechanism* — locate a critical parameter against a pre-registered target, deterministically, with a
cited hash — is proven now with `autotune` (read-only ORRERY, CWD `C:\ORRERY`). These stand in for X2 before
`ekg_spectral` exists; the *located value* is an analytic target, not γ (labelled accordingly).

| run | command (abbrev) | located | verdict/exit | declared blake2b |
|---|---|---|---|---|
| anchor | `autotune.py --golden` | peak@0.37 (tool prints **GOLDEN OK**) | pass / 0 | `c79002f23cf236baab5ecdb5753603a7a3853199f750b0139771d4e6cdd55bbe` |
| crossing | `--objective threshold --obj-center 0.5 --obj-width 0.05 --lo 0 --hi 1 --locate crossing --level 0.5 --target 0.5 --tol 0.02 --seed 0` | **0.500000** | pass / 0 | `db490a3111b770996fcec05b9bc59981f35395980be5f1ad617101dd13552c80` |
| band@γ | `--objective peak --obj-center 0.374 --obj-width 0.05 --lo 0.2 --hi 0.55 --locate argmax --target 0.374 --tol 0.02 --seed 0` | **0.374014** (err 1.4e-5) | pass / 0 | `27c8cb4719f8a6e2dd455e386cd3da63ab1a9630896ab11a084a29c004d04eeb` |

`autotune --selftest`: 6/6 PASS, incl. "wrong pre-registered target → G-OFF-TARGET fires → exit 1" (the instrument
*can* fail honestly). The band@γ run is deliberately centered on the true 0.374 to show the locate resolves the
Choptuik value's neighborhood to 1.4e-5 — the search is not the bottleneck; the *physics tool feeding it* is, and
that tool is `ekg_spectral`, still to be built.

---

## 5 · Literature (verified; unsure items flagged)

- **Choptuik, M. W., "Universality and scaling in gravitational collapse of a massless scalar field," Phys. Rev.
  Lett. 70, 9–12 (1993).** — the origin; γ ≈ 0.37, discrete self-similarity (echoing), "structure on arbitrarily small
  spatiotemporal scales," AMR required to see it. *Verified: PRL 70, 9–12 (1993); link.aps.org/doi/10.1103/
  PhysRevLett.70.9.* This is the ground truth and the reason a static uniform grid cannot see γ.
- **Hamadé, R. S. & Stewart, J. M., "The spherically symmetric collapse of a massless scalar field," Class. Quantum
  Grav. 13, 497 (1996)** (arXiv gr-qc/9506044). — **CORRECTION to the brief's framing, stated honestly:** their
  method is **double-null (characteristic) coordinates with characteristic *adaptive mesh refinement*** (the first
  characteristic-grid AMR; both the standard Berger–Oliger algorithm and their own simplified version, giving
  indistinguishable results), confirming Choptuik's γ and echoing. It is **not** primarily a pseudospectral scheme.
  *I do not claim H&S as spectral precedent* — that would be a fabricated pedigree. I cite them as (i) independent
  confirmation of γ via a *characteristic* formulation (relevant to F1/double-null, not to me) and (ii) evidence that
  *coordinate choice + adaptivity* is the axis on which resolution is won. *Verified via Winicour, "Characteristic
  Evolution and Matching," Living Rev. Relativity (multiple eds.; the H&S characterization is explicit there),
  arXiv gr-qc/0102085; and the CQG/arXiv record.* Method detail `[CITE-UNVERIFIED-EXACT-QUOTE]` only in that I could
  not extract a verbatim methods sentence from the mangled PDF; the double-null-AMR characterization itself is
  corroborated by two independent review sources.
- **Boyd, J. P., "Chebyshev and Fourier Spectral Methods," 2nd revised ed., Dover, New York (2001)**
  (ISBN 9780486411835). — the spectral-methods authority underpinning §1: CGL collocation and differentiation
  matrices (§5), rational Chebyshev `TB` functions for semi-infinite domains (ch.17), exponential filtering (ch.9),
  and the exponential-convergence-for-analytic-functions / Gibbs-for-non-smooth dichotomy that §2–§3 rest on.
  *Verified: Dover 2001, 2nd revised ed.; full text also at depts.washington.edu/ph506/Boyd.pdf.*
- **Gundlach, C. & Martín-García, J. M., "Critical Phenomena in Gravitational Collapse," Living Rev. Relativity 10,
  5 (2007)** (arXiv gr-qc/0001046 is the earlier Gundlach 1999/2000 Living Review; the 2007 update is lrr-2007-5).
  — the review giving γ = 0.3737 (massless scalar), the DSS echo period Δ ≈ 3.44, the self-similar-attractor /
  codimension-one picture, and the similarity-coordinate (CSS/DSS) reduction that §1.1 uses. *Verified: Living Rev.
  Relativity 10 (2007), 5.*

**Additional genuine spectral / fixed-grid precedent (verified titles/venues; specific claims flagged):**
- **`bamps`** — a **nodal pseudospectral** code using **Gauss–Lobatto–Chebyshev** collocation with **domain
  decomposition (multiple patches)**, used to study the neighborhood of the spherical Choptuik critical solution and
  its non-spherical deformations (e.g. Rinne/Hilditch/Bugner et al.; "Aspherical deformations of the Choptuik
  spacetime," arXiv 1807.10342; "Twist-free axisymmetric critical collapse of a complex scalar field," arXiv
  2402.06724). — **This, not Hamadé–Stewart, is the real spectral-critical-collapse pedigree**: it proves Chebyshev
  pseudospectral methods *do* reach the Choptuik regime. *Code identity and Chebyshev-collocation characterization
  verified via search abstracts; exact per-paper method details `[CITE-UNVERIFIED]` pending full-text read.*
- **Garfinkle, D., "Choptuik scaling in null coordinates," Phys. Rev. D 51, 5558 (1995)** (arXiv gr-qc/9412008).
  — null/similarity coordinates adapted to the DSS "beautifully avoid the need for mesh refinement," with "perfect
  numerical tuning in double-precision arithmetic in spherical symmetry." *Directly supports layer-1 of my map (§1.4).
  Venue PRD 51, 5558 (1995) `[CITE-UNVERIFIED-PAGE]`; the arXiv id gr-qc/9412008 and the "no mesh refinement" claim
  are verified via review search.*
- **"Type II critical collapse on a single fixed grid: a gauge-driven ingoing boundary method," arXiv 2008.12726.**
  — reproduces Choptuik mass-scaling, DSS, and universality on a **single fixed grid, no AMR**, by a spatial-gauge
  choice that "concentrates resolution in the central region as collapse occurs" while "the numerical grid remains
  fixed... no points need to be removed or added." *Independent proof that the fixed-grid + coordinate-concentration
  strategy (my thesis) reaches γ.* *Title/venue verified via abstract; γ/Δ values `[CITE-UNVERIFIED]`.*

---

## 6 · Honest risks + STEAL list

### Risks (ranked; the first is the crown risk)

1. **The map may not hold smoothness through the critical regime (the ORACLE-HONESTY kill).** If residual
   non-self-similar content keeps sharpening on the fixed grid, the spectrum stops decaying geometrically and Gibbs
   returns — `cheb_convergence.py` shows Chebyshev's algebraic −0.68 order on a bare C¹ cusp, i.e. *un-mapped spectral
   is no better than FD near a true singularity.* The whole crown rests on the similarity map keeping the evolved
   field analytic-on-grid. **This is measurable (X4 spectrum decay) and pre-registered as a falsifier — it is where I
   most expect to be attacked, and I would rather be killed here honestly than fake γ.**
2. **Layer-2 adaptivity vs determinism (the DETERMINISM strain).** If `T*` drift forces more than a single
   closed-form scalar map update, the "spectral is trivially deterministic" advantage erodes toward the AMR fork's
   problem. Mitigation: keep to layer-1 + *fixed* domain decomposition; accept a wounded claim if layer-2 is
   unavoidable.
3. **Conditioning / stiffness.** `‖D²‖~O(N⁴)`; if the formulation exposes it, `dτ` shrinks. Mitigation: first-order
   system (binding `D`, `N⁻²`), implicit-split the stiff linear part (deterministic dense LU), or cap N via domains.
4. **"Differently hard, not simpler" (the brief's own challenge).** CSS source terms, `T*` determination, and DSS
   log-periodic bookkeeping are genuine complexity — spectral is not obviously *simpler* than double-null. Honest
   position: spectral is not simpler to *derive*; it is dramatically cheaper to *resolve* (§2's 9,700× DOF), and it
   offers a spectrum-decay smoothness receipt no FD/AMR scheme gives. I sell resolution-per-DOF and a self-diagnosing
   convergence, not derivational simplicity.
5. **Precedent honesty.** My strongest pedigree is `bamps` and the fixed-grid gauge paper, **not** Hamadé–Stewart
   (whom I explicitly decline to misclaim). A tribunal that expected H&S-as-spectral should note I corrected it.

### STEAL list (harvest regardless of my fate — Charter §steal)

- **The exponential filter as principled dissipation** — a declared diagonal coefficient-space operator; cleaner and
  more tunable than KO dissipation; *any* fork (even the FD incumbent) can adopt it as a better-behaved dissipator.
- **Spectrum-decay as an honesty gate** — "is the coefficient tail geometric?" is a *free, self-diagnosing*
  resolvability receipt. Bolt it onto the double-null or AMR tools: if their solution is truly resolved, a spectral
  transform of it must show geometric decay. This is a metamorphic check any winner should carry.
- **`autotune --locate crossing/argmax` as the p\* search for every fork** (ledger F5) — replaces the hand bisection
  that let p\* wander; deterministic, pre-registered target, cites its own blake2b. Proven now (§4.3). Steal into the
  winning tool whatever the formulation.
- **The similarity/compactified map itself** — even if spectral loses, "evolve in coordinates where the DSS solution
  is log-periodic so a *fixed* grid suffices" (Garfinkle; arXiv 2008.12726) is the deepest determinism story on the
  board and belongs to whichever formulation ships.
- **The 9,700×-DOF convergence table** — a portable, citable argument (blake2b `b6ad2eba…`) for *why* fixed-order
  uniform grids hit D-021's ceiling; useful in the phase-3 adjudication regardless of the chosen scheme.

---

### Reversibility lemma (house discipline)

If spectral is deleted in phase 3, its **pre-registered reinstatement trigger** is: *any evidence that the winning
formulation's solution, when spectrally transformed, shows a geometric coefficient tail through the near-critical
regime* (i.e. the field is analytic-on-grid there) — because that is precisely the condition under which spectral's
9,700×-DOF advantage becomes decisive and its only crown risk (risk #1) is retired. Re-arm spectral the moment that
smoothness receipt appears.

*Every γ figure herein is `[ARGUMENT-GRADE]`; the only evidence-grade artifacts are the convergence law
(blake2b `b6ad2eba341d427ffa8c825058cb6a383c4500660fabbf3f576b9920d63dcaaf`) and the three ORRERY autotune locates
(hashes in §4.3). The spectral EKG instrument that would make γ itself evidence-grade is specified in §4.2 and does
not yet exist.*
