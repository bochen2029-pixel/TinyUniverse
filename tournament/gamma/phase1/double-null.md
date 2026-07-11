# APPROACH — double-null / characteristic (phase 1, the γ fork)

**Persona:** the double-null / characteristic advocate (F1 attack stance).
**One-line thesis:** the uniform grid failed *for polar-areal* because polar-radial slicing crowds the resolution it needs into a coordinate region the grid cannot follow; **change the coordinates so the grid rides the imploding pulse and the self-similar structure stays resolved at fixed cost.** D-021 is a *gauge* artifact, not a grid law.
**Central exhibit:** Garfinkle (1995) recovered Choptuik's critical phenomena — including the discrete self-similarity — from a **null-coordinate code that does *not* use adaptive mesh refinement**. That is an existence proof that a *non-adaptive* grid can see this physics, provided the coordinates are null.

> **Scope honesty up front.** The real double-null EKG integrator does not exist yet. Every γ figure below is labelled **`[ARGUMENT-GRADE]`**; the only evidence-grade numbers in this doc are the ORRERY `autotune` runs (§4, §ORRERY), which certify the *search* mechanism, not γ. The deciding experiment (§4) is what would move the crown claim to evidence-grade.

---

## 1 · The approach, fully specified

### 1.1 Geometry & coordinates
Spherical symmetry, one radial null structure. Two standard, near-equivalent formulations; I name the primary and the fallback and pick the primary.

- **Primary — Bondi / single-null (Misner–Sharp, "outgoing-null" a.k.a. the Christodoulou/Garfinkle line).** Coordinates `(u, r)` with `u` a retarded (outgoing-null) time and `r` the areal radius. Metric
  `ds² = −g(u,r) ḡ(u,r) du² − 2 g(u,r) du dr + r² dΩ²`,
  scalar field `φ(u,r)`, auxiliary `h ≡ ∂_r(rφ)` (the natural characteristic variable). Einstein's equations reduce to **radial ODEs along each `u`-slice** for the metric functions `g, ḡ` (hypersurface/constraint equations, integrated outward in `r` from the regular origin), plus **one evolution equation** advancing `φ` from slice `u` to `u+du`. This is the Goldwirth–Piran / Garfinkle structure. *(This is the single-null Misner–Sharp / Hamadé–Stewart variant the charter names.)*

- **Fallback — double-null.** Coordinates `(u, v)` both null,
  `ds² = −α²(u,v) du dv + r²(u,v) dΩ²`,
  `r` and `α` are evolved; the field obeys a wave equation that in these coordinates is a simple `∂_u ∂_v (rφ) = …` diamond-update. Marching is on a null diamond grid `(u_i, v_j)`.

Both share the decisive property (§2): **the grid is laid along the characteristics of the field**, so an imploding pulse does not "run off" a fixed radial grid — it stays sampled by the same rays as it focuses.

### 1.2 The characteristic marching scheme
- **Primary (single-null, marching in `u`):** on each outgoing slice `u = uₙ`, given `φ(uₙ, ·)`, integrate the hypersurface ODEs **outward in `r`** (a clean radial quadrature — the same elegance polar-areal has for its constraints) to get `g, ḡ`, hence the characteristic speed field. Then take one **explicit step in `u`** for `h ≡ ∂_r(rφ)` using the reduced wave equation, RK-in-`u` with 2nd-order centred `r`-differences. Origin regularity is imposed at `r=0` by the local expansion `rφ → 0`, `∂_r(rφ)` finite. Outer boundary: excise inside the apparent horizon (see §3) and track the mass aspect there.
- **Fallback (double-null diamond):** the field update is the exact null-diamond rule
  `(rφ)_{i+1,j+1} = (rφ)_{i+1,j} + (rφ)_{i,j+1} − (rφ)_{i,j} − ∫∫ (source) du dv`,
  with `r, α` advanced by their constraint/evolution pair on the same diamond. Second-order, unconditionally causal (the stencil *is* the light cone), trivially deterministic.

Fully deterministic fp64, single-threaded ⇒ `(params) → byte-identical declared JSON → blake2b` (the tiny_nexus / substrate_nexus idiom, lifted verbatim).

### 1.3 The observable (the crown, and its redundancy)
- **Primary observable — mass scaling.** Bisect (via ORRERY `autotune`, §4) the initial amplitude to `p*`. For a supercritical sequence `p > p*`, the black-hole (apparent-horizon) mass obeys the Choptuik law
  `M_BH ∝ (p − p*)^γ`,  fit `ln M_BH` vs `ln(p − p*)`; slope = γ.
- **Redundant recovery #1 — subcritical curvature scaling.** For `p < p*`, the maximum central curvature (equivalently `max |∂φ|²` at the origin) scales as
  `max R ∝ (p* − p)^(−2γ)`  (Garfinkle–Duncan). A **second, independent** estimate of γ from the *other* side of the critical point — the DE-redundant gate demands they agree.
- **Redundant recovery #2 (qualitative, cheapest) — the echo.** In similarity time `τ = −ln(p* − p)` (or along the null cone approaching the accumulation point) the critical solution is **periodic with period Δ ≈ 3.44** (discrete self-similarity). Detecting the log-periodic echo is a *far coarser-grid* signature and is Garfinkle's actual reported result. It is the fallback crown if the exponent stays noisy (see §6, STEAL).

**Confirming γ stays the crown.** γ ≈ 0.374 is the number the null literature *actually recovered on a non-adaptive grid* (Garfinkle 1995). Unlike the reframe stance, this approach does **not** concede the crown — it argues the crown was lost to a coordinate choice, and the fix is a coordinate change, not a scope change.

---

## 2 · Why this beats the polar-areal incumbent D-021 killed

The incumbent (`substrate_nexus.cpp`, polar-areal, uniform radial grid) is honest and correct about *its own* failure. The failure is **diagnosably coordinate-specific**, and the diagnosis is the whole argument:

1. **Where the resolution has to be, and where the grid puts it.** Critical collapse concentrates all structure into a shrinking central region over successive echoes; the physical scale to resolve shrinks like `e^{−nΔ/2}` per echo. On a **uniform polar-areal `r`-grid**, cells are fixed in areal radius `r`. The self-similar feature slides *through* the grid toward `r=0` and, crucially, the **central lapse `α(0) → 0`** as the slice hangs up outside the forming horizon (polar slicing is singularity-avoiding — a feature for stability, a bug for resolution): coordinate time grinds nearly to a halt exactly where the physics is happening, so ever more evolution buys ever less proper progress. The grid caps the resolvable curvature at ≈`(φ/dr)²` — *measured* in D-021.

2. **Why refining doesn't converge (the chaos in D-021).** Halving `dr` on the polar-areal grid does not track the feature; it just moves the same fixed wall inward a little while the near-critical solution — which is exponentially sensitive to `p` near `p*` — is re-sampled differently. The bisected `p*` wandering 0.40 → 0.356 from N=800 → 1600 is the fingerprint of a *search over a grid-dependent classifier*, not of true `p*` motion. (Note this indicts **both** the grid and the hand-bisection classifier — see §6 STEAL and F5.)

3. **What null coordinates change.** The grid is laid along the **characteristics of the field**. An imploding shell of scalar radiation follows an ingoing null ray; in `(u, r)` or `(u, v)` the numerical rays **advect with the pulse**, so the same computational cells that carried the pulse at large `r` carry it as it focuses. The self-similar profile is (nearly) *static in the null-adapted frame that co-moves with the accumulation point* — periodic in `τ`, not translating through the grid. Hence a **fixed, non-adaptive** null grid keeps the feature resolved at fixed cost. This is not a hope: it is exactly what Garfinkle's non-AMR null code did.

4. **The clean part is preserved, the broken part is replaced.** The elegance polar-areal has — constraints as outward radial ODEs — is **retained** in the single-null formulation (hypersurface equations integrated outward in `r`). We keep the good structure and drop the coordinate that hides the collapse.

**Precise claim vs the graveyard corpse.** The graveyard holds *uniform-ultrafine polar-areal* (D-021). This approach is **not** a re-proposal of that corpse: it is a *different coordinate system* whose non-adaptive grid is placed along null rays. The corpse reason ("a fixed grid caps the self-similar central curvature; refining goes chaotic") is a statement about **spacelike-radial** cells hung on a singularity-avoiding slice. It does not transfer to null cells that ride the pulse — the reinstatement-trigger logic (below) is exactly a test of that transfer.

---

## 3 · Guarantee-audit against the CHARTER's fixed constraints

The box (CHARTER §"fixed constraints"): deterministic golden · golden-gated & honest · CPU ≤ ~5 min single-file · exit codes 0/1/2.

| constraint | status | where I am at risk (stated plainly) |
|---|---|---|
| **Determinism** (`(params)→byte-identical blake2b`) | **Structurally safe.** Explicit fp64, single-threaded, no RNG, no atomics; the null-diamond / RK-in-`u` update is a fixed arithmetic sequence. Same idiom already frozen 18× in this project. | Low. The only risk is transcendental-libm drift across compilers — mitigated exactly as substrate_nexus does (host libm pinned, golden re-run out-of-process). |
| **CPU ≤ ~5 min, single-file** | **Plausible with margin.** Garfinkle-class runs were cheap on 1995 hardware; a single-file C++17 fp64 characteristic integrator at a null resolution sufficient for γ is well inside 5 min on this machine. One evolution is O(N²) in a diamond of side N; a p*-bisect + a ~12-point scaling sequence is ~30–50 evolutions. | **Medium.** If the null grid needs to be *very* fine near the accumulation point to pin γ to tol, the diamond cost can grow. Costed in §4; the fallback is fewer echoes + the echo-Δ crown (§6). |
| **Golden-gated & honest** (anchor + metamorphic + redundant) | **Designed in.** Anchor = γ ≈ 0.374 / Δ ≈ 3.44 (literature). Metamorphic = γ invariant under null-grid refinement and under `(r0, σ)` of the initial pulse (it *must* be, being universal). Redundant = mass-scaling **and** curvature-scaling must agree (DE-redundant). | **This is the real battleground.** See the three named risks below the table. |
| **Exit codes 0/1/2** | **Trivial.** Same plumbing as substrate_nexus: 0 all gates pass, 1 a declared gate fired (e.g. `|γ−0.374| ≥ 0.05` → a real negative result, honestly reported), 2 error. | None. |

**The three places I am genuinely at risk (oracle-honesty specifics):**

- **R-a · Null-grid caustics / focusing breakdown.** Ingoing characteristics can cross (caustics) as the field focuses; the single-null `(u, r)` chart breaks down where `∂r/∂(grid) → 0`. Near-critical, the slice must be **terminated cleanly at the apparent horizon** (excision) *before* the caustic, and the mass read there. If the code instead NaNs at the caustic and that gets misclassified as "collapse," the classifier is corrupted — the exact D-021 failure mode in new coordinates. **Mitigation & gate:** the collapse criterion is a *resolved* apparent horizon (`2m/r → 1` with the horizon ≥ k cells wide), never "the run blew up"; runs that caustic before a resolved horizon are **excluded, not misclassified** (and counted — if too many are excluded, the fit is declared unreliable, gate fires).

- **R-b · The r→0 regularization.** Both charts need origin regularity (`rφ → 0`, `∂_r(rφ)` finite). A sloppy origin condition injects a spurious central source that *mimics* curvature growth — poisoning the very observable (curvature scaling) that is supposed to cross-check mass scaling. **Mitigation & gate:** origin handled by the analytic local expansion, and the **S1 flat-space oracle** (a weak pulse must disperse with the mass aspect conserved to `<1e-3` and `α`/`a`≈1) is retained verbatim as the regularization's canary — it fails loudly if the origin leaks.

- **R-c · Refinement of the null grid near collapse.** My whole pitch is "non-adaptive suffices" (Garfinkle). But Hamadé–Stewart (1996) needed **characteristic AMR** to get γ to *high* accuracy. So there is a real possibility that **non-adaptive null resolves the echo and a rough γ, but not γ to the CHARTER's `|γ−0.374|<0.05`**. If so, the honest outcome is: the *transition + echo* is a clean non-adaptive null crown, and *tight γ* still wants (characteristic) AMR — which would **wound, not kill**, and hands a STEAL to the AMR advocate. This is pre-registered in §4's falsifier.

---

## 4 · Pre-registered falsifier + costed deciding experiment (DE-γ)

### 4.1 What result kills me (pre-registered, before building)
> **KILL:** a single-file, deterministic, ≤5-min **non-adaptive** double-null / single-null EKG tool, driven by `autotune` to locate `p*`, yields a mass-scaling fit with **either** `|γ_fit − 0.374| ≥ 0.05` **or** `R² < 0.99`, **and** the failure does **not** improve toward 0.374 as the null grid is refined (i.e. it is not merely under-resolved but *converging wrong* or *not converging*). If refinement *does* march γ_fit toward 0.374 but can't reach tol within the 5-min budget, that is a **WOUND** (non-adaptive null resolves the phenomenon but not the exponent at N0 cost — the crown reframes to "transition + echo," tight γ → characteristic-AMR contract), **not** a kill.

> **Secondary KILL (oracle-honesty):** mass-scaling γ and subcritical curvature-scaling γ **disagree by > 0.03** on the converged grid (DE-redundant) — that means one observable is grid-noise dressed as a result.

### 4.2 The exact single-file EKG characteristic-integrator contract to build
`contracts/ekg_null.contract.md` v0.1.0 (semver from 1.0.0 on freeze). **Single file** `experiments/double-null/ekg_null.cpp`, C++17, fp64, stdlib only, no GPU.

- **Formulation:** single-null Bondi/Misner–Sharp primary (§1.1), metric `ds² = −g ḡ du² − 2 g du dr + r² dΩ²`, `φ(u,r)`, characteristic var `h = ∂_r(rφ)`; hypersurface ODEs integrated outward in `r`; explicit RK-in-`u`; 2nd-order centred `r`-differences; Kreiss–Oliger dissipation `ε_KO` declared; origin by local expansion; apparent-horizon excision `2m/r → 1`.
- **Units:** geometric `G = c = 1` (γ dimensionless — inherit substrate_nexus's units decision, D-020).
- **Battery (mirrors substrate_nexus S1/S2/S3, replaces S4):**
  - **N1** flat-space wave — weak Gaussian disperses, mass aspect conserved `<1e-3`, no horizon *(the R-b origin canary)*.
  - **N2** subcritical — disperses, no horizon; record `max R` for curvature scaling.
  - **N3** supercritical — resolved apparent horizon forms, `0 < M_BH < M_ADM`.
  - **N4 (crown)** critical scaling — `autotune`-located `p*`; mass-scaling γ over a supercritical sequence; **PASS = `|γ_fit − 0.374| < 0.05` at `R² > 0.99`**; redundant curvature-scaling γ agrees `< 0.03`; **plus** echo period `Δ = 3.44 ± 0.3` detected (log-periodicity of `max R`/`ln(p*−p)`).
- **Envelope:** `--json`/`--golden`/`--selftest`/`--only Nk`/`--dev p`; exit 0/1/2; declared JSON → blake2b; frozen config `(N_null, u_max, r_max, ε_KO, λ)` pinned when green.
- **Determinism gate:** battery run twice in-process must serialize identically; golden re-run out-of-process is the stronger check.

### 4.3 The `autotune` run that locates p* (the search, already de-risked)
The p*-search is **not** hand bisection (that was half of D-021's instability, per F5). It is a pre-registered level-crossing:
```
python C:\ORRERY\tools\autotune\autotune.py \
  --tool <path>\ekg_null.exe --sweep p --metric result.collapsed \
  --fixed "--N_null <N> --u_max <U> --r_max <R>" \
  --lo <p_lo> --hi <p_hi> --points <K> \
  --locate crossing --level 0.5 --target <p*_pre-registered> --tol <t> --seed 0 --json
```
`ekg_null.exe` emits `result.collapsed ∈ {0,1}`; `autotune` locates the 0.5-crossing = `p*`, and **G-OFF-TARGET fires (exit 1) if it lands off the pre-registered p***, so the search is self-certifying. *(This is exactly the mechanism I proved deterministic today with a stand-in objective — see §ORRERY and RESULTS.md.)*

### 4.4 Cost
One evolution: O(N_null²) diamond, N_null ~ few×10³ ⇒ ~10⁷–10⁸ fp64 ops ⇒ sub-second to seconds. `autotune` p*-locate at K≈21 points + a 12-point scaling sequence + curvature sweep ⇒ ~50 evolutions ⇒ **target a few minutes**, inside the 5-min box with the fallback (fewer echoes) if focusing cost bites (R-c).

---

## 5 · Literature (cited precisely; unverifiable items flagged)

- **Choptuik, PRL 70, 9 (1993)** — "Universality and scaling in gravitational collapse of a massless scalar field." Discovery of Type-II critical collapse: universal `p*`, mass scaling `M_BH ∝ (p−p*)^γ` with **γ ≈ 0.37**, and discrete self-similarity (echoing). *(Journal ref per the assumption ledger / charter; the numeric γ and echoing are the field's ground truth.)*
- **Garfinkle, PRD 51, 5558–5561 (1995)** — "Choptuik scaling in null coordinates" (arXiv gr-qc/9412008). **VERIFIED verbatim abstract:** *"The algorithm uses the null initial value formulation of the Einstein-scalar equations, but does not use adaptive mesh refinement… it is verified that the critical solution exhibits periodic self-similarity."* **My central exhibit:** a non-adaptive null code recovers Choptuik's critical phenomena. (Corroborating secondary sources describe it as a modified Goldwirth–Piran double-null code recovering the mass-scaling parameters and the echoing.)
- **Hamadé & Stewart, CQG 13, 497 (1996)** — "The spherically symmetric collapse of a massless scalar field." **Double-null** formalism, but developed the **first characteristic grid *with* adaptive mesh refinement** to reach the accuracy needed to confirm Choptuik's exponents. **Honest use:** they are the double-null *pioneers*, and evidence that *high-accuracy* γ in null coordinates may want AMR — which is precisely my R-c risk, not a point in my favour. `[CITE-UNVERIFIED: the specific γ value they quote — IOPscience returned 403; do not attribute a number to them until the paper is read.]`
- **Gundlach & Martín-García, Living Rev. Relativity 10, 5 (2007)** — "Critical phenomena in gravitational collapse" (arXiv gr-qc/0001046 / 0711.4620). The authoritative review; source for **γ = 0.3737** and echo period **Δ ≈ 3.44** for the 4D massless real scalar. *(Values per the assumption ledger, which cites this review; the PDF is binary and I could not re-extract the exact digits in-session — the anchor digits are `[CITE-UNVERIFIED against the primary PDF this session]` but are the project's already-adopted ground truth.)*
- **Garfinkle & Duncan, PRD 58, 064024 (1998)** — subcritical curvature scaling `max R ∝ (p*−p)^(−2γ)`, the basis of my redundant-recovery observable. *(Journal ref per the assumption ledger F3; not independently re-verified this session — `[CITE-UNVERIFIED]`.)*

I did **not** fabricate any γ or Δ figure. The only quantitatively-certain numbers I produced this session are the ORRERY `autotune` hashes below.

---

## ORRERY-arming — evidence for the *search*, honestly scoped

**Tool:** `autotune` v1.0.0 (read-only). Provenance: `--selftest` → SELFTEST PASS; `--golden` → `GOLDEN OK blake2b=c79002f23cf236baab5ecdb5753603a7a3853199f750b0139771d4e6cdd55bbe`. Driver + full record: `../experiments/double-null/` (`run_autotune.py`, `RESULTS.md`).

| run | purpose | outcome | declared blake2b |
|---|---|---|---|
| **charter-arming** (`--objective threshold --obj-center 0.5 --obj-width 0.05 --lo 0 --hi 1 --points 41 --locate crossing --level 0.5 --target 0.5 --tol 0.02 --seed 0`) | the charter's exact arming command: locate a level-crossing at a pre-registered target | `x_located=0.500000`, on_target, **exit 0**; re-hashes byte-identical | `db490a3111b770996fcec05b9bc59981f35395980be5f1ad617101dd13552c80` |
| **DN-locate-pstar** (threshold sigmoid @0.334, crossing, target 0.334, tol 0.005) | a Choptuik-shaped p*-locate: collapse fraction steps 0→1 across a pre-registered p* | `x_located=0.334000`, on_target, **exit 0** | `26ec96dea28874a60730c15628eea966d610667e5a56beae5074617cd8c25ac5` |
| **DN-locate-wrongtarget** (same curve, target 0.360) | **negative control / KILL guard**: the search must refuse a fabricated p* | `on_target=false`, **G-OFF-TARGET fired, exit 1** | `371b72f65d9a25b200886890d5f9bdc131f53302bddcccf5307394c79fd52838` |

**What this proves (and what it does not).** It proves the DE-search mechanism — `autotune --locate crossing` against a **pre-registered** target — is deterministic, byte-reproducible, and **self-certifying**: it fires a gate and exits 1 when handed the wrong critical point, so it cannot silently rubber-stamp a fabricated `p*` (the honest-oracle discipline made mechanical). It does **not** prove γ, and the sigmoid is a stand-in — `obj-center 0.334` is *not* a measured amplitude and *not* a γ figure. The real double-null EKG tool (§4.2) is what turns the crown claim evidence-grade. Any γ figure in this doc is **`[ARGUMENT-GRADE]`**.

---

## 6 · Honest risks + STEAL list

### Risks (fair-minded, worst-first)
1. **R-c is the one that could sink me:** non-adaptive null may resolve the *echo* and a *rough* γ but not γ to `<0.05` at N0 cost — Hamadé–Stewart needed characteristic AMR for accuracy. Outcome would be a WOUND (crown reframes to transition+echo; tight γ → AMR contract), not a clean win. **This is the honest most-likely partial outcome.**
2. **R-a caustics** could reintroduce a D-021-style misclassification if the collapse criterion is "it blew up" rather than "a resolved horizon formed." Mitigated by excision + a width gate, but it is real engineering risk near `p*`.
3. **R-b origin leak** could poison the curvature cross-check (the redundancy that makes γ honest). Mitigated by the analytic origin expansion + the retained S1 canary.
4. **F1↔F3 entanglement (ledger A1):** double-null may *make mass-scaling the natural clean observable* while CSS makes the echo natural — so the "which observable" question may not be independent of this fork. If curvature-scaling is unexpectedly cleaner than mass-scaling in null coords, I adapt the crown observable accordingly (still γ).
5. **Budget vs crown (ledger A4):** if only very fine null grids reach tol and they exceed 5 min, the operator must price budget-vs-crown — I flag it rather than silently overrun.

### Reversibility lemma (house discipline — the reinstatement trigger for the corpse this touches)
This approach does not delete a corpse; it *challenges* the D-021 corpse's scope. **Pre-registered reinstatement trigger for "non-adaptive null resolves γ":** if the built `ekg_null` tool shows γ_fit **not converging** toward 0.374 under null-grid refinement (the R-c KILL branch), then "non-adaptive suffices for the *exponent*" is retired, the *transition+echo* result is kept as the non-adaptive null crown, and **tight γ is handed to a characteristic-AMR contract** — with this doc's §2 diagnosis preserved as the reason AMR should be *characteristic*, not spacelike. Conversely, if it *does* converge to tol, D-021's deferral is overturned for the null gauge (the corpse stays valid *only* for polar-areal, as originally measured).

### STEAL list (harvest regardless of this approach's fate)
- **S-1 · `autotune`-located p*, not hand bisection.** The self-certifying level-crossing search (proven deterministic today, exit-1 on a wrong target) should replace substrate_nexus's 20-iter hand bisection **in every approach** — it removes the search-half of D-021's instability (F5) whether or not double-null wins.
- **S-2 · The resolved-horizon collapse criterion + exclusion accounting.** "Collapse = a resolved apparent horizon ≥ k cells wide; caustic/NaN runs are excluded and *counted*" is a clean honesty rule any formulation (AMR, CSS, spectral) should adopt.
- **S-3 · Redundant γ (mass-scaling ∧ curvature-scaling ∧ echo-Δ) as the oracle gate.** Requiring two independent observables to agree to `<0.03` is the metamorphic discipline that stops grid-noise-as-γ — steal it into whichever substrate wins.
- **S-4 · The echo-Δ crown as a *cheaper, robust* deepest-GR receipt.** Even if the reframe advocate wins on the exponent, "we see the period-3.44 discrete self-similarity" is a coordinate-agnostic receipt worth banking; null coords make it especially clean (Garfinkle's actual result).

---

*Thesis in one breath: D-021 killed a **gauge**, not the crown. Garfinkle recovered this physics on a non-adaptive **null** grid; put the grid on the light cone so it rides the collapse, drive p\* with a self-certifying `autotune` crossing, and gate γ against two independent observables. Build `ekg_null.cpp`, freeze its golden, and let the number — not the argument — decide.*
