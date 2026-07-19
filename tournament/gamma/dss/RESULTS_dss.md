# RESULTS — scalar-field DSS crown (γ = 0.374 via Δ-periodic critical solution)

*The γ-thread's canonical results file (analog of `../phase4/RESULTS_hka_beta.md`).
Contract: `contracts/similarity_nexus.contract.md` v0.1.0 (operator-approved 2026-07-19).
Primary source: Gundlach gr-qc/9604019 (`../GUN96_src/chop2.tex`, gitignored).*

## 2026-07-19 · Session 1 (autonomous, post-fluid-β-close) — the machine built, honest walls named

**Target ladder:** Stage-A = the DSS background (fields on the (ζ, τ-periodic) cylinder;
**Δ = 3.4453 ± 0.0005** emerges as an eigenvalue) → Stage-B = Floquet perturbations
(λ₁ = −2.674 ± 0.009 → **γ = 0.374 ± 0.001**), then the C++ `similarity_nexus` port.

### Built + validated (assets, all committed)

1. **`nr_evolve.py` — near-critical collapse evolver** (polar-areal massless scalar,
   Gundlach-compatible variables; the lapse CANCELS in X = √(2π)(r/a)Φ and
   Y = √(2π)(r/a)Π; central proper time t_G accumulated for τ = ln(t\*−t_G)).
   Debugged by measurement: the spurious `+wΠ/2r` origin term (Φ's equation is a plain
   gradient) caused a self-sustaining origin blowup; Π's equation moved to the
   uniformly-accurate ∂_r(wΦ) + 2wΦ/r form; polar-slicing freeze detected at 2m/r>0.9 ∨
   α(0)<0.02; reflection-safe dispersal window. **The central mode φ(0,t) must be EVOLVED**
   (∂_tφ|₀ = (α/a)Π|₀) — reconstructing φ from ∫Φ pins φ(0)≡0 and erases the echo signal
   (measured 0.015 → 0.50 swing after the fix).
   **Measured:** p\* = 0.03751655962597 (bisection to rel 4.5e-14, N=800);
   t\* = 14.913 ± 0.05; **in-house Δ_echo ≈ 3.2** (~1.2 grid-resolvable periods — the
   D-021 uniform-grid resolution physics; honest bound, not a claim);
   max|φ(0)| ≈ 0.50 (the Choptuik universal amplitude scale in these units).
2. **Seed pipeline** (cylinder sampling in the gauged frame + τ-gauge rotation + harmonic
   projection): seed RMS(Y₁) = 1.8–2.2, independently agreeing with the shooting-side
   amplitude-continuation dip at RMS(Y₁) ≈ 1.9 — two methods, same amplitude scale.
3. **`nr_dss.py` — two-sided-shoot BVP** (17× fast pseudo-Fourier residual, exact
   dense-DFT matmuls; boundary ladders per app A with the D₂ term).
   **VERDICT (measured, structural):** shooting is hostile here — off the solution
   manifold the a-constraint (f = a⁻², linear closed form) has NO real solution
   (DAE off-manifold non-existence) ⇒ cliff landscapes defeat LM regardless of softening
   (sloped cliffs, graded SSH steps, clamps — all tried). The SSH X₋-ladder's DC is a
   **parity SYMMETRY ZERO** (ε-regularizing it breaks the SSH regularity condition —
   measured). This is WHY Gundlach used global relaxation.
4. **`nr_relax.py` — the paper's global relaxation** (harmonic-coefficient unknowns per
   ζ-node, implicit midpoint, center BCs as pointwise node-0 relations, SSH conditions at
   ζ=0 exactly; the system closes SQUARE: 31·Nz+11). **Vacuum control: |r| = 0.0 exactly.**

### The two vacuum-watershed catches (RAYFORMER, same-day discipline)

- **The unpinned relaxation converged to |r| = 2e-28 with "Δ = 3.425, Nz-independent" —
  IT WAS THE VACUUM** (all field harmonics 1e-16; Δ = meaningless path leftover; vacuum
  solves the system for any Δ). **The tell: machine-zero residual at finite truncation.**
  A genuine solution floors at the harmonic-truncation scale (~1e-4 at KE=10); 1e-28 means
  truncation-exact = trivial. (Same error class as the fluid crown's Friedmann-as-EC,
  caught in-session by the spectral-tail check now specced in the contract.)
- The **amplitude-pinned** solve (Gundlach's own normalization device as an appended
  residual) at the seed's c = 0.072 lands NONTRIVIAL but **weak-field-perturbative**
  (gdev ~ 0.002, r_phys ~ 1.6e-2 ~ amplitude²): the seed's SSH X₊ extraction
  underestimates the true amplitude several-fold.
- **Pin-continuation upward in c** (0.10 → 0.55): Δ is WEAKLY DETERMINED at low amplitude
  (wandered 3.40 → 3.32 → slammed the 3.7 search bound, poisoning warm starts); c ≥ 0.38
  solves diverge. The period only locks at the nonlinear amplitude — carrying a free Δ
  through the weak-field regime is fragile ⇒ freeze Δ at the literature value during
  continuation, release at the end.

### State + next steps (for the next session)

- **Δ NOT measured, γ NOT measured, none faked.** The machinery (evolver + seed + both BVP
  formulations + the vacuum/weak-field diagnostics) is real, validated, committed.
- Next levers, in order: (1) **Δ-frozen amplitude continuation with full Newton budgets**
  (probe below); (2) **richer harmonics** (KE 10→14+, M 32→48; Gundlach used content up to
  k ~ 21 — the strong-field branch may be deformed/absent at KE=10); (3) **better SSH-side
  seeding** (the cylinder's X₊₀ extraction: iterate the ξ₀-frame; sample deeper echoes at
  N=1600 with a re-bisected p\*); (4) two-parameter arclength continuation (c, Δ) if the
  frozen-Δ path misbehaves.

### Probe verdict (same session): TRUNCATION-BINDING confirmed

Delta-frozen (3.4453) pinned (c=0.25) full-budget solve at (M,KE,KO) = (32,10,9):
r_phys plateaus at 0.029 and the **X+ harmonic tail RISES into the truncation edge**
([0.13, 0.20, 0.21, 0.22, 0.44] for k = 1,3,5,7,9) — the strong-field solution's spectral
content piles up against KO = 9. The binding constraint is the tau-truncation, not the
optimizer (Gundlach's own content reaches k ~ 21 at 2N = 128). Next: (48,14,13), then
(64,20,21) if the k=13 edge still rises.

### (48,14,13) continuation verdict: the OPTIMIZER is now the limiter (diagnosis inverted)

Delta-frozen pinned continuation at the richer truncation: every step stalls at xtol within
13-43 s (trust region collapse; r_phys 0.66 -> 13.8 rising through the c-ladder) — BUT the
c=0.10 solve formed genuinely PHYSICAL structure: gdev = 0.26 (strong-field scale) and a
DECAYING spectrum peaked at k=5 ([0.11, 0.19, 0.30, 0.15, 0.043, 0.022, 0.037] for
k=1..13). The richer basis can hold the solution; generic scipy trf+lsmr cannot converge
it. Seed representation also improved (|r| 7.64 -> 3.38 at the same seed).

**Next session's build (ranked):**
1. **Exact-Newton relaxation solver** (what Gundlach actually used): finite-difference
   Jacobian on the known block-tridiagonal sparsity, DIRECT sparse solve (spsolve/banded
   LU — no lsmr regularization), damped line search, both pins; coarse-Nz first.
2. Better seed: N=1600 evolver + re-bisected p* (cleaner SSH amplitudes, deeper resolvable
   echo) — the current SSH RMS(X+)=0.07 extraction is ~3-5x low.
3. Then: release amplitude pin -> release Delta -> Delta vs 3.4453 with the which-solution
   battery (nonzero decaying tails, no machine-zero, Delta vs Delta/2 impostor).

### Exact-Newton + conditioning probes (same session): the SEED is the common denominator

- Exact-Newton (direct sparse solves + backtracking): line search never accepts alpha > 0.03;
  |r| creeps 3.38 -> 1.85 and plateaus at every c — the Jacobian's predicted descent does
  not materialize. The a-clamp is EXONERATED (f in [0.32, 1.0] at the seed; no kink).
- Jacobian conditioning at the seed (Delta-frozen, pinned): cond ~ 1.3e6; the smallest
  singular directions (1.7e-3..1e-2) are MID-CYLINDER field rearrangements (z ~ -1.5..-3.2,
  mixed g/X+-; not xi0, not boundary modes) — the causally-weakly-pinned interior valley,
  characteristic of advection-dominated BVPs seeded far from the solution.
- ALL FIVE solver campaigns (LM-shoot, trf cliffs, trf-relax, pinned continuations,
  exact-Newton) now agree by measurement: the binding constraint is SEED QUALITY (the
  cylinder extraction from the N=800 run underestimates SSH amplitudes ~3-5x and carries
  frame/interpolation inconsistency). N=1600 re-bisection launched (needed by every future
  path). Secondary levers for session 2: pseudo-transient continuation (Psi-tc damped
  Newton), complex-step Jacobian, k~21 truncation (Gundlach 2N=128).

### N=1600 chain (final probe of session 1): the EXTRACTION pipeline binds, not the data

- The finer evolver is richer: 15 clean crossings (vs 8), t* = 15.192 +- 0.04, and the
  in-house **Delta_echo improves 3.216 -> 3.334** (N=800 -> N=1600; Gundlach 3.4453) —
  a resolution-convergent honest estimate now within 3.2%. This trend is itself a
  deliverable (the evolver-side redundant-recovery observable of the contract).
- BUT the relaxation seed built from it is WORSE (|r| = 13.7 vs 3.38; SSH RMS(X+) = 0.048
  vs 0.071): the cylinder-extraction pipeline (SSH locator with xi0-prime neglected,
  frame approximations, nearest-snapshot sampling) — not the raw evolution — limits seed
  quality. Newton descends further initially (3.40 -> 1.46) then stalls in the same
  plateau class (gdev frozen 0.084).

## SESSION-1 CLOSE — the wall, named exactly

**Delta and gamma are NOT measured (none faked).** The binding constraint, measured across
six solver campaigns and two evolver resolutions: **seed/extraction quality into a
weakly-conditioned interior valley** (cond(J) ~ 1e6; near-null mid-cylinder modes).
In-hand assets: the full evolver + p* at two resolutions (0.03751655962597 @ N=800,
0.03732817692976 @ N=1600), Delta_echo 3.33 (resolution-trending to 3.4453), the validated
relaxation system (vacuum EXACTLY 0, square, truncation-configurable), exact-Newton +
sparse-Jacobian tooling, and the two vacuum-watershed catches as house law.

**SESSION-2 RECIPE (ranked):**
1. **Extraction rebuild**: solve xi0(tau) from the FULL coordinate condition (with xi0'),
   locate the SSH per-tau on the fine grid, sample by 2-D interpolation in (tau, zeta)
   from MANY snapshots (not nearest), and fit boundary functions by least squares over a
   whole period band. Target: seed |r| < 1 with SSH RMS(X+) ~ 0.2-0.4.
2. **Psi-tc (pseudo-transient) damped Newton** on the relaxation system + complex-step
   Jacobian (kills FD noise in the flat valley).
3. **k ~ 21 truncation** ((64,20,21)) once descent works.
4. Then the pin-release ladder + which-solution battery (unchanged).

## SESSION 5 (2026-07-19 afternoon) — three catches, the first healthy state, endgame v1's honest failure

*(recorded in session 6; detail in git `def560c`/`c4ed641`/`e7e174b` + `docs/CONTINUATION_2026-07-19_session5_gamma.md`)*

- **WHICH-SOLUTION CATCH #3 — the g<0 junk valley:** the un-logged (64,20,21) grind
  fabricated a low-|r| valley with g.min = −0.457 (f ~ 3e5), exposed by Nz=60 evaluation
  (|r| = 18.2). Fix: the g-slot now carries **W = ln g** (positivity by parametrization;
  the g-equation becomes LINEAR, W,z = 1 − a²; the SSH condition goes log-clean).
- **IMPOSTOR #4 — the high-k linear wave:** the log-g chain descended to |r| = 8.8e-3 with
  ZERO low harmonics (k=1–5 ≈ 0, g ≈ 1) — a near-linear high-frequency wave satisfying the
  total-RMS pin (the pin doesn't select WHICH harmonics carry the amplitude).
  **Battery addition: LOW-K DOMINANCE.**
- **ROOT CAUSE of every diluted seed:** the extraction τ-window sat in the PRE-CRITICAL
  INFALL (t_G ∈ [1.2, 14.6] but echoes start ≈ 11.8 ⇒ ~70% transient). The clean window is
  one period inside the echo regime: T ∈ (t\*−t_first_crossing)·[e^{−Δ}, 1].
- **The center band is unmeasurable during the echo epoch** (the attractor is strong-field
  at ζ = −2: measured g(−2) = 0.27 vs asymptotic 0.94 — the D-021 resolution wall at
  extraction level) ⇒ **seed v3 drops Y₁ entirely**; the BVP's node-0 BCs own the center.
- **FIRST FULLY-HEALTHY STATE** (`v3_48.npy`): |r| = 0.0787 at (48,14,13),
  g ∈ [0.521, 1.064], LOW-K-DOMINANT decaying tail [0.0499 … 0.015], pin RMS 0.0320 —
  ALL battery checks green.
- **ENDGAME v1 VERDICT (`endgame.log`):** pad → (64,20,21) grind reached |r| = 6.9e-3 but
  DRAINED under the total-RMS pin — g → [0.978, 1.000], tail INVERTED to edge-piling:
  impostor #4 recurring at higher truncation (the grind emptied k<7 while holding the pin).
  Pin-release then converged to **THE VACUUM** (|r| = 1.5e-13 machine-zero, all-zero
  tails — the house tell), and the released Δ sat at its frozen initialization: **the
  "PRELIMINARY Delta = 3.445300" banner was the INPUT echoed back to 7 digits (zero
  Δ-gradient in the vacuum), NOT a measurement.** Δ still not measured, none faked.
- **Dying edits at the context wall** (committed session 6): `pin=('lowk', c, w)` — RMS of
  the k≤5 odd SSH X₊ COEFFICIENTS, immune to high-k relocation — and `lm(lam0=…)`.

## SESSION 6 (2026-07-19, post-rehydration) — endgame v2: battery-gated low-k ladders

- **Driver `endgame2.py`** (smoke-validated: vacuum exactly 0 at both truncations; the pad
  reproduces v1's padded |r| = 1.1385 EXACTLY against the lost inline pad; lowk pin row
  exactly 0 at the measured c = 0.017332; `v3_48` battery GREEN as-is). Ladder = sanity
  regrind → gated grind → pin-weight RAMP 30→10→3→1→0.3→0.1→0 (drain at any rung =
  PIN-SUPPORTED/off-manifold verdict) → Δ release (a Δ that does not MOVE from 3.4453 is
  flagged ZERO-INFORMATION — the suspiciously-exact law).
- **v2-64 MEASURED VERDICT (`endgame2.log`): the mid-flight gate fired at round 0.** With
  the lowk pin ON, the (64,20,21) grind took |r| 1.139 → 0.0403 while the battery went
  UNHEALTHY: the pin held (SSH k≤5 amplitude fixed; low-k interior only −40%; g still
  [0.611, 1.002]) but the FOUR newly-opened odd harmonics **k = 15–21 blew up in the
  interior** (0.085 at k=17 vs 0.030 at k=1 — lowk_dominant + tail_decays FAIL). Same
  low-|r|-but-wrong class as v1, caught immediately instead of at the end. ⇒ **K must be
  opened GRADUALLY (M-first-then-K), not four harmonics at once** — the optimizer answers
  fine-collocation residual rows with high-k junk when given the freedom.
- **run48 VERDICT (`run48.log`, 26 s): the pin releases CLEANLY at the native (48,14,13).**
  Every rung GREEN — g/tails byte-stable through the full ramp w: 30 → 0 and the state
  holds |r| = 0.0786 with NO pin at all ⇒ **the healthy state is NOT pin-supported; the
  vacuum does not swallow it.** BUT the released-Δ banner (3.4453001, drift +8.8e-8) is
  NOT a measurement: |r| sits on the KNOWN 7.9e-2 LM plateau (~3 orders above the
  truncation floor) — the optimizer cannot move in ANY direction including Δ, so Δ reports
  its initialization. **Battery hole found + patched: no absolute-|r| floor check** (it
  said GREEN on a non-solution). New check: `at_floor(|r| < 5e-3)` required on FINAL
  before any Δ claim. **The wall at 48³ is now cleanly isolated: the flat-valley plateau
  itself** (cond ~1e6, mid-cylinder near-null modes — the session-1 diagnosis), not the
  pin, not the vacuum, not high-k.
- **nzprobe VERDICT (`nzprobe.log`) — the diagnosis that unifies every failure this
  session:** the plateau is NOT Nz-independent. At Nz=60 the LM descends 7.86e-2 →
  2.65e-2 (still descending at maxiter) — **but by piling content against the KO=13
  edge** (tail 0.0156 → 0.1327 RISING; battery FAIL), the session-1 truncation-binding
  signature. Nz=80 stalls earlier (6.19e-2), battery marginal. Reading: **finer z
  resolves the sharp SSH-adjacent structure whose τ-content genuinely reaches k~21
  (Gundlach's own content) ⇒ either axis refined ALONE distorts — the solution needs the
  JOINT target (KO=21) × (Nz≥60)** — and v2-64 measured that opening that much freedom
  at once gets filled with junk. Also measured: the Nz=40 state carries node-scale
  z-roughness (raw upsampled |r| = 8.3 at Nz=60) — the valley's near-null mid-cylinder
  modes ARE node-scale wiggles.
- **runA-A0 VERDICT (`runA.log`): the lowk pin's OWN blind spot, measured.** With the
  annealed Tikhonov clamping high-k (ε=10) the padded-state grind STILL drained —
  g → [0.993, 1.000], everything ≈ 0 except **k=5 = 0.042: the optimizer relocated the
  entire pinned amplitude into the single k=5 coefficient** (an RMS over 6 coefficients
  is ONE scalar; the weak-field wave lives inside its level set). Battery caught it via
  strong_field. **Unifying fact across v1 / v2-64 / A0: from a FAR padded state
  (|r| ≈ 1.1) the descent direction is ALWAYS drain**, in whatever flavor the
  pin/penalty leave open — while from the near-manifold 48³ state it never is (run48).
  Corollary: the pad-jump 0.079 → 1.14 says the 48³ "solution" is flattered by coarse
  τ-collocation (aliasing-hidden inter-collocation error).
- **Pin evolution: `pin=('vec', c6, w)`** — the 6 SSH X₊ coefficients clamped
  INDIVIDUALLY (not relocatable). Bug-catalog addition: **RMS-form pins (total AND
  lowk) are relocatable; only a vector pin is not.**
- **runMK design: the converged-rung ladder** — never be far: rung 0 = the pure M-jump
  48→64 at (14,13) (same basis, same u, finer collocation — directly measures the
  aliasing-hidden error and whether the old basis can close it); then K opens TWO
  harmonics per rung from a CONVERGED state (annealed Tikhonov on just the newest;
  vector pin at w=30 throughout); then Nz climbs 40→50→60; then ramp + Δ release with
  `at_floor`.
- **runMK rung-0 MEASUREMENTS (`runMK.log`):** (a) **the VECTOR pin defeated the drain**
  — first far-state grind (raw |r| = 0.83) that kept g strong-field ([0.505, 1.101])
  instead of collapsing: the drain pathology of v1/v2-64/A0 is CLOSED by clamping the 6
  SSH coefficients individually. (b) The M-jump quantifies the aliasing-hidden error:
  the 48-collocated |r| = 0.079 becomes 0.83 raw / **0.117 re-ground at M=64** — the
  KO=13 basis floors at ~0.1 under honest collocation, with the missing k>13 content
  aliased into its edge harmonics (tail flat/edge-elevated). (c) **Gate-design lesson:
  tail-shape checks are INVALID mid-ladder** — a KO-truncated basis at honest
  collocation is EXPECTED edge-loaded until K opens; the first runMK aborted itself on
  exactly this. Mid-rungs now gate on catastrophe only (drain/vacuum/Δ-band, warns
  printed for the rest); tail shape + `at_floor` gate the FINAL. runMK v2 relaunched
  (`runMK2.log`).
- N=3200 p\* bisection hedge in flight (`pstar3200.log`) → third point for the in-house
  Δ_echo resolution series [3.216, 3.334, …] → 3.4453.
