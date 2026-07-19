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
