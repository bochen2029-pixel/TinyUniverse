# TINY UNIVERSE — SESSION-5 MEMORY DUMP (γ-crown campaign, 2026-07-19 ~13:00)

*Perpetual-memory dump (operator-requested at 87% context, precaution against compaction).
Self-distilled by the session that lived it. If you are a fresh/compacted session: this file
IS your rehydration; follow §0 exactly. Trust files over prose; verify against git.*

---

## 0 · RECONSTITUTION POINTER (read in this order, then verify)

1. **This file, whole.**
2. `C:\Users\user\.claude\projects\C--TinyUniverse\memory\MEMORY.md` (the index hook holds
   the compressed running history of ALL sessions incl. today's β-crown close).
3. `C:\TinyUniverse\tournament\gamma\dss\RESULTS_dss.md` — the γ-thread canonical results
   (every measured verdict, session-2 recipe).
4. `RUN_STATE.md` (top sections) → `DECISIONS.md` D-032 → `contracts/similarity_nexus.contract.md`.
5. **Verify:** `git -C C:\TinyUniverse log --oneline -8` (HEAD ≥ `6858280`, tree clean or
   only dss checkpoints); the γ work lives in `tournament/gamma/dss/` (Python; CPU-only;
   run from that directory). Full harness NOT needed (40/40 green this morning; don't re-run).
6. The raw transcript backstop: this session's `.jsonl` under
   `C:\Users\user\.claude\projects\C--TinyUniverse\` — grep it for specifics, never re-read whole.

---

## 1 · CORE (never drop)

- **GOAL:** the scalar-field DSS crown — measure **Δ = 3.4453 ± 0.0005** (echoing period,
  Stage-A background) and **γ = 0.374 ± 0.001** (= −1/λ₁, Stage-B Floquet perturbations) on
  the Choptuik critical solution, per **operator-APPROVED** `contracts/similarity_nexus.contract.md`
  v0.1.0 (Gundlach gr-qc/9604019 = `tournament/gamma/GUN96_src/chop2.tex`, gitignored).
  Python research first, then the C++ `substrate/similarity_nexus.cpp` port (the proven fluid pipeline).
- **POSITION (the furthest point):** the global-relaxation solve at truncation (M,KE,KO) =
  **(48,14,13)**, amplitude-pinned (c = 0.04 on SSH odd-RMS(X₊)) and Δ-FROZEN at 3.4453, has
  reached **|r| = 0.135** (LM, checkpoint `tournament/gamma/dss/lm_v1.npy`), gdev = 0.153
  (strong-field), spectrum peaked at k = 9 (0.228) — still pressing the truncation edge.
  Trajectory today: seed 782 → 40 → 4.3 (extraction fixes) → 0.93 (hand-rolled-Jacobian
  Newton) → **0.135 (LM)**. Descending, not stalled.
- **THE SINGLE NEXT ACTION:** continue LM from `lm_v1.npy` with a bigger budget (it exited
  'slow', λ ~ 2e-2 — not converged, grinding); THEN escalate truncation to **(64,20,21)**
  (`R.configure(64,20,21)`; re-extract the seed OR upsample lm_v1 by zero-padding harmonic
  slots per node — upsampler must map coefficient layouts NE15/NO14 → NE21/NO22); grind to
  the truncation floor; THEN the release ladder: release the amplitude pin (Δ still frozen)
  → release Δ → **measure Δ** with the which-solution battery.
- **HARD CONSTRAINTS:** never fake a number (D-016/D-021/D-032); **the machine-zero tell**
  (|r| ≈ 1e-15 at finite truncation = the VACUUM — check spectral tails EVERY convergence);
  the Δ/2 = 1.72 impostor check; contract gates unchanged; commit early and often; push after
  commits (repo public, HKA99_src/GUN96_src must stay gitignored; local branch
  `backup/pre-excise-2026-07-19` must NEVER be pushed).
- **SOURCE OF TRUTH:** everything committed through HEAD `6858280`+; solver checkpoints
  (`*.npy`, gitignored? — NO, .npy not ignored; they are small and committed? CHECK: they are
  untracked scratch — treat `lm_v1.npy` as the live checkpoint, regenerate via the pipeline
  if lost: seed_v2.npy → newton → lm; pipeline is ~15 min total).

## 2 · RING 1 — the active thread, complete technical state

### Conventions (hardest-won; a wrong sign here cost hours)
- Gundlach time: **τ = ln(t\* − t_G)**, DECREASING toward the singularity (τ→−∞); t_G =
  central proper time from the evolver (αa=1 outer-normalized evolution, t_G = ∫α(0)dt).
- **ζ_G = ζ_raw − ξ₀(τ)** with ζ_raw = ln(r/T), T = t\*−t_G. SSH at ζ_G = 0 (the LAST
  relaxation node, z=0 exactly; the singular point is never evaluated — midpoint scheme).
- Parity-in-FREQUENCY (κ=0 symmetry): a, g, ξ₀ carry EVEN τ-harmonics; X± = X±Y carry ODD.
  Raw evolver data has even X±-content = transient junk — parity projection is CORRECT.
- τ-gauge: ξ₀'s k=2 SINE = 0 (all fields rotated by the common phase ph2 before fitting).
- Fields: X = √(2π)(r/a)Φ, Y = √(2π)(r/a)Π (the lapse CANCELS); g = a·α(0)/α (central renorm).
- Center asymptotics (Gundlach app A): Y = Y₁e^{ζ} + Y₃e^{3ζ}, X = X₂e^{2ζ},
  g = 1 − (Y₁e^ζ)²/3, X₂ = (e^{ξ₀}/3)(Y₁′ − (1+ξ₀′)Y₁); bounded modes ⇒ node-0 BC relations.
- SSH conditions: D₀ = (1+ξ₀′)e^{ξ₀}g = 1 (even) and C₋ − e^{ξ₀}gX₋,τ = 0 (odd) at z=0.

### File inventory (`tournament/gamma/dss/`)
- `nr_evolve.py` — near-critical evolver + echo analysis + v1 seed pipeline.
  **p\* = 0.03751655962597 (N=800) / 0.03732817692976 (N=1600)** (bisections to ~4e-14);
  t\* = 15.1917 ± 0.04 (N=1600); **in-house Δ_echo: 3.216 (N=800) → 3.334 (N=1600)** —
  resolution-convergent toward 3.4453 (a real deliverable). max|φ(0)| ≈ 0.51.
- `nr_extract.py` — **extraction v2** (the seed maker): 2-D (t_G, r) interpolation (Cyl2D);
  ξ₀ from the FULL coordinate condition (1+ξ₀′)e^{ξ₀}g = 1 iterated (stabilized: clamped
  interpolation-only roots in window ζ_raw∈[−0.5,0.9], smooth_even(kmax=6), 0.5-damping);
  Y₁ per-phase band-LSQ (ζ∈[−4.5,−2], r ≥ rmin = 4dr masked) + harmonic reconstruction over
  resolvable phases (~13% of phases have NO center data at any period — hole filled by
  periodicity); `build_seed`: EVERY node row = resolvability-weighted harmonic LSQ over
  phases with per-PHASE synth/data blending (deep center synthesized from Y₁-asymptotics).
  RMS(Y₁) ≈ 1.79 (consensus band); odd-RMS(X₊(SSH)) ≈ 0.04.
- `nr_relax.py` — the global relaxation system (Gundlach app D): harmonic-coefficient
  unknowns per ζ-node; `configure(M,KE,KO)`; square system NPF·Nz + NE; **vacuum control
  EXACTLY 0** at (32,10,9) and (48,14,13); pin=(c,w) appends w·(RMS_odd(X₊[last]) − c).
- `newton_drive.py` — `jac_fd` = **HAND-ROLLED sparse FD Jacobian** (explicit NE dense
  columns + own 3-coloring over nodes). ⚠ **scipy's grouped `approx_derivative` silently
  returns EXACT-ZERO columns for near-zero dense entries** (measured: manual dR/du = 1.8e3
  where scipy gave 0) — NEVER go back to it. `newton` (line search; plateaus in valleys),
  **`lm` (adaptive-λ, direct spsolve on JᵀJ+λD; THE working solver — 0.87 → 0.135)**.
- `nr_dss.py` — the retired two-sided shoot (structural wall: off-manifold the a-constraint
  has no real solution ⇒ cliffs; SSH X₋ DC = parity symmetry-zero). Keep for the record.
- `dss_drive.py` — trf-continuation driver (superseded by newton_drive).
- `RESULTS_dss.md` — canonical results/verdicts. Logs: `newton_v3.log`, `lm_v1.log`, etc.
- Checkpoints (untracked .npy, regenerable): `seed_v2.npy` (extraction v2 seed, |r|=4.3),
  `newton_v2/v3.npy`, **`lm_v1.npy` (BEST: |r|=0.135)**.

### The bug catalog (each measured; do not re-hit)
1. φ(0,t) must be EVOLVED (∂_tφ|₀ = (α/a)Π|₀) — ∫Φ pins φ(0) = 0, erases echoes.
2. Spurious +wΠ/2r in the Φ-equation (plain gradient!) = origin blowup; Π-eq in
   ∂_r(wΦ)+2wΦ/r form with parity ghosts.
3. Polar-slicing freeze: detect collapse at 2m/r > 0.9 ∨ α(0) < 0.02 (0.995 unreachable).
4. Vacuum watershed ×2: unpinned relaxation → flat space at |r| = 2e-28 (the TELL);
   low-pin → weak-field perturbative pseudo-solution (r ~ amplitude²).
5. Δ weakly determined at low amplitude — FREEZE Δ during amplitude work.
6. ξ₀ fixed-point without smoothing/clamps diverges (ξ₀ → +29).
7. Sub-grid sampling (r < 4dr) poisons fits (Y₁ 3× inflated); mask + synth-fill.
8. Per-node binary synth/data switches put seam spikes in the interior residual.
9. scipy grouped-FD dead-column artifact (above).
10. The evolver's outer BC reflects — dispersal windows must beat the reflection return.

### Which-solution battery (MANDATORY before believing any converged Δ)
Nonzero DECAYING spectral tails (not machine-zero, not edge-piling) · gdev = O(0.1) ·
Δ ∈ [3.40, 3.49] and NOT ≈ 1.72 (the Δ/2 impostor) · release-order: pin off (Δ frozen) →
Δ free · convergence: Δ stationary under truncation (48,14,13)→(64,20,21) and Nz 40→60.

## 3 · RING 2 — earlier today (terse; all verifiable in git/docs)

- **β crown CLOSED** (D-032): the 3-session wall was the BACKGROUND — Friedmann-as-EC
  (V₀=−c_s IS the Friedmann point of the §IV sonic line; V≡−√(1−1/A) to 1.9e-10). True EC:
  V₀\*=0.112439401388. β measured 3 ways; C++ `fluidcss_nexus` v1.0.0, goldens `27af7920` +
  `9f8587fd` (κ₀=2.810577 → β=0.3557988 vs lit 0.35580192). Harness = 40 goldens ALL GREEN.
- QC + public-push day: history excised of HKA99_src per operator (hash map
  `docs/HASHMAP_2026-07-19.md`); `backup/pre-excise-2026-07-19` LOCAL-ONLY.
- similarity_nexus contract v0.1.0 drafted → operator-APPROVED ("proceed with the gamma crown").

## 4 · RING 3 — roadmap beyond the γ crown (operator queue)

1. γ Stage-B after Δ: Floquet perturbations (λ₁ = −2.674 → γ) on the converged background —
   linear problem, same machinery + shifted operator (A+λB); gauge control at λ=0 exactly.
   Then C++ `similarity_nexus` port + goldens + contract v1.0.0 (mirror the fluid close).
2. Renderer Axis-B (R0 Vulkan⇄CUDA) — blocked on Vulkan SDK install (measured gate, D-002).
3. AMR contract → N4-GPU `horizon` (D-028 reinstatement; recipe in tournament/gamma/phase1/).
4. v1 polish tail: GARGANTUA Kerr art pass.

## 5 · Cadence + operator notes

Operator: Bo Chen, terse green-light directives ("go", "keep at it", "keep going until it
lands or hard wall"), values honest walls over fake numbers, quota/machine-sleep can cut
mid-flight → commit+push every increment, keep RUN_STATE/RESULTS current. Machine slept
once today (~06:50–09:40); background tasks survive but wakeups drift. Context was at 87%
when this dump was written — if you are reading this fresh, assume compaction happened.

## 6 · VERBATIM TAIL (most recent last)

> **Operator:** "contract approved, proceed with the gamma crown etc etc (or whatever else
> you think is best/ your best recommendation etc ) go!"
> **Operator:** "okay... keep at it, thanks!!"
> **Operator:** "keep going until the gamma crown lands or hard wall etc"
> **Operator (this trigger, at 87% context):** "also, as a precaution whenever you get a
> chance, start writing to disk as file your current memory just in case of compaction so
> that i can manually rehydrate you or a new session in worst case if and as needed,
> basically starting logging things to file, your memory, things of this session and your
> recommendation paths and roadmaps forward etc etc, somewhere obvious and that is written
> to disk etc thanks."
> **Me (state at dump time):** LM broke the Newton plateau — |r| = 0.135 at (48,14,13),
> pinned+Δ-frozen, checkpoint `lm_v1.npy`; next = longer LM → (64,20,21) → release ladder
> → Δ measurement with the battery. Δ/γ NOT measured yet, none faked.

*Continue the grind. Check the tails. Never trust a perfect zero.*
