# TINY UNIVERSE — continuation prompt (handoff from the 2026-07-12 session)

*Paste this into a fresh session. It tells you what to read first, where things stand, the recommendations, and the bigger-picture roadmap.*

You are resuming **TINY UNIVERSE** (`C:\TinyUniverse`) — a scale-compressed physics universe: a native Windows C/C++/CUDA app where quantum → Newtonian → relativistic → black-hole → cosmological physics **emerge from one substrate** under retuned constants (ħ↑, c↓, G↑), built under the ORRERY doctrine (contract-first, golden-gated, deterministic-or-it-doesn't-ship, every claim gated against an oracle or honestly walled). It is the *inhabitable* sibling of ORRERY.

## 0. Before you reason — SOAK context (read in this order, then VERIFY vs git)

1. `C:\Users\user\.claude\projects\C--TinyUniverse\memory\MEMORY.md` → the linked `tiny-universe-rehydration.md` (reconstitution pointer) + `local-disk-search-tools.md` (tooling, §1).
2. In `C:\TinyUniverse`: **`RUN_STATE.md`** (current position + the OPERATOR WORK QUEUE) → **`ROADMAP.md`** (the three-axis north-star plan — the bigger picture) → **`DECISIONS.md`** (read D-020…D-027 closely — the *why*) → `ARCHITECTURE.md` §5/§6/§11 → `CLAUDE.md` (doctrine + bootstrap ritual).
3. **`git log --oneline -15`** = GROUND TRUTH. The docs have lagged git before — if RUN_STATE/TASKLIST disagree with git, **git wins; reconcile the docs** (spec never lags code). HEAD should be `5648168` (or later), tree clean (untracked `runs/*.png|.gif` are non-declared render artifacts).
4. Do NOT re-run the full golden harness just to "check" on rehydration — receipts are logged + two-passed; trust them. (`python harness\verify.py` for a genuine cold verify; GPU is shared, preflight exits 3 if <2 GB free; the ~15 CPU goldens run regardless.)

**Disk beats memory. Never resume blind from a summary.**

## 1. Tooling (the operator flagged this — do NOT use grep / Get-ChildItem -Recurse)

- files by name/date/size → `python C:\everything\search.py "QUERY"` (`ext: path: dm:today size: |=OR !=NOT "phrase"`; `-l --recent --in DIR`). Forward-slash the script path in bash to dodge `\e` mangling.
- content inside files → `C:\everywhere\build\Release\everywhere.exe [-i -n -F -e PAT ...] [-g GLOB] PATH` (the ripgrep replacement).

## 2. Where things stand (2026-07-12, HEAD 5648168)

**v1 (M0–M7): COMPLETE** — the regime ladder + the world (21 goldens): beauty → gravity → arrow → relativity → black holes → quantum → cosmos (torus, little-planet, a(t), roaming bubbles).

**v2 substrate ladder** (the intellectual core — *one field, one clock, one lattice*):
- **N0 `substrate_nexus`** (CPU fp64 spherical EKG oracle) — Choptuik Type-II *transition* demonstrated; precise γ DEFERRED to AMR (D-021, honest wall).
- **N1 `field`** (GPU Schrödinger–Poisson) — the PM+ψ gravity WELD, 6/6 goldens.
- **N2 `lapse`** — the clock: gravitational redshift, exact Schwarzschild (2/2).
- **N3 `curve`** — geometry: light bending (the 1919 factor-of-2 = the N2-time half + the N3-space half), perihelion precession, Shapiro delay (3/3). **→ all FOUR classical tests of GR now pass on the substrate.**

**v1 polish + crowns (this session):**
- **2.5PN inspiral** (`nexus/inspiral_nexus`) — the Peters GW chirp + circularization, machine precision.
- **Q-006 RESOLVED** (`nexus/precession_nexus`) — the SR+1PN **7π superposition is CORRECT**; the app's old 6.41π was a **normalization artifact** (nominal vs the orbit's force-distorted actual semi-latus rectum). *A RAYFORMER save: the conclusion FLIPPED only when the actual `p` was measured — the first draft had it backwards.*
- **galaxy little-planet orbit** — the "see the universe" beauty win (`runs/galaxy_orbit.gif`). The a(t) `expand` timelapse is an honest CINEMATIC wall (dim, near-uniform box).
- **fluid-β Stage-A BANKED** (`substrate/fluidcss_nexus`) — Evans–Coleman CSS radiation-fluid background: oi*=3/8, sonic (3/2,2/√3,3/4,−1/√3) EXACT, 2m/r=1/3 — all emergent. Golden `fluidcss_stageA`.

**The honest walls (do NOT grind or fake — each needs a specific unblock):**
- **fluid-β Stage-B → β=0.35580192**: the perturbation ∂ₛ-coupling fails its gauge-mode exactness gate. Resume = re-transcribe HKA (5.5)–(5.10) from the primary TeX **or** re-derive the coupling, then the armed eigenvalue solver reads the non-gauge zero (discard the κ≈0.357/1 gauge modes). See `tournament/gamma/phase4/RESULTS_hka_beta.md`.
- **scalar-γ = 0.374** (DSS crown, G97 gr-qc/9604019): not started; harder (τ-periodic hyperbolic BVP + Floquet). The HKA regular-center template that landed Stage-A is the playbook.
- **precise Choptuik γ** (N0 / N4 horizon): a uniform grid caps the self-similar curvature (D-021); needs **AMR** (Berger–Oliger; a concrete recipe is in `tournament/gamma/phase1/amr-berger-oliger.md`).

## 3. The discipline (non-negotiable) + hard-won lessons

- **Contract-first** (no source before `contracts/<name>.contract.md`), **golden-gated** (`(dials,seed,steps)→blake2b`, sm_89-pinned for GPU), **deterministic-or-it-doesn't-ship**, exit 0/1/2/(3=GPU-contended) never conflated.
- **Every module names its oracle**; gate against exact/analytic truth. **Honesty over progress**: a real wall is reported with its exact resume state, never tuned away (D-016/D-021). Two-pass any citable claim.
- **The RAYFORMER lesson (D-027, D-022):** measure/normalize against what the object *actually is*, not its nominal label — it flipped Q-006 and caught the N1 Poisson 1/N³ bug. **Measure before concluding.**
- **Don't grind a wall** — log it honestly and pivot. **Don't re-verify recorded goldens** on rehydration. **Don't poll your own background tasks.** Reorient minimally.
- Each `_nexus` tool is a self-contained single file (blake2b + envelope + golden). CPU oracles build with `cl`; GPU with `nvcc -arch=sm_89` (no `--use_fast_math` in declared paths; links cufft). Registration per milestone: ARCHITECTURE §11 + TASKLIST + DECISIONS + RUN_STATE + MODULE, then an atomic commit.

## 4. Recommendations — what's next (ranked)

The physics-oracle track is in excellent shape; the clean-landable oracle wins are largely done. Highest-leverage directions now:

1. **Axis B — the RENDERER (the biggest single gap to the north star).** The universe is honest but you can't yet *fly through* it. Start **R0 `interop`** (Vulkan⇄CUDA shared buffer → live window) per ROADMAP §4. **Prerequisite: the Vulkan SDK is NOT installed (measured, D-002/Q-005)** — the operator must install it first, or you unpark the existing GL-blit / `--shot` path. This turns the instrument into an *inhabitable* universe (the whole thesis). Highest strategic value.
2. **Axis C — the AMR contract (unlocks the crowns).** A well-scoped build: Berger–Oliger AMR wrapping the existing EKG solver (recipe in `amr-berger-oliger.md`; zero-refinement must reproduce the N0 golden bit-for-bit). Unlocks the precise **Choptuik γ** (N0) AND **N4 `horizon`**. The most valuable *physics* next step.
3. **N4 `horizon` (transition-only).** Emergent trapped surfaces on the GPU lattice (reproduce N0's Type-II transition, BH mass→0 at threshold), deferring the exponent to AMR. Ladder progress without the AMR wall.
4. **fluid-β Stage-B** — closest to "armed," but only with a focused research session (re-transcribe HKA 5.5–5.10). Uncertain payoff; the overnight run couldn't crack it.
5. **Kerr art pass** — GARGANTUA per-pixel geodesic BH render (needs Axis B or the lift from `C:\blackhole\files\blackhole.cu`). Beauty, subjective gate.

## 5. The bigger-picture roadmap (from ROADMAP.md — the north star)

**Three axes converging into one binary:** **A** substrate depth (physics, N0→N7), **B** beauty (renderer, R0→R5), **C** rigor (the crowns β/γ + RAYFORMER discipline). *The thesis:* a black hole in the final demo is (A) an emergent horizon, (B) a path-traced Kerr shadow, and (C) a configuration whose formation threshold reproduces the published critical exponent — the same object passing all three gates at once.

**Ambition tiers:** 🥉 physics real & unified (N0–N7 headless) · 🥈 you can see it (R0–R2, CINEMATIC-gated) · 🥇 playable + one crown landed (RTX 60 fps + β or γ) · 🏆 the north star (both crowns + the "Collapse" moonshot + a publishable RAYFORMER isomorphism/speedup).

**Where we are on the arc:** **Axis A is strong through N3** (full weak-field metric + all 4 GR tests + the relativistic-binary physics + the fluid crown's Stage-A). **Axis C has Stage-A banked**; the crowns (β/γ) are walled pending AMR / primary-source re-derivation. **Axis B is the least-advanced and the biggest gap** — the substrate is honest but not yet inhabitable.

**Strategic recommendation: pick a lane and go deep — make it *seeable* (Axis B renderer; needs the Vulkan SDK) or *provable* (Axis C AMR → the Choptuik crown).** The oracle ladder (Axis A) can keep climbing (N4–N7) but the top rungs are research-grade; each should be one focused contract-first module.

## 6. Cadence

The operator works directive-by-directive ("go N4", "proceed", "keep going"), expects a full milestone per directive then a pause, wants to *see* visual results (use `SendUserFile`), and quota can cut mid-flight — so commit early, keep RUN_STATE current, atomic-commit every milestone close. Suspend and ask plainly if the operator returns confused or frustrated.

*Build one rung right. Freeze its golden. Prove its number. Make it beautiful. Let the universe compound.*
