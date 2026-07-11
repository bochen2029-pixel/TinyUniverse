# CLAUDE.md — TINY UNIVERSE build harness

You are building **TINY UNIVERSE**, a scale-compressed physics universe: a native Windows C/C++/CUDA app where quantum → Newtonian → relativistic → black-hole physics emerge from particle count and density under retuned constants (ħ↑, c↓, G↑). Read this first, every session. `ARCHITECTURE.md` is the spec; this file is the operational entry point.

## Bootstrap (every session start)

1. Read `ARCHITECTURE.md` — dials (§6), invariants (§5), frame contract (§7), milestones (§11).
2. Read `RUN_STATE.md` — current task + next concrete action.
3. Read `TASKLIST.md` and `DECISIONS.md` (open/active).
4. `git status` + `git log --oneline -15` — ground truth. Once the harness exists, run it cold (`harness/`) to verify reality matches claimed state before proceeding.

Time-to-productivity target: under 2 minutes.

## The doctrine (inherited from ORRERY, adapted for an interactive app)

- **The spec is the product. The contract is sacred. The golden is load-bearing. The code is ephemeral.** You are never working on "the app"; you are working on one module against its fixed contract.
- **Contract-first.** No module source until its `contracts/<name>.contract.md` (+ schema where applicable) exists and is operator-reviewed.
- **Determinism or it doesn't ship** — TU form: fixed timestep; same dials + seed + **input trace** + step count ⇒ byte-identical declared state hash. Visuals, wall-clock, and frame pacing are explicitly **non-declared**. All randomness via counter-based RNG (liborrery `rng.cuh` idiom): value = f(seed, particle_id, tick) — never stateful.
- **Every physics module names its oracle** — the analytic/CPU-double reference its correctness is gated against (`tiny_nexus` is oracle-grade). Oracle gate failure ⇒ run-state SUSPECT, never a silent fallback.
- **Regimes emerge; nothing hard-switches.** Representation handoffs (wavefunction ↔ point particle) are allowed but must be physics-triggered (Ratchet redundancy threshold), declared, and deterministic.
- **Beauty is law.** Any user-facing visual must pass the `CINEMATIC.md` checklist before it ships. "Looks fine" is not a gate; the checklist is.
- **Discipline tightens with capability.** Golden-gate everything; two-pass any claim the project will cite publicly. RAYFORMER (`C:\RAYFORMER`, ADR-007) is why: every isomorphism/speedup claim gets an honest baseline and a measurement before it is believed.
- Never wait, always log; durable over fast; future-me is a stranger; measure resource size before consuming (`C:\chunker\` for big files, `C:\everything\search.py` to locate).

## The build loop (per module)

1. **Contract** — write/review `contracts/<name>.contract.md` (+ `.schema.json` for headless faces). Semver from 1.0.0.
2. **MODULE.md** — purpose, contract link, invariants, oracle, build command, known issues.
3. **Implement** — C++/CUDA default (single-file where possible); headless face always present.
4. **Build** — per `BUILD.md` (vcvars64 + nvcc `-arch=sm_89`; CMake preset for multi-file).
5. **Golden** — freeze `(dials, seed, trace, steps) → declared-state hash` in `goldens/<name>/`.
6. **Verify** — harness green (compile + selftest + golden). Two-pass for citable claims.
7. **Register** — update `ARCHITECTURE.md` §11 row + `TASKLIST.md`. Spec never lags code.
8. **Save point** — atomic commit (module + contract + golden + MODULE.md + RUN_STATE + TASKLIST).

## Build order (decided; gates in TASKLIST.md)

M0 `nexus` → M1 `canvas` (Vulkan⇄CUDA interop + CINEMATIC port) → M2 `newton` → M3 `arrow` (thermo + Ratchet) → M4 `einstein` → M5 `gargantua` (BH) → M6 `planck` (quantum bubbles) → M7 `cosmos` (torus, little-planet, light history). RTXDI / DLSS / Slang / UE-TextureShare shell: backlog, measured adoption only.

## Hard rules (never violate)

- Never write module source before its reviewed contract.
- Never ship non-determinism in the declared path; never let visuals leak into a golden.
- Never make a breaking contract change without a MAJOR bump + migration note.
- Never claim a speedup/isomorphism without an honest baseline (ADR-007 protocol).
- Never use `--use_fast_math` in declared paths; no float atomics in declared reductions (liborrery shapes only).
- Never edit `C:\ORRERY` from this repo. liborrery files are lifted **verbatim** (D-005); any divergence is a fork, named and justified in DECISIONS.md.
- Exit codes: 0 pass · 1 declared-gate-fired (real negative result) · 2 error. Never conflate 1 and 2.

## Capabilities

- **CUDA 13.1 + RTX 4070 Ti SUPER 16 GB** (`-arch=sm_89`), MSVC 2022 via vcvars64 — see `BUILD.md`.
- Quarries: `C:\blackhole\files\blackhole.cu` · `C:\ASTRA-7\proto\astra_nexus.cpp` · `C:\Buddhabrot_CUDA` · `C:\ORRERY\lib` · `C:\RAYFORMER` (lesson) · `C:\buddhabrot-main` (CMake/ckpt patterns).
- Web search for primary sources; cite precisely.

## Definition of done (never "finished", always compounding)

A windowed universe that is beautiful under the CINEMATIC gate and honest under the physics gates — every regime claim ("this is real quantum mechanics") backed by a green oracle test a stranger can run cold from `BUILD.md` in one pass.

*Build one module right. Freeze its golden. Let the universe compound.*
