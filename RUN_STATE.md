# RUN_STATE.md

**As of:** 2026-07-11 · **Milestone:** M3 `arrow` · **State:** M2 CLOSED — the universe gravitates. 1M particles with live PM gravity: **347 fps avg / 159 min** @1080p. Four scenario goldens frozen and reproduced (`GOLDEN OK` ×4): kepler `448847ec` · threebody `f67abce4` · cloud `975389db` · galaxy `ff17431a`.

## What M2 established

- **Deterministic self-gravity**: 128³ PM (fixed-point CIC deposit per Invariant 4 → cuFFT Poisson → FD force grid), KDK + Kahan drift, fp64 tiny solver for precision-reference scenarios (D-014). Periodic box = the 3-torus, natively — M7's topology already lives in the gravity solver.
- **Physics receipts**: kepler energy drift < 1e-6 over 10⁶ ticks (fp32-boundary/fp64-internal); the figure-8 three-body choreography bounded over 29.5 periods; cloud relaxes at 0.46% drift; galaxy 1.9% @10⁴ ticks with momentum conserved to 2e-6. Regime masks live: the galaxy disk runs 97.5% REL (v ~ 0.16c — the compressed dials at work, M4's opening argument).
- **The envelope face is real**: `--scenario X --json|--golden` emits liborrery `full_envelope` and uses `golden_check` — the app's headless face is now contract-grade (D-003 fulfilled end-to-end).
- Evidence: `runs/newton_galaxy_25s.png` (differentially sheared arms — dynamics, not decoration), `runs/newton_cloud_collapse.png`, `contracts/newton.contract.md`.

## Current task — M3 `arrow` (thermodynamics + inscription)

1. Contract first: T/P/S meters (dS/dτ display only — the Mirror's ledger rule); mixing scenarios that make the arrow of time visible.
2. **The Ratchet mechanic**: record events into environment DOF, per-bubble redundancy R, collapse past `(1−p)ρ = p` — mapping contract anchored to ORRERY `ratchet` receipts (nexus N6 verified the closed form in-repo).
3. Short-range forces (spatial hash + P³M correction) if the collision/thermal physics needs them — deferred from M2.
4. harness/verify.py (carried from M0/M2 — overdue; nexus + 4 scenario goldens make it worth automating now).

**M3 gate:** entropy meter shows dS/dτ ≥ 0 emerge from reversible microdynamics in a mixing scenario; the Ratchet in-sim matches nexus N6 at a small config; a detector scenario visibly writes records before collapse (prep for M6's double slit).

## Chores carried

P1 Vulkan presentation (SDK) · ImGui at P1 · TAA · clang/g++ nexus parity · art pass (dust lanes; galaxy de margin — D-014 note) · cufft64 dll in dist packaging (D-014) · harness/verify.py.

## Standing context

Oracle: tiny_nexus v1.0.0 (`ad64f810`). Dials v0 + frame contract v1.0.0 frozen. Build now links `core\lib\envelope.cpp` + `cufft.lib` (BUILD.md golden path needs this update — do at M3 contract commit). Repo docs authoritative over agent memories.
