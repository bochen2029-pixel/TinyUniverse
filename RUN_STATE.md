# RUN_STATE.md

**As of:** 2026-07-11 · **Milestone:** M1 `canvas` · **State:** M0 CLOSED — nexus v1.0.0 green (11/11), golden `ad64f810` frozen, 24.1 s runtime, determinism verified in-process ×2 and out-of-process ×3.

## Current task

M1 `canvas` — first light. Order of work:
1. Resolve Q-001 to its default (god-hand v1) and freeze `contracts/frame.contract.md` v1.0.0 (regime-mask canonical table from nexus N10: QUANTUM 0x1 · CLASSICAL 0x2 · REL 0x4 · COMPACT 0x8 · BOUND 0x10 · MASSLESS 0x20).
2. liborrery lift into `core/lib/` (verbatim, pinned ORRERY commit recorded; swap nexus's self-contained BLAKE2b for the envelope — golden supersession expected, signed in `goldens/nexus/NOTE.md`).
3. Vulkan + CUDA external-memory interop skeleton (cuda-samples `simpleVulkan` / Mímir as reference); two-stream frame loop.
4. CINEMATIC stack port from GARGANTUA (HDR → mip bloom → auto-exposure → AgX + ACES parity mode → dither), ImGui HUD shell, damped camera.

**M1 gate:** 1M inert drifting particles at 60 Hz 1080p with every CINEMATIC §7 box checked — the "prettiest N-body screensaver". If that frame doesn't beat three.js, stop and fix before physics.

## Chores carried

- clang++/g++ parity build of nexus (owed from M0; no non-MSVC toolchain installed).
- Pin Vulkan SDK version in BUILD.md when installed.

## Standing context

M0 evidence: `runs/nexus_v1.0.0_freeze.json` · `goldens/nexus/` · `nexus/MODULE.md`. Notable frozen physics: 1PN precession 2.3% off analytic (higher-order PN, inside 5% gate); P–W ISCO to 2.7e-9; Ratchet MC within 0.45% of closed form; the γ=10⁸ β-representability ceiling demonstrated deterministically (D-011). Repo docs are authoritative over agent memories where they disagree.
