# RUN_STATE.md

**As of:** 2026-07-11 · **Milestone:** M2 `newton` · **State:** M1 CLOSED — the beauty gate is settled. 1M particles at 1080p: **499 fps avg / 226 min** (SSAA 2×: 178/152); CINEMATIC §7 **10/10**. Evidence: `app/MODULE.md`, `runs/firstlight_cinematic.png` / `_physical.png` / `_ssaa.png`.

## What M1 established

- `app/tinyuniverse.cu` v0.1.0: all rendering in CUDA (HDR → mip bloom → auto-exposure → AgX/ACES → astro stretch → dither → CUDA HUD), two-stream frame loop per frame.contract v1.0.0, damped camera, cinematic/physical honesty toggle, one binary two faces (`--shot` headless / windowed / `--bench`).
- liborrery lifted verbatim (ORRERY d56c4c7; `core/lib/LIFT.md`); BLAKE2b byte-compat green → nexus golden `ad64f810` stands.
- D-012: presentation = P0 GL blit (Vulkan SDK absent); P1 Vulkan presentation owed, gated on SDK + DLSS milestone. D-013: CUDA bitmap HUD; ImGui rides P1.

## Current task — M2 `newton`

1. Contract first: `contracts/newton.contract.md` — spatial hash + leapfrog KDK + PM gravity (cuFFT Poisson) + P³M near-field; conservation gates (deterministic reductions via `core/lib/reduce.cuh`); scenario contracts `kepler` / `cloud-collapse` / `three-body`; oracle = nexus N2 + drift bounds.
2. The app grows its first *forces* — the galaxy stops being inert and starts to swirl for real.
3. First scenario goldens: (dials, seed, empty trace, steps) → declared state hash via liborrery envelope.

**M2 gate:** 1M gravitating particles at 60 Hz; energy/momentum/L drift within declared bounds over 10⁶ ticks; scenario goldens frozen and two-pass verified.

## Chores carried

- P1 Vulkan presentation (SDK install first) · ImGui at P1 · TAA · clang/g++ nexus parity · art-direction pass (dust lanes, interarm contrast, flare sprites).

## Standing context

Physics oracle: `tiny_nexus` v1.0.0 (golden `ad64f810`). Dials v0 frozen. Frame contract v1.0.0 frozen (god-hand, regime hex table). Repo docs are authoritative over agent memories where they disagree.
