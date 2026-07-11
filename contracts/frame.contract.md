# CONTRACT — frame · v1.0.0 · status: FROZEN 2026-07-11 (Q-001 resolved: god-hand v1)

**The sacred seam.** Everything the renderer may read, everything the input system may write, and the sync semantics between them. The renderer is a *client*: it never mutates physics state; physics never depends on frame rate. Breaking changes here are MAJOR bumps and ripple into every golden.

## Presentation layer (D-012)

All **rendering** is CUDA kernels (accumulation, bloom, exposure, tone map, HUD) ending in an 8-bit RGBA image. Only *presentation* of that finished image touches a graphics API:

- **P0 (current):** thin OpenGL blit — `cudaGraphicsGLRegisterImage` surface/array copy + textured quad (GARGANTUA-proven, zero SDK dependencies). ~80 lines; nothing above the final image touches GL.
- **P1 (committed target):** Vulkan swapchain with external-memory images + timeline semaphores (`cudaImportExternalMemory`/`cudaImportExternalSemaphore`), required for DLSS/Streamline. Gate to swap: Vulkan SDK installed + the DLSS milestone (backlog). Q-005's Vulkan-vs-DX12 door remains closed on Vulkan.

## Ownership & streams

Two CUDA streams: `simStream` (low priority, fixed-dt ticks) and `presentStream` (high priority) — **presentation never waits on accumulation** beyond the published-tick event. Particle state is double-buffered ping-pong: a tick reads buffer A, writes B, then publishes B (event + index flip). The renderer reads the last *published* buffer only. No device-wide syncs in the frame loop.

## Buffers (SoA; capacity N_max = 4,194,304 default; live count published per tick)

| name | type | writer | semantics |
|---|---|---|---|
| pos | float4[N] ×2 (ping-pong) | sim | xyz = position (su); **w = emitT** (blackbody emission temperature; 0 = non-emissive) |
| mom | float4[N] | sim | xyz = relativistic momentum p = γmv (sm·su/s); **w = mass** (sm; 0 = photon) |
| tau | float[N] | sim | proper time (s) — drives ALL visible animation rates |
| regime | u32[N] | sim | derived-only bitmask (canonical table below) |
| bubble | u32[N] | sim | ψ-handle, 0 = none (allocated from M6; bubbleVis K = 4096/bubble) |
| phi | 3D tex r32f | sim | PM potential (from M2) — renderer uses for lensing/redshift; sim for τ/1PN |
| events | append u64[] | sim | declared log: inscriptions, collapses, ignitions, captures, mergers (from M3) |
| inputTrace | append records | input | THE determinism boundary — sim consumes only this, never raw input |

## Canonical regime bitmask (frozen; values proven in nexus N10)

| bit | name | derivation |
|---|---|---|
| 0x01 | QUANTUM | has ψ (bubble ≠ 0) |
| 0x02 | CLASSICAL | default for massive, v ≤ 0.05c |
| 0x04 | REL | v > 0.05c, or massless |
| 0x08 | COMPACT | within 10 r_s of a compact object |
| 0x10 | BOUND | pair/system energy < 0 |
| 0x20 | MASSLESS | mass = 0 (photons) |
| 0x40 | RECORDED | reserved (M3: redundancy R past threshold) |
| 0x80 | INSCRIBING | reserved (M3: record events in flight this tick) |

Derived-only (liborrery `regime.h` pattern): computed from state every tick, never stored authoritative.

## Input trace v1 (god-hand; Q-001 resolved to default)

Record = `{ tick u64, tool u16, flags u16, pos float3, dir float3, mag float }` (40 bytes, packed).
Tool enum v1: `0 NONE · 1 SPAWN · 2 FLING · 3 COMPRESS · 4 HEAT · 5 COOL · 6 DETECTOR (M3+) · 7 CAMERA_MARK (non-declared bookmark)`.
Replay of (dials, seed, trace, steps) is bit-exact (Invariant 2). Camera motion is NOT in the trace (visuals are non-declared); CAMERA_MARK exists only for replay authoring convenience.

## The clock rule

`dt` is a dial (1/240 s). Sim advances in whole ticks; the render frame may interpolate *visually only* (never feeds back). Pause/slow-mo scale tick *dispatch rate*, never dt. Tick budget per frame is bounded (≤ 8) — a slow render slows sim time rather than exploding the tick queue (non-declared behavior; replay is per-tick, unaffected).

## Declared vs non-declared

Declared (golden-hashable): pos/mom/tau/regime at tick T, conserved-quantity meters, event log, live count. Non-declared: every pixel, bloom/exposure state, HUD, frame timing, bubbleVis sampling order, camera.

## Changelog

- **v1.0.0 (2026-07-11)** — frozen at M1: Q-001 → god-hand v1; regime canonical hex table (nexus N10 values); emitT packed into pos.w, mass into mom.w; presentation-layer clause added (D-012: P0 GL blit interim, P1 Vulkan committed); tick budget rule.
- **v0.1.0-DRAFT (2026-07-11)** — shape review draft.
