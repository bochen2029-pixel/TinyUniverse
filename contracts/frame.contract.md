# CONTRACT — frame · v0.1.0-DRAFT · status: SHAPE REVIEW (freezes as v1.0.0 at M1)

**The sacred seam.** Everything the renderer may read, everything the input system may write, and the sync semantics between them. The renderer is a *client*: it never mutates physics state; physics never depends on frame rate.

## Ownership & interop

All buffers Vulkan-allocated with `VK_EXTERNAL_MEMORY_HANDLE_TYPE_OPAQUE_WIN32_BIT`, imported into CUDA via `cudaImportExternalMemory` (zero copy). Two timeline semaphores (`simDone`, `frameDone`) imported both ways. Two CUDA streams: `simStream` (low priority, fixed-dt ticks) and `presentStream` (high priority) — **presentation never waits on accumulation**; the renderer reads the last *published* tick's buffers (double-buffered publish index, atomically flipped at tick end). No device-wide syncs in the frame loop.

## Buffers (SoA; count = N_max, live count published per tick)

| name | type | writer | semantics |
|---|---|---|---|
| pos | float3[N] | sim | position, su, box-relative |
| mom | float3[N] | sim | relativistic momentum p = γmv, sm·su/s |
| mass | float[N] | sim | sm; 0 = photon |
| tau | float[N] | sim | proper time, s — drives ALL visible animation rates |
| emitT | float[N] | sim | emission temperature (blackbody input), K-analog; 0 = non-emissive |
| regime | u32[N] | sim | derived-only bitmask (QUANTUM·CLASSICAL·REL·COMPACT·BOUND·RECORDED…) — canonical hex table fixed at v1.0.0 |
| bubble | u32[N] | sim | ψ-handle, 0 = none; bubble grids are sim-internal (renderer gets density splats via `bubbleVis`) |
| bubbleVis | float4[B·K] | sim | per-bubble visualization samples (|ψ|² point cloud), renderer-consumable |
| phi | 3D tex (r32f) | sim | PM potential — renderer uses for lensing/redshift; sim uses for τ/1PN |
| events | append[u64] | sim | declared event log: inscriptions, collapses, ignitions, captures, mergers |
| inputTrace | append | input system | tick-stamped interventions (tool, position, magnitude) — THE determinism boundary: sim consumes only this, never raw input |

## The clock rule

`dt` is a dial (1/240 s). Sim advances in whole ticks; the render frame interpolates *visually only* (never feeds back). Pause/slow-mo scale the tick *dispatch rate*, never dt itself.

## Declared vs non-declared

Declared (golden-hashable): pos/mom/mass/tau/regime/events + conserved-quantity meters at tick T. Non-declared: everything the renderer derives (colors, bloom, exposure), bubbleVis sampling order, frame timing.

## Open until v1.0.0 (M1)

Input-trace schema detail (blocked on Q-001) · regime bitmask canonical hex values · N_max and bubbleVis K defaults · df64 ladder reservation (D-010) · exact semaphore protocol writeup.
