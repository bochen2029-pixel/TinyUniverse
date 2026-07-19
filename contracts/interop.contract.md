# interop.contract.md — `interop` (R0: the Vulkan⇄CUDA presentation rung)

**Version 1.0.0 — FROZEN** (approved by operator directive 2026-07-19: "keep going and
proceed to next" on the presented draft — the same green-light convention as the
choptuik contract; open questions resolved as recommended, below; frozen with golden
`interop_r0` `4ba7fbcb…`, two-passed same day; windowed face smoke-verified live). Prerequisite CLEARED: Vulkan SDK 1.4.350.0 pinned, device enumeration verified
(D-034). Basis: ROADMAP §4 R0 · D-002 (Path A: CUDA renders, Vulkan presents, renderer
swappable at the frame contract) · D-012 (P1 presentation: swapchain + external memory
+ timeline semaphores, replacing the P0 GL blit).

## Purpose

The smallest honest rung of Axis B: **a live Vulkan window presenting a CUDA-rendered
image zero-copy** (external memory), with the synchronization done right (timeline
semaphores) — "a window showing the sim breathe." R0 proves the MECHANISM; it wires no
substrate. Everything the north star needs later (CINEMATIC port at R1, DLSS/Streamline,
RTXDI) rides this presentation path.

## Module

`render/interop.cu` — single file, CUDA + Vulkan (via `%VULKAN_SDK%`), win32 surface.
No GLFW/ImGui at R0 (raw Win32 window, ~the P0 pattern); FetchContent stack arrives
with R1. Fixed 1280×720 at v0.1 (resize is declared OUT of R0 scope).

**The mechanism (the thing under contract):**
1. CUDA allocates the render target via external-memory-exportable allocation
   (`cudaExternalMemory` ⇄ `VK_KHR_external_memory_win32`); Vulkan imports it as a
   `VkImage` (or buffer→blit fallback if tiling import is unsupported — declared which).
2. Per frame: a CUDA kernel writes the image (R0 test field: a deterministic animated
   pattern — tick-parameterized wave/gradient, NO wall-clock) → timeline semaphore
   signals → Vulkan acquires, blits/copies to the swapchain, presents → semaphore
   returns the slot to CUDA. No CPU round-trip of pixels, no `vkQueueWaitIdle` in the
   steady loop.
3. **Device selection: explicitly the discrete NVIDIA device matching the CUDA device
   LUID** (`VK_KHR_device_group_creation`/`vkGetPhysicalDeviceProperties2` LUID ==
   `cudaDeviceProp::luid`) — the Intel iGPU also enumerates on this machine (D-034
   measured); index-0 selection is a named defect.

## Faces

- **windowed (default)** — the live window; ESC quits; title bar shows fps + frame
  counter (wall-clock display only — non-declared).
- **`--headless --frames N`** — no swapchain; the SAME external-memory path renders N
  frames offscreen and runs the gates below; prints the declared JSON + hash;
  exit 0/1/2 per house convention.
- **`--selftest`** — enumeration + capability probe (< 5 s): discrete device found by
  LUID match; `VK_KHR_external_memory_win32` + timeline semaphores supported; prints
  the chosen device name/driver.

## Gates (headless face; all must PASS for exit 0)

- **G-ROUNDTRIP** — the CUDA-written pattern, imported to Vulkan, copied back through
  the Vulkan side, is **byte-identical** to the CUDA source buffer (frame 0 and frame
  N−1). Zero-copy means zero corruption.
- **G-SYNC** — with a per-frame counter embedded in the image by the CUDA kernel, the
  Vulkan side observes strictly monotonic counters across N=240 frames (no torn/stale
  frame under semaphore ordering). Run with **validation layers ON
  (`VK_LAYER_KHRONOS_validation`): zero errors, zero warnings gated**.
- **G-DEVICE** — the selected Vulkan physical device LUID equals the CUDA device LUID
  (the discrete RTX card), asserted in the declared output.
- **G-DET** — the declared state (pattern-buffer hash at frames {0, N/2, N−1} +
  config) is byte-identical across two cold headless runs (two-pass at freeze).
  Golden: `goldens/interop/golden.hash`.
- **Recorded, NOT gated:** windowed fps (expect swapchain-limited/vsync; perf gates on
  presentation are flaky by design — the frame contract declares pacing non-declared).

## Oracle

The Vulkan validation layers (zero-diagnostics gate) + byte-identity round-trip. No
physics; no analytic oracle applies at R0.

## Honest boundary (declared limits)

- R0 presents a synthetic CUDA pattern — **no substrate wiring** (R1 `cinematic` ports
  the HDR→bloom→AgX stack; R2 wires the sim). No resize, no HDR swapchain, no DLSS.
- Windows + NVIDIA discrete only (`VK_KHR_external_memory_win32`); the iGPU is
  explicitly rejected, not supported.
- The CINEMATIC checklist does NOT apply to the R0 test pattern (it is a mechanism
  proof, not a shipped visual); it binds from R1.

## Determinism

Declared path = the headless face only: tick-parameterized pattern, fixed frame count,
no wall-clock, single queue, hash via the house blake2b idiom. The windowed face is
non-declared (presentation timing, vsync, monitor).

## Build (BUILD.md addendum on approval)

```
cl /std:c++17 /EHsc /O2 /W4 /I"%VULKAN_SDK%\Include" render\interop_host.cpp ... "%VULKAN_SDK%\Lib\vulkan-1.lib"
nvcc -O3 -arch=sm_89 ... render\interop.cu
```
(exact single-file vs two-TU split decided at implementation; recorded in MODULE.md.)

## Open questions — RESOLVED at approval (as recommended)

- Q-R0-1 → `B8G8R8A8_UNORM` swapchain at v0; HDR10 deferred to R1.
- Q-R0-2 → the shared allocation is a `VkBuffer` (exported `OPAQUE_WIN32`), presented
  by `vkCmdCopyBufferToImage` — the canonical CUDA-sample interop shape; the path is
  DECLARED in the golden string (`"path":"buffer-copy"`). Optimal-tiling image import
  is an R1 refinement, not an R0 requirement.
- Q-R0-3 → own binary `build\interop.exe` (smallest honest rung; zero blast radius on
  `tinyuniverse.exe`).
