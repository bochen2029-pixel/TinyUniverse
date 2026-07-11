# MODULE — app (tinyuniverse v0.2.0 "newton")

## M2 addendum (2026-07-11)

v0.2.0 adds self-gravity per `contracts/newton.contract.md`: PM 128³ (fixed-point CIC → cuFFT → FD force grid), KDK + Kahan drift, fp64 tiny solver (D-014), conservation meters, four scenarios with frozen goldens (kepler `448847ec` · threebody `f67abce4` · cloud `975389db` · galaxy `ff17431a`, all `GOLDEN OK` on re-run). Headless face is envelope-conformant: `--scenario X [--steps S] [--seed N] --json | --golden`; `--warm T` pre-evolves before `--shot`. Gate: **347 fps avg / 159 min** @1M/1080p with live gravity. Build now needs `core\lib\envelope.cpp` + `cufft.lib`:

```
nvcc -O3 -arch=sm_89 -Xcompiler "/O2" -o build\tinyuniverse.exe app\tinyuniverse.cu core\lib\envelope.cpp user32.lib gdi32.lib opengl32.lib cufft.lib
```

(Runtime dependency: cufft64 dll from the toolkit redist — D-014.) Everything below documents v0.1.0 first light and remains true.

---

# MODULE — app (tinyuniverse v0.1.0 "first light")

**Purpose:** M1 canvas — the windowed client over the frame contract, and the CINEMATIC stack's reference implementation. 1M inert drifting particles (galaxy scenario), rendered entirely in CUDA; presentation is a thin GL blit (D-012 P0).

**Contracts:** `contracts/frame.contract.md` v1.0.0 · gated by `CINEMATIC.md` §7.
**liborrery consumer #1:** `core/lib/rng.cuh` — all init randomness is counter-keyed `f(seed, particle, dim)`.
**Status:** M1 GATE MET 2026-07-11. Visuals are non-declared output — no golden; the gate evidence is below.

## Usage

```
build\tinyuniverse.exe [--w 1920 --h 1080] [--n 1000000] [--seed 20260711]
                       [--ssaa] [--phys] [--aces]
                       [--bench FRAMES]                 # windowed, vsync off, fps stats
                       [--shot PATH.bmp [--frames N]]   # headless render-to-BMP
```
Keys: drag orbit · wheel zoom · 1-3 presets · O auto-orbit · P pause · C cinematic/physical · T AgX/ACES · B bloom · A auto-exposure · +/− EV · F1 HUD · V vsync · Esc.

Build (BUILD.md; **no fast-math** — Invariant 6):
```
nvcc -O3 -arch=sm_89 -Xcompiler "/O2" -o build\tinyuniverse.exe app\tinyuniverse.cu user32.lib gdi32.lib opengl32.lib
```

## M1 gate evidence (2026-07-11, RTX 4070 Ti SUPER)

| gate | result |
|---|---|
| 1M particles, 1080p, 60 Hz | **499.4 fps avg (2.00 ms), min 225.9** — 8× over gate |
| quality ladder (SSAA 2×, 4K internal) | **177.8 fps avg, min 151.8** |
| beauty | `runs/firstlight_cinematic.png` (AgX + astro stretch) · `runs/firstlight_physical.png` (honest T⁴) · `runs/firstlight_ssaa.png` |

CINEMATIC §7 checklist: **10/10.** HDR float4 linear end-to-end, single sRGB encode ✓ · AgX default + ACES-Narkowicz parity, after bloom ✓ · mip-chain threshold-free bloom (13-tap down, 3×3 tent up, 6 mips) ✓ · EV exposure, manual + auto (log₂-lum histogram, 2–98% clip, ~0.8 s adaptation) ✓ · blackbody-only color (Tanner-Helland; g-factor arrives with M4 physics) ✓ · SSAA 2× dev quality (TAA/DLSS backlog) ✓ · triangular dither pre-quantization ✓ · critically-damped spring camera (k=90, c=19), damped preset transitions, no snapping ✓ · cinematic/physical toggle with HUD-enumerated deviations ✓ · typed HUD meters, physics cyan / render amber (CUDA 5×7 bitmap font; ImGui deferred, D-013) ✓.

## Frozen art constants (cinematic mode; changing these is an art decision, log it)

- Auto-exposure key **0.06** · bloom mix **0.035** · saturation **1.20**
- Astro stretch: piecewise luminance compression, identity for L ≤ 1, extended Reinhard on (L−1) with W=40 above — the CINEMATIC §4 extreme-range clause; enumerated on HUD as `STRETCH`
- Luminosity law: cinematic L ∝ (T/5800)^2.2 · physical L ∝ (T/5800)^4, no stretch, sat 1.0 (`NO CLIP` on HUD)
- Scenario: 68% spiral disk (log-spiral arms, pitch 0.24, inner cutoff 14 su, 7% blue giants 9500–15000 K), 12% bulge (σ=10), 20% halo shell (40–260 su); camera presets wide/close/top

## Internal design notes

- **Two streams** (frame contract): simStream low-priority runs fixed-dt ticks (dt = 1/240, ≤ 8/frame) with ping-pong pos buffers + published-index flip; presentStream high-priority waits only on the tick-done event. No device-wide syncs in the loop (the one `cudaStreamSynchronize(prsS)` is the exposure-histogram readback, on the present stream itself).
- Splat: perspective project + 2×2 bilinear additive `atomicAdd` (float atomics are fine here — rendering is non-declared; Invariant 4 governs declared reductions only).
- Presentation (D-012 P0): `cudaGraphicsGLRegisterImage` → `cudaMemcpy2DToArray` → textured quad. Everything above the finished 8-bit image is CUDA and survives the P1 Vulkan swap untouched.
- Headless `--shot` never creates a window or GL context — the same binary is the instrument face (D-003).

## Owed / backlog (tracked in TASKLIST)

P1 Vulkan presentation (gated on SDK install + DLSS milestone) · ImGui HUD (D-013) · TAA · art-direction pass (dust lanes as extinction, deeper interarm contrast, star flare sprites) · torus wrap + little-planet projection (M7).
