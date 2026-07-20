# MODULE — interop (R0: the Vulkan⇄CUDA presentation rung)

**Purpose.** Axis B unparked: a live Vulkan window presenting a CUDA-rendered image
**zero-copy** — shared device-local `VkBuffer` (exported `OPAQUE_WIN32`, imported via
`cudaImportExternalMemory`), presented by `vkCmdCopyBufferToImage`, ordered by two
exported **timeline semaphores** (cudaDone/vkDone; no `vkQueueWaitIdle` in the steady
loop). R0 proves the MECHANISM with a deterministic tick-parameterized test field; the
CINEMATIC stack (R1) and the substrate wiring (R2) ride this path, as will
DLSS/Streamline/RTXDI later.
**Contract:** `contracts/interop.contract.md` v1.0.0 (FROZEN; approved + implemented
2026-07-19, D-034 — the Vulkan SDK gate cleared the same session). **Oracle:** the
Khronos validation layer (zero errors AND zero warnings, gated) + byte-identity
round-trip (the CUDA view of the shared bytes memcmp'd against the Vulkan-copied
bytes). Device selection by **LUID match** against the CUDA device — the Intel iGPU
also enumerates on this machine (D-034 measured); index-0 selection is a named defect.
**Golden (1/1, GPU):** `interop_r0` `4ba7fbcb` (`--headless --frames 240 --golden`,
~1.2 s: G-ROUNDTRIP · G-SYNC monotonic embedded counters · G-DEVICE LUID ·
G-VALIDATION clean; declared path `buffer-copy`). `--selftest` = capability probe
(< 2 s). Windowed face (default) smoke-verified live (window up + presenting, FIFO
vsync; pacing/fps non-declared per the frame contract).

```
nvcc -O3 -arch=sm_89 -o build\interop.exe render\interop.cu -I"C:\VulkanSDK\1.4.350.0\Include" "C:\VulkanSDK\1.4.350.0\Lib\vulkan-1.lib" user32.lib
```

**Honest boundary:** synthetic pattern only (no substrate wiring — R2); fixed 1280×720,
no resize, no HDR; Windows + discrete-NVIDIA only; the windowed blit puts RGBA bytes in
a BGRA swapchain image (R/B swap visible in the test pattern — an accepted R0
cosmetic, exact swizzle lands with the R1 CINEMATIC port; the DECLARED headless path
is byte-exact and swapchain-free). CINEMATIC.md binds from R1, not R0.

# MODULE — cinematic (R1: the CINEMATIC post-chain, SUPERNOVA scene)

**Purpose.** The full CINEMATIC §1 pipeline on the R0 path: fp32 linear accumulation
(≥RGBA16F floor) → mip-chain threshold-free bloom (CoD 13-tap down / 3×3 tent up, 6
mips) → **energy-weighted log₂-luminance auto-exposure** (EV, tick-smoothed τ=90; the
astro meter — see the contract's freeze amendment for the two measured failures of the
percentile meter) → minimal-AgX tone map (ACES-Narkowicz parity on `T`) → cinematic
mode: astro-stretch W=40 + grain 0.3% → triangular dither → the ONLY sRGB encode →
BGRA pack (the R0 swizzle fixed). Scene: ~2500 counter-hashed blackbody stars
(Tanner–Helland, 2.2-linearized) + the SUPERNOVA (star 0: ×10⁵ hot-blue flash →
cooling decay, 1200-tick loop, pure function of the tick). Splats are analytic
gaussians σ=1.35 px (band-limited — the R1 AA resolution). NO float atomics in the
declared path (fixed-order per-pixel gather; histogram = integer atomics).
**Contract:** `contracts/cinematic.contract.md` v1.0.0 (FROZEN). **Oracle:** KAT
selftest (blackbody triples · AgX monotone/range · sRGB roundtrip 2e-6 · triangular
dither statistics at 4σ) + validation layers + the structural gates.
**Golden (1/1, GPU):** `cinematic_r1` `4962558c` (`--headless --frames 240 --golden`,
~2.4 s: G-RANGE 2^25.75 ≥ 2^10 · G-EXPOSURE EV moved 4.30 ≥ 1.5 across the flash ·
G-VALIDATION 0/0 · LUID). `--shot` frames: `runs/r1_supernova_flash.png` (blinding) +
`runs/r1_supernova_decay.png` (after-flash blindness re-adapting) — the judged half.
Windowed face: drag-orbit (critically-damped spring k=10, §5), `T` AgX/ACES/raw-clamp,
`C` cinematic/physical, `D` dither, HUD (minimal 5×7 font), title shows SN phase + EV.

```
nvcc -O3 -arch=sm_89 -o build\cinematic.exe render\cinematic.cu -I"C:\VulkanSDK\1.4.350.0\Include" "C:\VulkanSDK\1.4.350.0\Lib\vulkan-1.lib" user32.lib
```

**Honest boundary:** synthetic (physical-blackbody) scene — the substrate arrives at
R2 with pipeline + golden retained; fixed 1280×720; AA = band-limited splats (TAA/DLSS
at the DLSS milestone); LUT grading deferred; wheel-zoom and crossfade deferred (single
scene, fixed distance); manual-EV keys deferred to the R2 polish (auto + the C/T/D
toggles ship now). The judged §7 boxes await operator sign-off (recorded in TASKLIST).
