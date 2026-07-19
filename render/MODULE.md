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
