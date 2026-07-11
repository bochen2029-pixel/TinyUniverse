# BUILD.md — pinned toolchain

**Pins (2026-07-11):** CUDA 13.1 (V13.1.80) · MSVC 2022 Community (vcvars64) · GPU RTX 4070 Ti SUPER, 16 GB, `sm_89` · Windows 11.

## Golden path — single-file module (nexus, early CUDA modules)

CPU (nexus — any C++17, MSVC shown; clang++/g++ parity required like ASTRA-7):

```
cmd /c '"C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat" >nul 2>&1 && cl /std:c++17 /EHsc /O2 /W4 nexus\tiny_nexus.cpp /Fe:build\tiny_nexus.exe'
```

CUDA — the app (from M2 it links liborrery's envelope + cuFFT; cufft64 runtime dll comes from the toolkit redist, D-014):

```
cmd /c '"C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat" >nul 2>&1 && nvcc -O3 -arch=sm_89 -Xcompiler "/O2" -o build\tinyuniverse.exe app\tinyuniverse.cu core\lib\envelope.cpp user32.lib gdi32.lib opengl32.lib cufft.lib'
```

## CMake preset (multi-file; app from M1)

Mirrors ORRERY D-021: **fat binary** — sm_89 + sm_90 SASS, compute_120 PTX (forward JIT) · static cudart + static MSVC runtime (self-contained dist, Buddhabrot pattern) · `--use_fast_math` **banned** in declared paths · goldens stay hardware-pinned to sm_89 (re-baseline = operator-signed note in goldens/).

App dependencies via FetchContent (Buddhabrot-proven set): GLFW 3.4, GLAD v2, Dear ImGui 1.91, stb. Vulkan SDK: required from M1 — pin exact version in this file at M1 (D-002/Q-005).

## OptiX (M5)

Function-table dispatch via `optix_stubs.h` (`optixInit()` pulls `nvoptix.dll` from the driver — **no .lib link**). CMake auto-detects the SDK at `C:/ProgramData/NVIDIA Corporation/OptiX SDK 9.1.0/include`; if absent, the OptiX path is compiled out and the renderer falls back to the CUDA-compute geodesic path. The app must never *require* the SDK. (Buddhabrot v4 recipe, verbatim.)

## Rules

- Every module's MODULE.md records its exact build command; a fresh agent must compile any module from this file + MODULE.md alone, one pass.
- Selftests < 30 s per module; full golden suite < 5 min.
- Host≠device libm: pin transcendental KAT values separately per side (liborrery lesson).
