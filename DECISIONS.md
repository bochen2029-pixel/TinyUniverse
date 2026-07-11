# DECISIONS.md — architecture decision records

Format: D-NNN · date · status · decision / why / consequences. Reversals get new entries, never edits.

---

**D-001 · 2026-07-11 · ACCEPTED — Sibling of ORRERY, not part of it.**
ORRERY Invariant 3 mandates headless tools; TU opens a window and runs on an interactive clock. TU adopts the doctrine (contract/golden/oracle/two-pass) as a sibling repo. ORRERY is referenced, never edited, from here. *Consequence:* TU carries its own harness, goldens, and contracts.

**D-002 · 2026-07-11 · ACCEPTED — Path A rendering stack: CUDA + Vulkan external-memory interop + hand-rolled CINEMATIC stack; renderer swappable at the frame contract.**
Custom beauty is a solved problem on this machine (GARGANTUA). UE 5.8 = TextureShare showcase shell only, someday (display-only; Beta since 4.27; UE6 migration clock). Falcor = reference quarry, not foundation (research-framework churn vs single-file doctrine). *Consequence:* CINEMATIC.md is law; frame contract is the swap seam.

**D-003 · 2026-07-11 · ACCEPTED — One binary, two faces (Buddhabrot v4 pattern).**
`tinyuniverse.exe` windowed by default; `--scenario/--steps/--seed/--json` (+`--selftest/--golden`) = headless instrument face with the universal envelope. *Consequence:* no separate tool exe; goldens run the same binary CI runs.

**D-004 · 2026-07-11 · ACCEPTED — Decoherence = Ratchet inscription, not generic CSL.**
Collapse when a bubble's state is redundantly recorded into environment DOF; unwrite bound `min(1, p/((1−p)ρ))^R`, crossover `(1−p)ρ = p` — the closed form ORRERY's `ratchet` tool froze (golden 91fce3c4; MC↔analytic 0.06%). Physically principled (quantum-Darwinism-shaped), measurable in-sim (R is visible!), parameter-anchored. *Consequence:* nexus N6 verifies the math; the sim mapping (what counts as a record) is M3's contract. No metaphysical commitment implied (Q-003).

**D-005 · 2026-07-11 · ACCEPTED — liborrery lifted verbatim.**
`rng.cuh` (counter-based stateless RNG), `reduce.cuh` (deterministic reductions), envelope (blake2b goldens), `ckpt.h` (save format), `regime.h` (derived-only bitmask) copied byte-identical from `C:\ORRERY\lib` at a pinned commit, KAT selftest included. *Consequence:* any divergence is a named fork here; TU inherits the host≠device libm pins.

**D-006 · 2026-07-11 · ACCEPTED — Determinism via fixed timestep + recorded input traces.**
dt = 1/240 s is a dial; sim ticks decouple from render frames; user interventions are recorded per-tick and replayable bit-exact. Golden = (dials, seed, trace, steps) → declared-state hash. Visuals/timing non-declared. *Consequence:* replays, goldens, and netplay-grade reproducibility come from one mechanism.

**D-007 · 2026-07-11 · ACCEPTED — OptiX via function-table dispatch, compile-out fallback (Buddhabrot recipe).**
No .lib link; `nvoptix.dll` from driver; SDK auto-detected, never required. *Consequence:* the M5 RT-vs-compute geodesic measurement (ORRERY's parked `lens` spike, ADR-007 protocol) can ship either verdict without a hard dependency.

**D-008 · 2026-07-11 · ACCEPTED (adoption deferred to M4+) — Slang for shared math.**
Metric/potential/blackbody/redshift math written once, emitted as CUDA (kernels) and SPIR-V (shaders). Khronos-hosted; Falcor-proven. *Consequence:* until adoption, shared math lives in a header included by both sides; the header is written Slang-portable (no CUDA-isms in shared functions).

**D-009 · 2026-07-11 · ACCEPTED — CPU nexus before any CUDA (ASTRA-7 pattern).**
`tiny_nexus.cpp`: fp64, single file, ordered PASS/FAIL battery proving the four regimes compose and freezing dial defaults. It then serves as the standing oracle for every CUDA module. *Consequence:* M0 blocks everything; that is the point.

**D-010 · 2026-07-11 · ACCEPTED — fp32 core, df64 ladder only on measured need.**
Box is 512 su; fp32 holds ~7 digits — fine at world scale. If zoom-as-mechanic exceeds it, adopt Buddhabrot's df64 ladder (and only then consider AstraCoord sectoring). *Consequence:* no premature 128-bit coordinates; the gate is a measured artifact, not a hunch.

**D-012 · 2026-07-11 · ACCEPTED — presentation layer: thin GL blit now (P0), Vulkan external-memory presentation when the SDK lands (P1).**
Measured fact: no Vulkan SDK on this machine (`C:\VulkanSDK` absent, `VULKAN_SDK` unset). All *rendering* is CUDA kernels ending in a finished 8-bit image either way; only presentation touches a graphics API. P0 = GARGANTUA's proven `cudaGraphicsGLRegisterImage` blit (~80 lines, zero new dependencies); P1 = Vulkan swapchain + external memory + timeline semaphores, required for DLSS/Streamline, adopted when the SDK is installed (backlog, with the DLSS milestone). This is a presentation deferral, not a renderer decision reversal — D-002 stands; Q-005's Vulkan-vs-DX12 door stays closed on Vulkan. Also: GARGANTUA's build used `-use_fast_math`; TU builds without it (Invariant 6), accepted perf cost ~nil (workload is memory/atomic-bound). *Consequence:* frame.contract.md v1.0.0 carries the P0/P1 clause; M1's gate is judged on the CUDA pipeline, which survives the P1 swap untouched.

**D-013 · 2026-07-11 · ACCEPTED — HUD v1 is a CUDA bitmap-font kernel, not ImGui.**
The typed-meters CINEMATIC box needs on-screen text; ImGui wants GLFW + a real GL/Vulkan pipeline + CMake — none of which the single-file P0 app has or needs yet. A 5×7 bitmap font drawn by a CUDA kernel (44 glyphs, 2× scale, drop shadow, physics-cyan/render-amber typing) satisfies the box in ~120 lines with zero dependencies. ImGui arrives with the P1 Vulkan presentation swap (it rides the same pipeline work). *Consequence:* HUD strings are uppercase-limited (glyph set); fine for meters.

**D-011 · 2026-07-11 · ACCEPTED — nexus contract errata: two tolerances corrected by fp64 analysis at implementation.**
(1) N2 |ΔL/L| 1e-12 → 1e-11: leapfrog conserves L exactly in exact arithmetic, but 6.2×10⁷ steps of fp64 roundoff random-walk to ~10⁻¹² (measured 5.8e-13) — the original gate tested the FPU, not the integrator. (2) N7: the 1e-12 cosh↔naive agreement claim is unattainable at γ=10³ (the naive path's own cancellation is ~10⁻¹⁰ there) — domain corrected to γ ≤ 10²; and the "> 1e-3 at γ=10⁷" cancellation demonstration was rounding-luck-dependent (first run measured 4.0e-4 — lucky libm) — replaced by a deterministic form: ladder extended to γ=10⁸ where β rounds to exactly 1.0 and the naive path diverges, plus a 10⁶× error-ratio gate. *Why logged:* both errata are exactly the class of confident-wrong the doctrine exists to catch — the spec claimed precision fp64 cannot deliver; the fix is in the contract changelog, pre-golden, MINOR-scope. *Consequence:* contract v1.0.0 frozen with them applied; golden `ad64f810`.
