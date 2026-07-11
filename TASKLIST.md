# TASKLIST.md — the build plan

Status legend: ☐ planned · ◐ in progress · ☑ done (golden frozen) · ✖ killed (DECISIONS entry).

## M0 · `nexus` — the composition proof *(current)*
- ☑ Contract drafted: `contracts/nexus.contract.md` + `nexus.schema.json` (v1.0.0-rc1)
- ☐ **Operator review of contract + dial defaults** ← gate, nothing proceeds past it
- ☐ Implement `nexus/tiny_nexus.cpp` (C++17, fp64, single file, MSVC/clang/g++ parity)
- ☐ Battery N1–N11 green; `--selftest`; golden frozen; 3× byte-identical
- **Gate:** all four regimes step in one scene without contradiction; dial defaults frozen; nexus becomes the standing oracle.

## M1 · `canvas` — first light
- ☐ frame.contract.md v1.0.0 (freeze the DRAFT)
- ☐ Vulkan swapchain + CUDA external-memory interop skeleton (simpleVulkan/Mímir reference)
- ☐ liborrery lift into `core/lib/` (verbatim, KAT selftest green, pinned ORRERY commit recorded)
- ☐ CINEMATIC stack: HDR accumulation → mip bloom → auto-exposure → AgX (+ ACES parity mode) → dither
- ☐ Two-stream frame loop (presentation never waits); ImGui HUD shell; damped camera
- **Gate:** 1M inert drifting particles at 60 Hz, 1080p, CINEMATIC §7 all boxes checked. The "prettiest N-body screensaver" — if this frame doesn't beat three.js, stop and fix before physics.

## M2 · `newton` — the classical tier
- ☐ Spatial hash; leapfrog KDK; PM gravity (cuFFT Poisson) + P³M near-field
- ☐ Conservation gates (energy/momentum/L drift bounds, deterministic reductions)
- ☐ Scenarios: `kepler` (two-body vs nexus N2), `cloud-collapse`, `three-body`
- **Gate:** 1M gravitating particles 60 Hz; drift within declared bounds over 10⁶ ticks; goldens frozen.

## M3 · `arrow` — thermodynamics + inscription
- ☐ T/P/S meters (dS/dτ display only); mixing scenarios show the arrow emerge
- ☐ Ratchet mechanic: record events into environment DOF, per-bubble redundancy R, collapse past `(1−p)ρ = p` (mapping contract; params anchored to ORRERY receipts)
- **Gate:** nexus N6 parity in-sim at small config; the double-slit detector visibly *writes records* before collapse.

## M4 · `einstein` — the relativistic layer
- ☐ Rapidity-form p = γmv integrator; per-particle τ driving all animation; photons at c
- ☐ 1PN precession from Φ; 2.5PN drag (binaries inspiral); Kepler-at-t_emit observation
- **Gate:** nexus N3/N7/N8/N9 parity at declared tolerance; a clock pair visibly desyncs across a potential well.

## M5 · `gargantua` — compact objects
- ☐ Collapse condition → BH entity; Paczyński–Wiita disk dynamics; Hawking evaporation (T ∝ ħ/M)
- ☐ GARGANTUA Kerr renderer lifted for near-BH view (dynamic M, a); OptiX-vs-compute geodesic spike (D-007, ADR-007 protocol — result reported back to ORRERY `lens`)
- **Gate:** cloud → star → remnant → BH ladder runs unscripted; a small BH evaporates on camera; lensing oracle-checked against analytic Kerr.

## M6 · `planck` — quantum bubbles
- ☐ Split-step cuFFT bubbles (64³ first), spawn-on-isolation, inscription collapse
- ☐ Scenarios: `doubleslit` (build-your-own-detector), `tunneling`, `sho-eigenstates`
- **Gate:** nexus N5 parity; double-slit fringes emerge from single collapses; 16 GB budget measured (Q-004 resolved by data).

## M7 · `cosmos` — the tiny planet
- ☐ 3-torus wrap (forces + light); light-history ring buffer ("see your own past")
- ☐ Stereographic little-planet projection; scale factor a(t) cosmology mode
- **Gate:** the signature screenshot — the universe as a globe, lensed, self-visible.

## Backlog (measured adoption only; each needs an honest baseline)
- ☐ DLSS 4.5 via Streamline · ☐ RTXDI/ReSTIR emissive-particle lighting · ☐ NRD · ☐ Slang port of shared math (D-008) · ☐ UE 5.8 TextureShare showcase shell · ☐ MCP surface for scenario driving · ☐ df64 zoom ladder (D-010 gate)
