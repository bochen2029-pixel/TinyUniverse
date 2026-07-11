# TASKLIST.md — the build plan

Status legend: ☐ planned · ◐ in progress · ☑ done (golden frozen) · ✖ killed (DECISIONS entry).

## M0 · `nexus` — the composition proof ☑ **CLOSED 2026-07-11**
- ☑ Contract drafted + operator-approved (v1.0.0; errata D-011 applied at implementation)
- ☑ `nexus/tiny_nexus.cpp` implemented (C++17, fp64, single file, MSVC-built)
- ☑ Battery N1–N11 green in 24.1 s; golden frozen `ad64f810`; N11 + 3× out-of-process byte-identical
- ☑ Dial defaults v0 frozen; nexus is the standing oracle
- ☐ *(owed, moved to M1 chores)* clang++/g++ parity build — no non-MSVC toolchain on this machine
- **Gate: MET** — all four regimes stepped in one scene, zero contradictions (N10 masks 0x12/0x24/0x1C stable, drifts ≤ 2.4e-13).

## M1 · `canvas` — first light ☑ **CLOSED 2026-07-11** (P1 presentation owed)
- ☑ frame.contract.md v1.0.0 frozen (Q-001 → god-hand; regime hex table; D-012 P0/P1 clause)
- ☑ liborrery lifted verbatim (`core/lib/`, ORRERY commit d56c4c7, per-file sha256 in LIFT.md); BLAKE2b byte-compat verified (`harness/hash_compat.cpp` green) — **nexus golden stands**
- ☑ CINEMATIC stack in CUDA: float4 HDR → 13-tap mip bloom → log-lum auto-exposure → AgX + ACES parity → astro stretch (cinematic-only, HUD-declared) → triangular dither
- ☑ Two-stream frame loop (sim low-prio ping-pong publish, present high-prio); damped spring camera; CUDA bitmap HUD (D-013; ImGui deferred to P1)
- ☑ Presentation P0: thin GL blit (D-012 — Vulkan SDK absent, measured)
- ☐ *(owed → backlog)* P1 Vulkan swapchain + external-memory presentation (gate: SDK install + DLSS milestone)
- ☐ *(carried)* clang++/g++ parity build of nexus
- **Gate: MET** — 1M particles 1080p: **499 fps avg / 226 min** (SSAA 2×: 178/152); CINEMATIC §7 **10/10** (evidence: `app/MODULE.md`, `runs/firstlight_*.png`).

## M2 · `newton` — the classical tier ☑ **CLOSED 2026-07-11**
- ☑ PM gravity: 128³ CIC (fixed-point uint64 deposit — Invariant 4) → cuFFT Poisson → FD force grid → gather; periodic box = the torus, natively
- ☑ KDK leapfrog + Kahan drift; tiny solver (N≤32) fp64-internal (D-014); direct solver ≤4096
- ☑ Conservation meters via `fixed_atomic_add` (device) / fp64 (tiny); HUD drift readouts
- ☑ Scenario contracts + envelope face (liborrery `full_envelope`/`golden_check`): kepler · threebody · cloud · galaxy — **all gates green, goldens frozen + GOLDEN OK on re-run** (kepler de<1e-6 @10⁶ ticks; figure-8 bounded 29.5 periods; cloud 0.46%; galaxy 1.9%, p_drift 2e-6)
- ☐ *(deferred → M3)* P³M short-range correction (PM cell softening is the declared v1 resolution); spatial hash (arrives with short-range forces); harness/verify.py
- **Gate: MET** — 1M gravitating @1080p: **347 fps avg / 159 min**; 10⁶-tick drift evidence in goldens (kepler/threebody run 10⁶ ticks as their golden params).

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
