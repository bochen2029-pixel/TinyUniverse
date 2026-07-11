# RUN_STATE.md

**As of:** 2026-07-11 · **Milestone:** M5 `gargantua` · **State:** M4 CLOSED — relativity is always-on. **12/12 goldens GREEN in 66 s.**

## What M4 established

- Every particle moves relativistically: v = p/√(m² + |p|²/c²) (v < c structural), photons at c with factor-2 bending, proper time dτ = dt√(1 + (2Φ − v²)/c²) — **photons don't age** (the clamp). Tiny solver fp64 u-form.
- Oracles: keprel vs exact Sommerfeld 0.50% (`f985e473`) · clocks 6.3e-4 (`330c86a7`) · photons 0.83% (`c4c565de`).
- **D-016**: the 1PN field term was withdrawn after measurement contradicted the 7π superposition claim (6.41π measured; Q-006 holds the derivation). Strong-field precession arrives with M5's oracles. All 8 prior goldens superseded + re-frozen; the echo still reverses exactly under relativistic KDK.

## Current task — M5 `gargantua` (compact objects)

Contract first (`contracts/gargantua.contract.md`), then:
1. Collapse condition (cluster inside its own r_s) → BH entity (absorbing, mass-accumulating); Paczyński–Wiita Φ = −GM/(r − r_s) for near-BH dynamics — ISCO at 6GM/c² is the oracle (nexus N4).
2. Hawking evaporation: T ∝ ħ/M — small BHs visibly evaporate (emit particles/photons; energy bookkeeping declared).
3. GARGANTUA lift: the validated Kerr geodesic renderer (C:\blackhole\files\blackhole.cu) as the near-BH view; OptiX-vs-compute geodesic spike (D-007, ADR-007 protocol — result reported back to ORRERY `lens`).
4. Scenarios: `collapse` (cloud → BH), `isco` (test orbits at/inside 6GM/c², P–W oracle), `hawking` (small BH evaporates on camera); 2.5PN inspiral if in budget.

**M5 gate:** cloud → star → remnant → BH ladder runs unscripted; ISCO golden vs P–W analytic; a small BH evaporates with declared energy bookkeeping; 15/15-ish goldens green.

## Chores carried

P1 Vulkan (SDK) · ImGui · TAA · clang/g++ nexus parity · art pass · cufft64 dll packaging · P³M/spatial hash · BUILD.md app-line update (nvcc + envelope.cpp + cufft.lib) · Q-006 derivation.

## Standing context

Oracle: tiny_nexus `ad64f810`. 12 goldens (hashes in D-016). Contracts: nexus/frame/newton/arrow v1.0.0, einstein v1.0.1. App v0.2.x: 11 scenarios, envelope face, harness green. Repo docs authoritative over agent memories.
