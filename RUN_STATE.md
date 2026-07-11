# RUN_STATE.md

**As of:** 2026-07-11 · **Milestone:** M4 `einstein` · **State:** M3 CLOSED — the arrow of time is a frozen golden. **9/9 goldens GREEN in 55 s via `harness/verify.py`.**

## What M3 established

- **Entropy meters** (coarse 32³ x/v histograms, integer-atomic, declared machinery) with the ledger rule: gates on ΔS, never absolute level.
- **The arrow, three ways**: `merger` (S rises 0.90 nat through violent relaxation, golden `34a2db77`); `echo` (Loschmidt: momentum flip at tick 6000 retraces the rise EXACTLY — S 14.084043 → 14.570027 → 14.084043, golden `2f02d94f`); `ratchet` (in-sim inscription engine matches the Mirror/N6 closed form to 0.33% worst-class, golden `ccf4a3f8`).
- **The detector writes records**: transit → deterministic inscription → redundancy walk → RECORDED latch (regime 0x40, blue-tinted in render). Golden `a0c31f74`. This is M6's double-slit collapse mechanic, live and gated, waiting for ψ.
- D-015: merger mono gate corrected (real transient compressions, not jitter); SOLV_NONE; latch semantics.

## Current task — M4 `einstein` (the relativistic layer)

Contract first (`contracts/einstein.contract.md`), then:
1. Rapidity-form p = γmv integration (the galaxy already runs at 0.16c — 97.5% REL — Newtonian p=mv is now the *approximation being retired*); τ gains the potential term via Φ_PM: dτ = dt(1 − v²/2c² + Φ/c²).
2. Photons (massless particles at c); 1PN correction from the Φ grid (nexus N3 parity gate); Kepler-at-t_emit observation for the renderer.
3. 2.5PN drag scenario (binary inspiral) if budget allows — else M5 with the BHs.
4. Oracle gates: nexus N3/N7/N8/N9 parities at declared fp32 tolerances.

**M4 gate:** relativistic kepler scenario (compactness ≈ 5e-3) precesses at the nexus N3 rate within tolerance; a clock-pair scenario shows τ desync across a potential well; photons hold |v|=c through the PM field; goldens frozen, harness green.

## Chores carried

P1 Vulkan (SDK) · ImGui · TAA · clang/g++ nexus parity · art pass · cufft64 dll packaging · P³M/spatial hash (close encounters).

## Standing context

Oracle: tiny_nexus `ad64f810`. 9 goldens: nexus + kepler `448847ec` · threebody `f67abce4` · cloud `975389db` · galaxy `ff17431a` · merger `34a2db77` · echo `2f02d94f` · ratchet `ccf4a3f8` · detector `a0c31f74`. Contracts: nexus v1.0.0, frame v1.0.0, newton v1.0.0, arrow v1.0.0. Build: nvcc + envelope.cpp + cufft.lib (BUILD.md app line still needs updating — fold into M4 commit). Repo docs authoritative over agent memories.
