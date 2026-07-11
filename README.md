# TINY UNIVERSE

**A scale-compressed universe: one particle substrate where quantum, Newtonian, relativistic, and black-hole physics emerge a screen-zoom apart — native Windows, C/C++/CUDA, RTX-leveraged, no browser.**

> The Insta360 "tiny planet" effect bends a flat neighborhood into a whole globe. Tiny Universe does that to physics: retune ħ (up), c (down), G (up) so the crossovers that nature spreads across 60 orders of magnitude all happen between 1 and ~10⁶ particles — a range one GPU actually simulates. One lone particle is visibly a wave. A hundred decohere into billiards. Ten thousand ignite into a star. A hundred thousand collapse into a black hole that lenses the sky and evaporates in Hawking fireworks, because here ħ is big. Nothing switches modes: the regimes *emerge* from particle count, density, and speed, the way they do in nature — just closer together.

**Status:** founded 2026-07-11 · spec-first (no tool source yet, by design) · current milestone **M0 `nexus`** — the CPU composition proof, contract under review.

## The two faces, one binary

`tinyuniverse.exe` is a windowed real-time app by default and a headless, contract-bounded instrument with flags (`--scenario … --steps … --seed … --json`), on the Buddhabrot v4 pattern. The seam between the simulation core (**libtinyverse** — deterministic, golden-gated, oracle-checked) and the renderer (a client, gated by `CINEMATIC.md`) is the **frame contract** (`contracts/frame.contract.md`). Visuals and timing are non-declared output; physics is sacred.

## Layout

| path | contents |
|---|---|
| `ARCHITECTURE.md` | the spec: dials, regime ladder, invariants, frame contract, milestones |
| `CLAUDE.md` | operational entry point for build sessions |
| `CINEMATIC.md` | the rendering law — no visual ships without this checklist green |
| `BUILD.md` | pinned toolchain and build commands |
| `DECISIONS.md` | ADRs (D-001…) |
| `TASKLIST.md` / `RUN_STATE.md` / `QUESTIONS.md` | plan, current state, open design questions |
| `contracts/` | sacred module/scenario contracts + schemas (contract-first) |
| `goldens/` · `harness/` | frozen behavior records · compile-as-verification CI (from M0) |
| `core/` · `app/` · `nexus/` | libtinyverse · windowed client · M0 composition proof *(created at build time)* |

## Quarries (verified on this machine)

| source | what Tiny Universe lifts |
|---|---|
| `C:\blackhole\files\blackhole.cu` (GARGANTUA) | Kerr null-geodesic renderer (BH-regime view + oracle), covariant redshift g, blackbody + post chain, cinematic/physical toggle |
| `C:\ASTRA-7\proto\astra_nexus.cpp` | the nexus pattern (regimes compose, PASS/FAIL battery), AstraCoord floating origin, rapidity numerics, Kepler-at-t_emit |
| `C:\Buddhabrot_CUDA` | one-binary-two-faces, two prioritized CUDA streams (presentation never waits), OptiX-via-stubs recipe, fat-binary self-contained dist |
| `C:\ORRERY\lib` (liborrery) | deterministic RNG / reductions / envelope / ckpt / regime — lifted under the verbatim rule |
| ORRERY `ratchet` receipts + THE UNFINISHED MIRROR §6 | the inscription (Ratchet) decoherence mechanic: collapse when record redundancy crosses `(1−p)ρ = p` |

ORRERY (`C:\ORRERY`) is the sibling instrument — same doctrine, different clock. Reference it; never edit it from here.

## Doctrine

**The spec is the product. The contract is sacred. The golden is load-bearing. The code is ephemeral.** Deterministic or it doesn't ship — same seed + same input trace ⇒ identical declared state. Every physics module names its oracle. Regimes emerge; nothing hard-switches. Beauty is a spec (`CINEMATIC.md`), not an accident.
