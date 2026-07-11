# RUN_STATE.md

**As of:** 2026-07-12 · **Milestone:** M6 `planck` (NOT STARTED — **operator-requested pause after M5**) · **State:** M5 CLOSED. **15/15 goldens GREEN in 75 s.**

## What M5 established

- **Black holes are real entities**: unscripted formation from the PM density argmax (collapse golden `5bcb5f58`, peak cell 182,720 vs M_FORM 10⁵), Paczyński–Wiita near-field + absorption with exact fixed-point ledgers, Hawking evaporation on the analytic cube-root clock (pop tick 2991 = expected, ledger closed to 4.3e-8), and screen-space lensing — the hole's own Hawking glow bends into an Einstein ring (`runs/gargantua_isco_lens.png`).
- **Measured physics** (D-017): the SR-inertia + P–W marginally-stable orbit lies in (3.3, 4.5)·GM/c² — declared as a band, nexus N4 pins the Newtonian limit at 3.0. A real seed-culling bug was caught by the collapse gate and the peak-cell probe.
- Goldens now 15: nexus `ad64f810` · kepler `2f93cdfb` · threebody `c2c572af` · cloud `8a0bc6e4` · galaxy `dcc87925` · merger `bfcfb003` · echo `40a84691` · ratchet `8f4b811e` · detector `83ea180d` · keprel `f985e473` · clocks `330c86a7` · photons `c4c565de` · collapse `5bcb5f58` · isco `5801ed2f` · hawking `6bd3faeb`.
- Contracts: nexus/frame/newton/arrow/gargantua v1.0.0, einstein v1.0.1. App v0.3.x: 14 scenarios, one binary two faces, `harness/verify.py` green.

## Next task when resumed — M6 `planck` (quantum bubbles)

Contract first (`contracts/planck.contract.md`), then:
1. Quantum bubbles: local comoving grids (64³ first), split-step cuFFT ψ evolution; spawn-on-isolation, collapse-on-inscription — **the Ratchet engine (M3, golden-verified) is the collapse mechanism, already waiting**; the detector scenario's record machinery binds to ψ.
2. Scenarios: `doubleslit` (build-your-own-detector — the M3 detector + fringes), `tunneling`, `sho-eigenstates`; oracles = nexus N5 (σ(t) 5.8e-15, E₀ 1.3e-13).
3. Measure the 16 GB bubble budget (Q-004 resolves by data).

**M6 gate:** double-slit fringes emerge from single collapses; nexus N5 parity in-sim; bubble budget measured; goldens + harness green.

## Chores carried

P1 Vulkan (SDK) · ImGui · TAA · clang/g++ nexus parity · art pass (incl. GARGANTUA Kerr geodesic view + OptiX spike → ORRERY `lens` report, D-017; collapse-scenario beauty pass) · cufft64 dll packaging · P³M/spatial hash · Q-006 derivation · 2.5PN inspiral.

## Standing context

Toolchain: CUDA 13.1, sm_89, MSVC 2022; build line in BUILD.md (envelope.cpp + cufft.lib). The universe now runs its full classical ladder: beauty → gravity → thermodynamics/inscription → relativity → black holes, every rung golden-gated. Repo docs authoritative over agent memories.
