# RUN_STATE.md

**As of:** 2026-07-11 · **Milestone:** M7 `cosmos` **CLOSED** — the universe is a globe. **State:** THE REGIME LADDER + THE WORLD ARE BOTH COMPLETE. **21/21 goldens GREEN** (~105 s idle card). App v0.5.0 "cosmos": 20 scenarios. Contracts: + cosmos v1.0.0; frame → v1.1.0.

## What M7 established

- **Topology is live and exact.** Space is a 3-torus in coordinates, not just in the force field: `kDriftK` wraps positions into `[−256,256)` (the ±512 shift is a power-of-two, bit-reversible — the Loschmidt echo stayed byte-identical through it, D-019), minimum-image separations everywhere. `circumnav` `3f18f02c`: a photon laps the universe twice and comes home to 3.1e-5 su; photons are now **structurally ageless** (τ≡0) after an fp32 clamp-leak fix the gate caught.
- **Cosmology is a frozen number.** `expand` `ce448f2b`: Einstein–de Sitter comoving expansion on the existing PM solver (Green ×1/a, drift ×1/a², Zel'dovich growing mode) reproduces **linear growth D ∝ a** — slope 0.988, amplitude ratio 0.7% at a_end=1.927.
- **Quantum went roaming.** `bubbles` `78b753f1`: 64³ batched split-step bubbles spawn on PM-isolation (8/8), spread at the free-packet σ(t) to **0.57%**, and collapse by M3 Ratchet intrusion (4 collapsed / 4 alive). RECORDED particles stay classical.
- **The namesake shot shipped.** Stereographic little-planet projection + 27-image torus splat + light-history retarded-time sampling + bubbleVis ψ-glow, all non-declared, CINEMATIC §7 10/10. `runs/cosmos_littleplanet.png` (the globe), `runs/cosmos_planet_lensed.png` (a BH's Einstein ring on the little planet), `runs/cosmos_bubbles.png`. Perf: **225 fps avg / 58 min** windowed 1M with the 27-image splat + history.

## Next — operator's call (both doors open, one directive away)

- **v2 SUBSTRATE** (`docs/PROPOSAL_2026-07-12_v2_substrate.md`, operator-review pending): the one-field rewrite. First act on "go substrate": `contracts/substrate_nexus.contract.md` + build N0 (CPU fp64 spherical EKG, **Choptuik γ** golden — no GPU needed).
- **v1 polish backlog:** a(t) little-planet cosmology timelapse video · Φ-coupled bubbles · P1 Vulkan/ImGui · GARGANTUA per-pixel Kerr geodesic art pass · 2.5PN · Q-006.

## Chores carried (from M0–M6)

P1 Vulkan (SDK) · ImGui · TAA · clang/g++ nexus parity · art pass (Kerr geodesic view + OptiX spike → ORRERY `lens`) · cufft64 dll packaging · P³M/spatial hash · Q-006 · 2.5PN.

## Standing context

The regime ladder is COMPLETE (beauty → gravity → arrow → relativity → black holes → quantum) **and the world is complete** (wrapped, expanding, self-visible, roaming-quantum). Every physics claim a stranger reproduces cold from BUILD.md. 21 goldens, 8 contracts, ~105 s harness. GPU is shared (preflight exits 3 if <2 GB free; timings on a contended card not citable). Repo docs authoritative over agent memories.

## What M6 established

- **The ψ engine** (2D 256² split-step cuFFT, real + imaginary time, absorbing edges, host-fp64 observables) and **collapse as counter-keyed sampling** — the Tonomura experiment computed once, sampled 4096 times.
- **The measurement problem is a golden** (`doubleslit` `47a67d66`): fringe contrast 0.83 unobserved; the which-way detector (M3 record semantics) collapses it to 0.052. Interference emerges one particle at a time and dies when watched.
- **Oracle exactness**: tunneling T vs same-grid fp64 oracle to 1e-6 (`f1e7a061`); SHO ground state E₀ = ħω to 6 digits, σ exact (`fa2e009e`). D-018 findings: observable-vs-evolution coefficient bug caught by σ-exactness; oracle-isolates-implementation principle; near-field fringe shift is real physics.
- Q-004 resolved: ~5 MB per 64³ bubble — hundreds coexist with the full universe.
- Goldens now 18 (new: doubleslit `47a67d66` · tunneling `f1e7a061` · shoq `fa2e009e`). Contracts: + planck v1.0.0. App v0.4.x: 17 scenarios.

## Next task — M7 `cosmos` (the tiny planet)

Contract first (`contracts/cosmos.contract.md`), then:
1. **3-torus wrap live** (forces already periodic via PM — positions/light wrap next); light-history ring buffer ("see your own past" — ASTRA-7 Kepler-at-t_emit machinery, nexus N9-proven).
2. **Stereographic little-planet projection** — the product's namesake visual.
3. Scale factor a(t) cosmology mode (hot start → structure → heat death timelapse).
4. Roaming 3D bubbles integration (deferred from M6): spawn-on-isolation in the live universe, live Ratchet-collapse, bubbleVis.

**M7 gate:** the signature screenshot — the universe as a globe, wrapped and lensed, self-visible; cosmology timelapse runs; goldens + harness green.

## Chores carried

P1 Vulkan (SDK) · ImGui · TAA · clang/g++ nexus parity · art pass (Kerr geodesic view + OptiX spike → ORRERY `lens`; doubleslit fringe display saturation; collapse-scenario framing) · cufft64 dll packaging · P³M/spatial hash · Q-006 · 2.5PN · live bubble-budget probe.

## Standing context

The regime ladder is COMPLETE: beauty → gravity → thermodynamics/inscription → relativity → black holes → **quantum mechanics** — every rung golden-gated, every physics claim a stranger can reproduce cold from BUILD.md. 18 goldens, 7 contracts, ~199 s harness. GPU is shared (preflight in verify.py; timings on a contended card are not citable). Repo docs authoritative over agent memories.
