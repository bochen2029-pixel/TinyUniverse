# RUN_STATE.md

**As of:** 2026-07-12 · **Milestone:** v2 **N3 `curve` — geometry curves: light bends + orbits precess at exact GR, on the substrate** (3/3 goldens GREEN — all four classical tests of GR, on master). **State:** v1 complete (M0–M7); v2 ladder N0 (oracle) + N1 `field` (SP weld, 6/6) + N2 `lapse` (the clock, 2/2) + **N3 `curve` (geometry, 3/3) all CLOSED on master**. N3 integrates geodesics through the substrate's weak-field metric ds²=−(1+2Φ/c²)c²dt²+(1−2Φ/c²)dl²: `curve_deflect` `4e6c33ca` (**light bends at exact GR 4GM/bc²** — the 1919 factor of 2 DECOMPOSED: lapse/time [N2] = one half, spatial curvature [N3] = the other; full/lapse=2.004) + `curve_precess` `67272705` (**perihelion precession** 0.72°/orbit, 0.52%) + `curve_shapiro` `20bfd4d2` (**Shapiro time delay** 0.33%). **With N2's redshift, the substrate now passes ALL FOUR classical tests of GR from its own metric** (redshift · light-bending · precession · Shapiro delay). CPU fp64 geodesic oracle (no GPU, like N0). Tool `substrate/curve_nexus.cpp`; contract v1.1.0; D-024. Honest boundary: static weak-field geodesic oracle — a dynamical GPU metric field with back-reaction is N4+. **Next = operator's call: N4 `horizon` per the ROADMAP is the Choptuik critical-exponent crown, which is the documented AMR research wall (D-021, and the overnight fluid-β/scalar-γ walls) — NOT a clean one-pass rung. The weak-field GR arc (N1 field → N2 lapse → N3 curve) is COMPLETE; genuine next options: the AMR contract (unlocks the crown), the renderer (Axis B), v1 polish, or the fluid-β/scalar-γ crown research.** — *(prior N2 milestone line, superseded:)* the clock: gravitational time dilation, on the substrate (2/2 goldens GREEN, on master). **State:** v1 complete (M0–M7); v2 N0 CLOSED; **N1 `field` CLOSED** (merged to master `8ace261`); **N2 `lapse` CLOSED** — the substrate now has a clock. N2 turns N1's Newtonian potential Φ into a per-cell lapse α(x)=√(1+2Φ/c²) and a declared, hashed proper-time field τ(x)=∫α dt (c goes LIVE). Goldens: `lapse_redshift` `e2c75be5` (**EXACT Schwarzschild redshift** — a clock at r≈2r_s ticks 40% slow, z=0.398, α-err 5.7e-6 vs 1/√(1−r_s/r)−1) + `lapse_redshiftPM` `3dddb950` (**the weld** — the substrate's own PM-Poisson gravity makes a well of Newtonian depth 3.6% → correct gravitational time dilation). Both two-passed cold. Tool `substrate/lapse_nexus.cu`; contract `lapse.contract.md` v1.0.0; D-023. Honest boundary: temporal metric only (redshift exact; light-bending/precession → N3). *Reconciliation note: RUN_STATE/TASKLIST had lagged git — N1 was merged to master at `8ace261`; both are now corrected to ground truth.*

## OPERATOR WORK QUEUE (set 2026-07-12, in priority order)

1. **#3 · v1 polish** ← **IN PROGRESS.** a(t) little-planet timelapse video · GARGANTUA Kerr geodesic art pass · **2.5PN gravitational-wave inspiral** ✅ **DONE** (`nexus/inspiral_nexus.cpp`, Peters 1964 — circular merger time to 1.3e-13, eccentric a(e) circularization to 5e-11; goldens `2eba79de`/`4578d3ac`; D-025) · Q-006 ✅ **DONE** (`nexus/precession_nexus.cpp`, D-026 — the **7π superposition is CORRECT**; the app's 6.41π was a normalization artifact, nominal vs force-distorted p; goldens `a0e180df`/`db0818f2`/`f9df648f`). **#3 remaining: a(t) timelapse video · Kerr art pass.**
2. **#4 · crown research** — fluid-β Stage-B (the HKA perturbation ∂ₛ coupling) · scalar-γ = 0.374 (DSS, G97 gr-qc/9604019).
3. **#2 · renderer (Axis B)** — R0 `interop` Vulkan⇄CUDA live window (NOTE: Vulkan SDK not installed — measured gate).
4. **#1 · AMR contract** — the named future work that unlocks the Choptuik critical exponent (β/γ) for N0/N4.

---

**Prior (master):** v2 **N0 `substrate_nexus` CLOSED** — the substrate ladder's foundation stone. v1 complete (M0–M7); v2 track begun (D-020). **22/22 goldens GREEN** (~121 s): nexus + substrate_nexus (CPU fp64, no GPU) + 20 GPU. Harness runs CPU oracles under any card contention.

## What N0 established (v2 foundation stone; D-020/D-021)

- **A spherical Einstein–Klein–Gordon solver** (`substrate/substrate_nexus.cpp`, CPU fp64, polar-areal constrained evolution, G=c=1) — the standing oracle for the v2 GPU ladder, built **before any GPU code** and needing no card. Golden `13aa73e5` (N=800, r_max=24, ~20 s).
- **Four rigorous oracles** a stranger runs cold: flat-space mass conservation (7.9e-4) · subcritical dispersal · supercritical horizon formation (2m/r→0.98, lapse→3.9e-5) · massive-KG stability+conservation (3.3e-6).
- **The Choptuik phenomenon, honestly**: S4 demonstrates the **Type-II critical transition** — a bracketed critical amplitude p*=0.404 where the black-hole mass falls continuously to zero (M_BH 0.412→0.118, ratio 0.29; arbitrarily small black holes), clean far-field power law R²=0.998. **The precise universal exponent γ ≈ 0.374 is DEFERRED to the AMR contract (D-021)**: a uniform grid caps the self-similar curvature (measured), refining turns near-critical chaotic — Choptuik needed AMR. Per RAYFORMER/D-016, N0 ships what it proves (the transition) and names what it cannot (the exponent) rather than fake a number that would poison the oracle farm.

## Next — the v2 GPU ladder, or v1 polish (operator's call)

- **v2 continues:** N1 `field` (Schrödinger–Poisson 256³–512³ GPU, the two v1 engines fused) → N2 `lapse` → N3 `horizon` (BSSN+scalar; may reach precise γ) → N4 `star` (the hydrogen-ball sentence). Each gated against N0 + the 21 v1 goldens (the oracle farm). Contract-first per module.
- **v1 polish backlog:** a(t) little-planet cosmology timelapse video · Φ-coupled bubbles · P1 Vulkan/ImGui · GARGANTUA Kerr geodesic art pass · 2.5PN · Q-006.
- **AMR contract** (would let N0/N3 reach the precise Choptuik γ) — named future work.

## Chores carried

clang/g++ parity (nexus + substrate_nexus) · P1 Vulkan · ImGui · TAA · art pass (Kerr view + OptiX → ORRERY `lens`) · cufft64 dll packaging · P³M/spatial hash · Q-006 · 2.5PN · bound oscillaton (N0 S5 deferred).

## Standing context

v1: the regime ladder + the world, complete (M0–M7, 21 goldens). v2: the substrate rewrite begun — N0 the fp64 oracle spine (22nd golden). Every physics claim a stranger reproduces cold from BUILD.md; every claim the tool can't back is named, not faked (D-016/D-021). GPU shared (preflight exits 3 <2 GB; CPU oracles run regardless). Repo docs authoritative over agent memories.

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
