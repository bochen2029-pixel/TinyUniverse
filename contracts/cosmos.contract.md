# CONTRACT — cosmos (M7 the tiny planet) · v1.0.0 · status: FROZEN 2026-07-12 (operator directive "go M7")

**Purpose.** Make the universe's topology real and visible: space becomes a live 3-torus (canonical coordinates, seam-aware physics), light's finite speed becomes *history you can see* (retarded-time rendering across box images), the namesake stereographic little-planet projection ships, an honest Einstein–de Sitter expansion mode runs on the comoving PM machinery, and the M6 quantum engine goes roaming — 3D bubbles that spawn on isolation and collapse by live Ratchet inscription.

## 1 · Topology goes live (declared; always-on)

- **Canonical torus coordinates:** after every drift, live-particle positions wrap into `[−L/2, L/2)` per axis (L = 512). Wrap arithmetic is `±512.0f` — a power of two, **exact in fp32** near the seam, hence bit-reversible (the Loschmidt echo survives). Kahan compensation is unaffected (exact-op shift). Max speed c = 20 su/s ⇒ ≤ 1/12 su/tick ⇒ at most one crossing per tick, enforced by construction.
- **Minimum-image convention** for all pairwise separations in the direct solver, BH force/absorption/dedup, and BH–BH forces: `d −= 512·round(d/512)` per axis. For every pre-M7 trajectory (|d| < 300) `round() = 0` and the arithmetic is bit-identical.
- **Exemptions (declared):** DEAD particles stay parked at x ≈ 10⁶ (outside canonical range; skipped by every kernel — unchanged M5 semantics). The fp64 tiny solver does not wrap; its scenarios are bounded ≪ L/2 by contract (max apoapsis 240).
- **Golden supersession:** the PM force field has been periodic since M2 — M7 only canonicalizes coordinates. Any prior golden whose particles cross a seam during its run re-freezes under this directive (empirically identified at implementation; recorded in D-019). Non-crossing goldens must reproduce **byte-identical** — any unexpected mismatch is run-state SUSPECT.
- Frame contract bumps to **v1.1.0** (MINOR): pos canonical-range clause; per-particle `bubble` handle buffer goes live.

## 2 · Declared engines

- **Expansion (scenario `expand`):** Einstein–de Sitter background on comoving coordinates. a(t) = (1 + (3/2)H₀t)^{2/3}, closed-form, a(0)=1; H₀ = √(8πG ρ̄₀/3) from the actual box mass (N=10⁶, m=1, L³ ⇒ H₀ = 1.11731×10⁻² s⁻¹). Momentum variable p = a²ẋ (m=1): kick ṗ = −∇Φc with the PM Green's function scaled ×1/a(t); drift ẋ = p/a²(t_mid). SR factors remain in the kernels but are inert (v_pec ~ 5×10⁻³ su/s, declared Newtonian-comoving); τ not gated in this mode. Zel'dovich growing-mode init: 100³ lattice, single mode k₁ = 2π·4/L along x, displacement amplitude A_d = 0.02/k₁ (δ₀ = 0.02), velocity v = H₀ψ. Growth observable: |ρ̂(k₁)| from a deposit+FFT measurement pass at checkpoints (every 2000 ticks), ratio-normalized (CIC window cancels).
- **3D quantum bubbles (scenario `bubbles` + `--bubbles` flag elsewhere, HUD-declared):** 64³ complex-float grid over a 16³ su window (dx = 0.25, M6 lineage), split-step at the dial dt, **free evolution v1** (V = 0; ⟨p⟩ is then conserved exactly — collapse re-localizes position only, momentum carries through; Φ-coupling is enumerated backlog). Absorbing cos² border (8 cells); norm decay through the border = declared window leakage (window fixed v1). NB_MAX = 16 slots, batched cuFFT.
  - **Spawn on isolation** (physics-triggered, Invariant 5): every 24 ticks, a massive, non-DEAD, non-bubbled particle whose 3³ PM-cell neighborhood holds no mass but its own (fixed-point deposit, quantum-exact compare) spawns a bubble: grid centered on the particle, ψ₀ = Gaussian σ₀ = 1 su with phase e^{i p·r/ħ}; regime |= QUANTUM; classical drift and τ suspend for the episode (declared τ-gap). Slot assignment in ascending particle id (deterministic).
  - **Collapse by inscription** (M3 semantics): a live bubble intruded by ≥1 foreign massive particle inside its window accrues one record per tick (strong write); at R ≥ 16 the position collapses: one counter-keyed draw `f(seed, owner, tick)` samples the 3D |ψ|² CDF (cell + uniform jitter), the owner resumes classical there, regime gains RECORDED|INSCRIBED, slot frees.
- **Light-history + little-planet (visuals, non-declared, CINEMATIC-gated):** ring buffer of past published positions (half4, snapshot every 24 ticks, depth adaptive to VRAM at alloc — never allocated in the headless JSON/golden path), sampled per periodic image at retarded time t − D/c (one fixed-point iteration, declared approximation); splats render the 27 nearest box images {−1,0,1}³ by minimum-image + offset. Stereographic little-planet projection as an alternate camera path (key L / `--planet`), same HDR→bloom→exposure→tonemap chain. HUD enumerates PROJ / HIST depth / A(t).

## 3 · Scenarios (seed 20260711)

| scenario | setup | steps | gates |
|---|---|---|---|
| **circumnav** | photon (\|p\|=1, +x) and 0.5c massive particle (u = 10/√0.75) launched from origin, no gravity (SOLV_NONE) | 12288 | photon laps twice, massive laps once, both return: \|Δx\|,\|Δy\|,\|Δz\| each < 0.02 su · τ_photon = 0 exactly · \|τ_massive/(T√0.75) − 1\| < 2×10⁻³ (fp32 accumulation floor, D-011 class) |
| **expand** | Zel'dovich lattice above, comoving PM | 24000 | LS slope of ln A vs ln a ∈ 1 ± 0.08 · \|A_end/A_0 / a_end − 1\| < 0.08 (linear growth D ∝ a; δ_end ≈ 0.039 stays linear) · p_drift < 10⁻³ · a_end declared in envelope |
| **bubbles** | 8 isolated loners (ring r = 150, p = 0.2) + 4 rigid intruder clusters (512 each) aimed to sweep bubbles 0–3 at t ≈ 8.7 s; deposit-only isolation checks (SOLV_NONE) | 4800 | n_spawned = 8 at first check · σ of the lowest-id untouched bubble at tick 1440 matches free-packet σ(t) = σ₀√(1+(ħt/2mσ₀²)²) to 2% (host-fp64 observable; 3.3σ window truncation ≈ 0.3%, floor-safe) · n_collapsed = 4, n_alive = 4 at end · declared hash includes bubble metadata block + surviving ψ grids |

Declared state hash (all three): particles + BH block (inactive) + bubble block where active. Oracle pedigree: circumnav is exact kinematics; expand's D ∝ a is the EdS linear-growth analytic; bubbles' σ(t) is nexus N5's free-packet law in 3D.

## 4 · Gate (M7, from TASKLIST)

The signature screenshot — the universe as a globe (little-planet), wrapped and lensed, self-visible via light-history — passes CINEMATIC §7 · circumnav/expand/bubbles goldens frozen · every pre-M7 golden either byte-identical or seam-crossing-superseded (D-019) · harness 21/21 GREEN.

## 5 · Changelog

- v1.0.0 (2026-07-12) — initial freeze under operator directive "go M7". GPU idle at authoring (14.4 GB free); wall-clocks citable.
