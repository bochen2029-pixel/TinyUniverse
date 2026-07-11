# CONTRACT — planck (M6 quantum bubbles) · v1.0.0 · status: FROZEN 2026-07-12 (operator directive "go M6")

**Purpose.** Real single-particle quantum mechanics in the universe: a split-step ψ engine with genuine interference, tunneling, and eigenstates — and **collapse by inscription**: the M3 Ratchet/detector machinery (golden-verified) is the measurement mechanism. The double slit becomes an experiment you run, one particle at a time.

## Declared engine (v1 scope)

- **ψ grid:** 2D 256×256 complex float over a 64×64 su window (dx = 0.25 su), cuFFT C2C, Strang split-step at the dial dt (real time) or τ = 0.002 (imaginary time). Absorbing edges (cos² border, 8 cells, declared). Potentials/masks: hard wall-with-slits (perfectly absorbing except apertures — Tonomura-style, declared), square barrier V₀, harmonic well.
- **Observables host-side fp64** from the downloaded grid (norm via GPU fixed-point reduce for imaginary-time renormalization — Invariant 4 shapes).
- **Collapse = position measurement at the screen window** (x ∈ screen ± 2 su at t*, y-marginal — declared approximation): each of N independent identical particles samples the *same* final |ψ|² via counter-keyed inverse-CDF draws (`f(seed, particle, dim)`) — the Tonomura electron-at-a-time experiment, computed once, sampled N times. Dots become real particles at the screen plane (p = 0, declared re-pointification).
- **Which-way mode:** a detector at the slits inscribes at slit-passage (M3 slab semantics: deterministic strong write, R → 64, regime 0x80|0x40 latched); each particle's branch is drawn from the branch weights (∫|ψ_slit|²), then evolves *that* masked branch to the screen. Interference must die.
- **v1 boundary (declared):** single-bubble engine scenarios, headless-first (windowed shows the final dot field statically); roaming 3D 64³ bubbles with spawn-on-isolation in the live universe are the M7+ integration. **Q-004 budget resolution:** 64³ complex float = 2.1 MB/bubble + scratch+plan ≈ ~5 MB/bubble → ~3.2 GB per thousand bubbles: hundreds of bubbles fit alongside the M2–M5 universe (~0.8 GB) with room to spare on 16 GB. Declared resolved by arithmetic; a live batched-FFT probe rides the M7 integration (GPU contended at freeze time — noted).
- **Declared state hash:** particles (dots) + BH block (inactive) + the final ψ grid bytes.

## Scenarios (seed 20260711; engine runs at init; golden `steps` = 0 sim ticks, engine steps declared below)

| scenario | engine | gates |
|---|---|---|
| **doubleslit** | packet (σ=1.5, p=2 → λ_dB = π/2 su) from x=−15 through a wall at x=0 (2 slits: width 1.5, separation d=4), screen at L=14, 3600 steps; 4096 dots sampled; both modes run: no-detector + which-way | fringe contrast **C(k*) = 2|Σe^(−ik*yⱼ)|/N at the analytic k* = d·p/(ħL)**: no-detector C_ψ > 0.35 and C_dots > 0.30 · which-way C_det < 0.12 · k-scan peak within 5% of k* · n_recorded = 4096 in detector mode · |C_dots − C_ψ| < 0.1 (sampling agrees with theory) |
| **tunneling** | packet (σ=2, p=1.5, E≈1.125) vs square barrier V₀=1.8, width 1 su, 2400 steps; T = transmitted fraction at t* | |T_2D − T_ref| < 1e-3 where **T_ref = in-scenario host fp64 1D split-step oracle**, same params (y-separable physics) · sanity 0.005 < T < 0.2 |
| **shoq** | 2D isotropic SHO ω=1, imaginary time τ=0.002 × 15000 | |E − ħω|/ħω < 1e-3 (2D ground = ħω) · |σ_x − √(ħ/2mω)|/σ < 1e-2 |

Oracle pedigree: nexus N5 proved the 1D split-step math in fp64 (σ(t) to 5.8e-15, E₀ to 1.3e-13); the tunneling oracle is computed fresh per run (host fp64); doubleslit's k* is analytic.

## Gate (M6, from TASKLIST)

Double-slit fringes emerge from single collapses and **die when the detector watches** · tunneling matches its fp64 oracle · SHO ground state at ħω · bubble budget declared (Q-004) · 18/18 goldens green.

## Changelog
- v1.0.0 (2026-07-12) — initial freeze under operator directive. GPU contended at authoring (100% util, 3 GB free): engine loads are trivial (256² FFTs); goldens remain valid under contention (order-invariant reductions, deterministic cuFFT); wall-clocks not cited.
