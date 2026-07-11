# PROPOSAL — TINY UNIVERSE v2: "SUBSTRATE"

**One field. One clock. One lattice.**

**Status:** PROPOSAL (operator review requested). Supersedes nothing — v1 (M0–M6, 18 goldens, HEAD `ca9e847`) remains the instrument of record and becomes v2's **oracle farm**. This document is the candidate `ARCHITECTURE.md` §-set for a v2 track (working name M8+ `substrate`, or a sibling repo `TINY UNIVERSE II`). Drafted 2026-07-12 by the M0–M6 build instance, from the operator's one-substrate brainstorm.

---

## 0 · Thesis (the sentence the whole proposal serves)

v1 proved the regimes; it did not unify them — the ψ engine runs *beside* the classical particles (~70% of the founding dream). The obstruction was never physics, it was compute: full GR + quantum field theory is uncomputable at nature's scales. **But the tiny universe doesn't run at nature's scales — and every dial chosen to make physics *visible* is exactly the dial that makes unified physics *computable*:**

| dial | visibility purpose (v1) | computability purpose (v2) |
|---|---|---|
| c = 20 su/s | relativity at walking speed | CFL: dt ≤ dx/c = 0.0125 s — the dial dt (1/240 s) fits with 3× headroom. **Numerical relativity becomes cheap.** |
| ħ = 0.5 | one particle is visibly wavy | λ_dB lands at grid scale — quantum pressure is resolvable on the same lattice as gravity |
| G = 2×10⁻³ | BHs at N ~ 10⁵ | horizons form *inside a 512-su box* at grid-resolvable radii — no AMR hierarchy needed to first order |

**The tiny universe is not a scale model of real physics. It is the scale at which real unified physics actually runs.** The gimmick is the method.

## 1 · The substrate, staged (each stage is a shippable milestone with its own oracles)

### S1 — Schrödinger–Poisson: one complex field replaces the particles
`iħ ∂ψ/∂t = [−(ħ²/2m)∇² + mΦ]ψ` · `∇²Φ = 4πG(m|ψ|² − ρ̄)` on a 3D lattice (256³ → 512³).

- **This is v1's two engines fused**: the split-step kinetic operator is M6's ψ engine; the FFT Poisson solve is M2's PM solver. Same cuFFT, one loop.
- **Regime emergence is literature-proven** (fuzzy dark matter: Schive, Mocz et al.): via the Madelung transform ψ = √ρ·e^{iS/ħ}, the field obeys Euler + a quantum-pressure term Q = −(ħ²/2m)∇²√ρ/√ρ that *switches itself off* where gradients are mild. Large scales: exact collisionless-matter behavior (galaxies, mergers, violent relaxation). Small scales: interference, tunneling, solitonic cores where quantum pressure halts collapse. **The quantum↔classical boundary is ħ/m — a dial, not a code path.** Nothing hard-switches, ever.
- Oracles: free-packet σ(t) and SHO (nexus N5 pedigree) · soliton mass–radius relation (analytic) · **v1 cross-checks**: galaxy rotation support, merger entropy rise, cloud collapse morphology reproduced *by the field*.
- Cost: ψ(2) + Φ(1) + scratch ≈ 8 floats/cell → 256³ ≈ 0.5 GB, **512³ ≈ 4.3 GB — feasible day one.**

### S2 — the lapse: proper time becomes a property of space
Promote v1's per-particle τ to a per-cell evolution rate α(x) = √(1 + 2Φ/c²): every cell advances its local phase by α·dt. Photons = a second massless field (or geodesic tracers) through the same α. Gravitational redshift and lensing stop being renderer features and become substrate behavior. Weak-field, declared (post-Newtonian field theory), one substrate.

### S3 — Einstein–Klein–Gordon: the real thing
The field coupled to actual spacetime evolution (BSSN/Z4c + complex scalar matter, moving-puncture gauge, tiny-c CFL). This is the *solved* literature class of scalar-field collapse — boson stars, oscillatons, **Choptuik critical collapse** — run at dials where it fits on one card.

- **The operator's dream sentence, verbatim, in one PDE system**: a ball of field-hydrogen contracts → pressure (S4 closure) supports it → fuel depletes → collapse resumes → **the lapse crashes toward zero and an apparent horizon forms** (the horizon finder replaces v1's density-argmax BH detector). Time near the hole slows because α does what Einstein's equations force it to do; light bends because it propagates through the actual metric. No patchwork.
- **Crown oracle: the Choptuik exponent.** A `substrate_nexus` CPU fp64 tool (nexus-II, contract-first, spherical 1D EKG) reproduces critical-collapse mass scaling M ∝ (p−p*)^γ, γ ≈ 0.374 — a number from the deepest end of gravitational physics, frozen as a golden, before any GPU code exists. Plus: oscillaton masses, Schwarzschild exterior, v1's isco/keprel/clocks/photons goldens as weak-field cross-checks.
- Cost: BSSN (~25 evolved fields) + scalar + radiation + RHS scratch ≈ 40–60 floats/cell with in-place discipline → 256³ ≈ 2.7–4 GB fp32; 384³ ≈ 9–13 GB; 512³ requires FP16 storage (§3).

### S4 — matter closures + the record layer (the Mirror ontology, load-bearing)
- **Fusion** = a *local source term in the same field equations*: amplitude → radiation-field energy above a density/temperature threshold, with fuel depletion (flux-limited diffusion for the radiation — another local PDE on the same lattice). Every sim below nuclear scale needs a closure; one substrate + local closures is still one substrate.
- **Hawking** = semiclassical horizon-flux law (declared, as in v1).
- **The Ratchet survives unification as the second layer**: an inscription lattice R(x) over the unitary substrate — records ratchet where the decoherence functional crosses threshold (v1's M3 machinery, generalized from particles to cells). v2's ontology is then explicitly the Mirror's: **the unitary field (the not-yet-written) + the record lattice (the written) + inscription between them.** Collapse remains a distinct stochastic process — that is a feature, and it is the design's deepest tie to the theory it grew beside.

## 2 · What does NOT unify — permanent, printed honesty

1. **Many-body entanglement.** A field on a 3D grid is a mean-field universe; N entangled particles live in 3N dimensions and no lattice escapes that exponential. Single-particle interference is exact (it *is* field interference); Bell/EPR correlations between separated lumps are faked or absent, forever, declared.
2. **Individual particle identity.** v2's "particles" are solitons/excitations of the field. Gameplay-visible dots become tracked excitation peaks (or v1-style tracer particles advected by the Madelung flow — a declared visualization layer, not dynamics).
3. Fusion/nuclear microphysics = closures (as in every stellar code ever written). Hawking = semiclassical. Measurement = Ratchet layer.

## 3 · The FluidX3D lessons (memory-bandwidth-first; the operator's instinct, formalized)

Substrate updates are **bandwidth-bound**, so the whole performance game is bytes/cell/step:

- **VRAM-resident everything** (v1 already lives this).
- **FP16 storage, fp32 arithmetic** (FluidX3D's FP16S/FP16C precedent) — halves bytes; constraint-drift impact measured per Invariant 6 (unbiased-acceleration rule), not assumed.
- **In-place streaming** (the "Esoteric Pull" idea the operator half-remembered — one lattice copy instead of two by interleaving reads/writes). For BSSN stencils: in-place update ordering or minimal double-buffering of only the fields that need it.
- **No atomics anywhere** — stencil gather only. Determinism gets *stronger* than v1's (already order-invariant) design, goldens get easier, and the Loschmidt echo becomes exact by construction (time reversal = conjugate ψ).
- **Arithmetic on this card** (4070 Ti SUPER, ~672 GB/s): at ~200 B/cell/step ≈ 3.4×10⁹ cell-updates/s → **256³ at ~200 steps/s — near-real-time at the 240 Hz dial**; 384³ at ~60 steps/s (¼-speed cosmology mode); 512³ with FP16. Sparse bricks (NanoVDB-style) when the universe is mostly empty.
- GPU is shared (v1 preflight discipline carries over verbatim).

## 4 · Relation to v1 (nothing is thrown away)

- **v1 is the oracle farm**: all 18 goldens become v2 cross-checks — rotation curves, the ISCO band, the echo, fringe contrast, the evaporation clock. v2 starts from eighteen receipts it must reproduce, not from trust. This is "the code is ephemeral, the golden is load-bearing" paying compound interest.
- **The renderer transfers whole**: the CINEMATIC stack reads density/flux fields instead of splatted particles (|ψ|² *is* an emission field); the lens warp is superseded by S2/S3 native bending; the HUD meters read the same fixed-point reductions.
- **The frame contract gets a MAJOR bump** (v2.0.0): field textures replace SoA particle buffers as the declared state; tracer particles become a visualization clause; goldens hash the field.
- Doctrine unchanged: contract-first, golden-gated, oracles named, errata logged, exit codes sacred.

## 5 · Proposed milestones & gates

| stage | deliverable | gate |
|---|---|---|
| **N0** `substrate_nexus` | CPU fp64 spherical-1D EKG tool (contract + battery: flat-space packet, oscillaton, subcritical dispersal, **Choptuik γ within 5%**) | battery green, golden frozen — *before any GPU substrate code* |
| **N1** `field` | SP 256³–512³ GPU; galaxy/cloud/merger/echo re-run *as field scenarios* | v1 cross-check gates + soliton oracle; echo exact by conjugation |
| **N2** `lapse` | α(x) evolution rates + photon field | clocks/photons v1 parity from the substrate |
| **N3** `horizon` | BSSN+scalar, collapse → apparent horizon | Schwarzschild exterior match; N0 Choptuik reproduced on GPU at declared tolerance |
| **N4** `star` | fusion closure + radiation + Ratchet lattice | **the sentence**: hydrogen ball → star → fuel-out → black hole, unscripted, one substrate, with the lapse and the light doing GR's work |

## 6 · Risks, priced

BSSN gauge/stability at toy dials is *plausible but unrehearsed* (mitigation: N0 CPU tool first, 1D spherical where everything is checkable) · fixed-grid horizon resolution bounds the smallest honest BH (declare it; AMR is a later contract) · FP16 constraint drift (measure, Invariant 6) · excitation-tracking for gameplay identity is genuinely open design · perf estimates are bandwidth arithmetic, not measurements (ADR-007: no claim survives contact with `--bench` unmeasured).

## 7 · The ask

Approve as the v2 track (M8+ or sibling repo). First concrete act on "go substrate": write `contracts/substrate_nexus.contract.md` and build N0 on the CPU — the Choptuik golden is the foundation stone, and it needs no GPU at all.

*v1 proved each rung of the ladder is real. v2 removes the ladder and lets the dials climb it alone.*
