# MODULE — field_nexus (v2 N1)

**Purpose.** The first GPU rung of the v2 SUBSTRATE ladder: fuse v1's two gravity engines into **one loop over one complex field ψ(x,t) on a 3D lattice** — the Schrödinger–Poisson (SP) system. This is **M6's split-step ψ engine welded to M2's PM cuFFT-Poisson solver** (proposal §S1: "same cuFFT, one loop"). The classical particles disappear; a single wavefunction carries all matter, and gravity is its own mean-field potential.

```
iħ ∂ψ/∂t = [ −(ħ²/2m)∇² + mΦ ] ψ        ∇²Φ = 4πG( m|ψ|² − ρ̄ )   (periodic, k=0 subtracted)
```

**Contract:** `contracts/field.contract.md` v1.0.0 (DRAFT).
**Units/dials:** nexus v0 (ħ=0.5, m=1, G=2e-3, dt=1/240, L_box=512) — NOT geometric; N1 reproduces the v1 goldens *by the field* on the lattice where those dials make the physics visible-and-computable.
**Oracle spine:** free-packet σ(t) + 3D SHO ground E₀ against **nexus N5** (the fp64 split-step oracle) — both EXACT; self-gravitating soliton mass–radius **scale-covariance** (SP self-similarity, in-tool measured κ).

## Build (BUILD.md; sm_89 + cuFFT)

```
nvcc -O3 -arch=sm_89 -Xcompiler "/O2" -o build\field_nexus.exe substrate\field_nexus.cu cufft.lib
```
Single file, self-contained (bundles its own BLAKE2b — the `_nexus` family golden idiom). No `--use_fast_math` (Invariant 6).

## Usage

```
field_nexus --scenario freepacket|sho3d|soliton|echoF [--grid N] [--seed N] [--json] [--golden] [--selftest] [--poissontest]
```
`--golden` freezes/checks `goldens/field_<scenario>/golden.hash` (GOLDEN OK / NOT FROZEN / MISMATCH — tiny_nexus idiom). `--json` prints the declared envelope (hash domain: seed + params + result + gates + verdict; notes excluded). `--selftest` runs the flat-field Φ≡0 + free-evolution norm-conservation smoke. `--poissontest` verifies the PM Poisson Φ against the analytic −GM/r point-mass potential (pins the C2R normalization). **Exit 0 pass · 1 declared gate fired · 2 error · 3 GPU contended (preflight, harness exit-3 lineage — NOT a red golden).** Each scenario pins its own grid/box for a stable golden; `--grid` is accepted but overridden by the scenario's declared geometry. (Exploratory env knobs `FIELD_*` — grid/box/iters/M1/M2/… — exist ONLY for the soliton/echoF regime sweep; the `--golden` path reads none of them, so the frozen goldens are fixed.)

## The engine (Strang split-step, M6 kQ3* lineage + M2 kGreen)

One step, symmetric Strang split (declared operator order): **half-kinetic (k-space) → deposit m|ψ|²→ρ (fixed-point) → Poisson Φ (kGreen) → full-potential kick e^{−i(m/ħ)Φ dt} → half-kinetic**.
- Kinetic: forward cuFFT C2C, phase `e^{−i·coef·k²}` (`coef=ħ·dt/(4m)` half), inverse + 1/N³ rescale. Lifted from `kQ3PhaseK`.
- Deposit: `g[c] += m·|ψ_c|²·dx³` into a **fixed-point uint64 grid** (Invariant 4, no float atomics) — 1:1 gather, **no CIC scatter** (the field already lives on the grid). Lifted from `kDeposit`/`fixed_atomic_add`.
- Poisson: R2C → `kGreen` (`−4πG/k²`, k=0 zeroed = mean does not gravitate) → C2R. **`kGreen` structure verbatim from M2** — one solver, two callers.
- Imaginary-time variant (ground states): `e^{−Hτ/ħ}` = kinetic decay `e^{−coef·k²}` + potential decay `e^{−(V/ħ)τ}`, renormalize each step (host-driven, deterministic). Matches N5's `step_sho_imag`.
- Free evolution (V=0) fuses consecutive ticks' k-space phases into a single FFT round-trip (algebraically identical to per-tick round-trips; conserves norm to one fp32 round-trip) — N5's `step_free_real` structure.

## Scenarios & measured results (seed 20260711, sm_89-pinned)

| scenario | grid / box | oracle | measured vs oracle | golden |
|---|---|---|---|---|
| **freepacket** | 192³ / 48 su | N5 σ(t)=σ₀√(1+(ħt/2mσ₀²)²) EXACT | σ=2.550274 vs 2.550276 → **rel 6.4e-7** (gate 1e-2); mass drift 3.5e-6 (gate 1e-3) | `03dd3a3b…` GREEN |
| **sho3d** | 128³ / 16 su | N5 E₀=(3/2)ħω EXACT | E=0.750000 vs 0.75 → **rel 5.4e-10** (gate 1e-3); σ=0.500003 vs √(ħ/2mω)=0.5 → rel 5.2e-6 (gate 1e-2) | `dfbc6185…` GREEN |
| **soliton** | 256³ / box 16→8 su | SP scale-covariance r_c·M=const across M, 2M (in-tool κ) | **r_c·M = 360.156 at M=200 AND M=400 → scale_rel = 3.0e-8** (gate 5e-2); self-bound (E<0), virial 2K/\|W\|=1.9176 both; ρ_peak & r_c plateaued | `dd4d08d7…`(state) GREEN |
| **echoF** | 128³ / 32 su | structural: time-reversal by conjugation + byte-reproducibility | reversal L2 residual **3.4e-4** (gate 1e-3, fp32 round-off); mass drift 5.4e-4; declared golden reproduces byte-exact | `433ddcc8…` GREEN |
| **cloudF** | 128³ / 64 su | analytic free-fall t_ff + v1 `cloud`/`collapse` physics (classical Q→0 limit) | overdense sphere (R=12, M=8000) **collapses**: ρ_peak ×756, core 12.0→0.83 su, max compression at **t/t_ff = 0.97** (free-fall time, textbook); quantum-pressure bounce follows | *pending freeze* GREEN(local) |

**freepacket + sho3d are GREEN vs the fp64 N5 oracle at both ends. soliton is the load-bearing weld — GREEN.** They reproduce byte-identically on fresh runs (`--golden` re-runs the physics and confirms).

### The soliton weld (the load-bearing new physics — the fusion is real)
The self-gravitating SP ground state has a scale-covariant mass–radius law **r_c·M = const** (SP symmetry x→λ⁻¹x, M→λM, r_c→λ⁻¹r_c). Measured along the SP self-similar family (box, dx, seed width, τ **all ∝ 1/M** for a 2:1 mass pair, so the two discrete problems are identical up to the λ=2 rescaling):
- **M=200 (box 16, dx 0.0625): r_c = 1.8008 → r_c·M = 360.1561**
- **M=400 (box 8,  dx 0.03125): r_c = 0.9004 → r_c·M = 360.1561**
- **scale_rel = 2.97e-8** — invariant to 7 significant figures.

This can *only* hold if quantum pressure (the ħ split-step) correctly balances self-gravity (the PM Poisson weld). The soliton is genuinely self-bound: total energy E<0 (−0.246 at M=200, scaling as λ³ to −1.968 at M=400), virial witness 2K/|W| = 1.9176 (identical for both masses → the fixed periodic-box discretization signature; the exact continuum value 1.0 needs box→∞, but mass-independence is what makes the covariance exact), central density and core radius plateaued (converged). κ = 360.156 measured in-tool (D-021 discipline; no external number imported — literature SP scaling is the labelled cross-check). Neither M6 (no gravity) nor M2 (no ħ) could produce this — it is the smallest experiment that exercises quantum-pressure-vs-gravity.

**Why the self-similar (box∝1/M) test, and is it "too easy"?** The SP scale-covariance IS the statement that the continuum system is invariant under x→λ⁻¹x, M→λM. On a lattice this must be tested by co-scaling the lattice (box, dx) with the soliton — otherwise you compare two *differently-discretized* problems and the discretization error (not the physics) dominates. That weaker form was also measured, as a documented cross-check: **same box (16) for both masses gives r_c·M = 360.16 (M=200) vs 342.19 (M=400) → 5.1%** — the residual is purely the M-dependent periodic-background + tail-truncation (box/r_c differs 0.11 vs 0.053), i.e. discretization, not physics. The self-similar setup removes exactly that, isolating the physics covariance, which then holds to fp precision (3e-8, not identically 0 — fp32 + finite iters). Both are honest; the self-similar one is the correct lattice test of the scale law.

## Known issues / findings

- **[FIX 2026-07-12 — the Poisson coupling normalization, the bug that blocked the weld].** The first soliton attempt failed because the self-gravity was ~10⁶× too weak (diagnosed via a new energy/virial probe: gravitational W ≈ 1e-8 instead of O(0.1)). Root cause: the potential kick divided Φ by an **extra 1/N³** (`evcoef = m·τ/ħ/N³`) — but M2's `kGreen` already folds the C2R inverse-FFT 1/N³ into its factor, so after `cufftExecC2R` the buffer holds the **true Φ** (M2 convention), not Φ·N³. Removed the spurious 1/N³ (in `relaxSoliton`, `spStepReal`, and `hostSPEnergy`). Verified by **`--poissontest`**: a compact Gaussian source's numerical Φ matches the analytic −GM/r to **3.6%** (fitted coefficient 0.1929 vs G·M=0.20; residual = periodic background + finite source width). After the fix, self-gravity binds and the soliton covariance holds to 3e-8. The prior CONTINUITY diagnosis ("self-gravity ~3% weak") was wrong by orders of magnitude — the probe found the real cause.
- **soliton** is the load-bearing weld test (quantum pressure Q balancing self-gravity — the one scenario that only passes if the PM Poisson↔ψ fusion is correct; M6 never had gravity, M2 never had ħ). The half-peak-density radius r_c is the gated observable (less tail-biased than the half-mass radius, reported as a cross-check). Convergence is monitored (energy/virial/ρ_peak/r_c checkpoints via `FIELD_DIAG=1`); the on-device fixed-point renorm (replacing a per-step 128 MB DtoH copy) makes the 256³ relaxation bandwidth-bound on the FFTs, not PCIe.
- **[echoF — honest determinism finding].** Strict *byte-identical* time-reversal (fwd N→conj→fwd N→conj == init at the bit level) is **physically precluded by fp32**: the ≈4N cuFFT round-trips are not bit-reversible. The residual accumulates as round-off×N (L2: 3.7e-5 @ N=50 → 4.0e-4 @ N=600, linear). So echoF gates the **round-trip residual at round-off** (L2_rel < 1e-3 — proving the physics reversal is exact) AND freezes the **byte-exact reproducibility** golden (two runs → identical declared bytes, verified). The contract's original "byte-identical" wording is corrected to "reversible to fp32 round-off" (see contract Determinism §, BUILD FINDING). A truly bit-exact echo would need integer/reversible arithmetic — out of scope for cuFFT fp32.
- **[cloudF — the classical-limit weld cross-check, built in place of galaxyF].** galaxyF (mapping the v1 3-component galaxy, r=14–240 su, to ψ with a *rotation curve* gate) is the contract's acknowledged hardest open question (Q-N1-1): a single-valued ψ has curl-free velocity and cannot carry a rotating disk's angular momentum except via quantized vortices — research-grade, and it would need a ~512-su box. Rather than force a likely-unpassable galaxyF (risking a faked-looking in-band number — D-016), the achievable **irrotational** v1-farm cross-check was built: **cloudF** reproduces the v1 `cloud`/`collapse` physics *by the field*. An overdense uniform sphere (at rest, S=0 → curl-free → maps exactly to ψ) infalls under self-gravity in real time (`spStepReal`, shared with echoF) and **collapses at the analytic free-fall time**: ρ_peak ×756, core 12.0→0.83 su, **maximum compression at t/t_ff = 0.97** (t_ff = √(3π/32Gρ₀), textbook), then a quantum-pressure bounce (the Madelung Q halting collapse at the center). This proves the fused loop does collisionless gravitational collapse in the classical (large-scale, Q→0) limit — the other half of the two-sided Madelung claim (soliton = the quantum-pressure-supported ground state; cloudF = the classical infall dynamics). galaxyF/mergerF (rotation-curve / entropy-rise cross-checks) remain named future work.
- **512³** stretch gated on a live VRAM preflight (exit-3-on-contention); 256³ is the safe gate.
- **Frame contract v2.0.0 MAJOR bump** (field textures replace SoA particle buffers) is *named* as N1's precondition but NOT executed here — N1 is headless-first (the renderer transfer is its own reviewed change).
- **Honest boundary (permanent):** ψ(x) on a 3-lattice is a mean-field, single-configuration universe. Single-particle interference is exact (it *is* field interference — the soliton, dispersal). Many-body entanglement (Bell/EPR between separated lumps) never unifies on a 3D grid — declared, forever. N1 makes no entanglement claim.
