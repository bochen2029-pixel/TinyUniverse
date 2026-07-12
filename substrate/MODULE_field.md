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
field_nexus --scenario freepacket|sho3d|soliton [--grid N] [--seed N] [--json] [--golden] [--selftest]
```
`--golden` freezes/checks `goldens/field_<scenario>/golden.hash` (GOLDEN OK / NOT FROZEN / MISMATCH — tiny_nexus idiom). `--json` prints the declared envelope (hash domain: seed + params + result + gates + verdict; notes excluded). `--selftest` runs the flat-field Φ≡0 + free-evolution norm-conservation smoke. **Exit 0 pass · 1 declared gate fired · 2 error · 3 GPU contended (preflight, harness exit-3 lineage — NOT a red golden).** Each scenario pins its own grid/box for a stable golden; `--grid` is accepted but overridden by the scenario's declared geometry.

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
| **soliton** | 256³ / 64 su | SP scale-covariance r_c·M=const across M, 2M (in-tool κ) | *see FINDINGS below* | *pending* |

**freepacket + sho3d are the first GREEN goldens — the 3D split-step ψ engine validated exact against the fp64 N5 oracle at both ends the design needs.** They reproduce byte-identically on fresh runs (`--golden` re-runs the physics and confirms).

## Known issues / findings

- **soliton** is the load-bearing weld test (quantum pressure Q balancing self-gravity — the one scenario that only passes if the PM Poisson↔ψ fusion is correct; M6 never had gravity, M2 never had ħ). Regime & measured κ documented above / in DECISIONS. The half-peak-density radius r_c is the gated observable (less tail-biased than the half-mass radius, which is reported as a cross-check).
- **512³** stretch gated on a live VRAM preflight (exit-3-on-contention); 256³ is the safe gate.
- **Frame contract v2.0.0 MAJOR bump** (field textures replace SoA particle buffers) is *named* as N1's precondition but NOT executed here — N1 is headless-first (the renderer transfer is its own reviewed change).
- **Honest boundary (permanent):** ψ(x) on a 3-lattice is a mean-field, single-configuration universe. Single-particle interference is exact (it *is* field interference — the soliton, dispersal). Many-body entanglement (Bell/EPR between separated lumps) never unifies on a 3D grid — declared, forever. N1 makes no entanglement claim.
