# N1 `field` — implementation plan (terse; design-only companion to `contracts/field.contract.md`)

**Status:** PLAN. Contract is DRAFT, unreviewed. No `.cu` written. This is the map a build session follows once "go N1" lands.
**One line:** fuse M2's PM cuFFT-Poisson + M6's split-step ψ into one Schrödinger–Poisson loop on a 256³ (→512³) lattice. Same cuFFT, one loop (PROPOSAL §S1).

## The loop (one step; Strang, normative order)
```
half-kinetic (k-space)  →  deposit |ψ|²→ρ (fixed-point)  →  Poisson Φ (kGreen)  →  full-potential kick e^{-i(m/ħ)Φ dt}  →  half-kinetic
```
Conjugation (ψ.y ← −ψ.y) between forward runs = exact time reversal → the echo golden.

## Kernels: lift / adapt / new
| kernel | source | change |
|---|---|---|
| `kPhaseK_field` | `kQ3PhaseK` (`app/tinyuniverse.cu:960`) | 3D whole-grid kinetic half; coef = ħ·dt/(4m). Near-verbatim (drop the per-bubble slot indexing). |
| `kPsiDeposit` | `kDeposit` (`:374`) idiom | **new but same fixed-point path**: cell c → `fixed_atomic_add(&g[c], m·(ψ.x²+ψ.y²)·dx³)`. 1:1 gather, **no CIC scatter** (field is on-grid). Invariant 4. |
| `kGreen` | `:396` | **VERBATIM** (−4πG/k², k=0 zeroed, same normalization). R2C→multiply→C2R. The whole point: one solver. |
| `kPhaseV_field` | `kQPhaseV` (`:915`) | potential half: `e^{-i·(m dt/ħ)·Φ}·ψ`, V=mΦ from the *live* Poisson solve (not a static array — the only real new wiring). |
| `kSponge` | `kQ3Edge` (`:976`) | absorbing BC only (freepacket); cos² 8-cell mask ×ψ. **Not** in periodic/echo scenarios. |
| `kNormAcc` | `kQ3NormAcc` (`:987`) | fixed-point norm for imaginary-time renorm (soliton/sho3d) + mass-conservation meter. |
| `kConj` | new (1 line) | ψ.y ← −ψ.y. echoF + any reversal receipt. |
| `kFixToReal` | `:392` | VERBATIM (uint64 ρ → float ρ for the FFT). |

Plans: one C2C (kinetic), one R2C+C2R pair (Poisson). Fixed for the run → cuFFT bit-stable on sm_89.

## First golden to target
**soliton** (256³, periodic, imaginary-time relax to SP ground state, gravity on). Why first: it is the **smallest scenario that can ONLY pass if the fusion is correct** — quantum pressure Q must balance self-gravity at the mass–radius scaling `r_c·M = κ`. M6 never had gravity; M2 never had ħ. If soliton passes, the weld holds. Then: freepacket/sho3d (N5 exactness, cheap) → galaxyF/mergerF (v1-farm cross-checks) → echoF (the byte-identical receipt).

## Order of build
1. Smoke: flat ψ → Φ≡0 to the ulp (proves `kGreen` reuse is honest); norm-conservation `--selftest`.
2. `kPhaseK_field` + FFT round-trip (free-packet dispersion vs N5 σ(t) — no gravity yet).
3. Add deposit→`kGreen`→`kPhaseV_field` = the full loop; soliton relax + hold.
4. Conjugation → echoF byte-identical (the determinism proof).
5. v1 cross-checks (galaxyF/mergerF) — outcome bands, not byte-match.
6. 512³ variant behind a GPU preflight (exit-3-on-contention, `d162ae4`).

## Risks (priced; all → measure, don't pre-tune — D-018)
- **R1 IC→ψ mapping** (rotational support has no exact wavefunction; phase is curl-free). Galaxy rotation is the exposed gate — measure vs the 8% band, widen only if the residual forces it. *(Q-N1-1)*
- **R2 imaginary-time renorm determinism at 256³** — M6 proved 256²; the fixed-point norm divide is unverified at 3D. Check in step 1. *(Q-N1-3)*
- **R3 512³ VRAM under a shared card** — 4.3 GB + two plans + scratch vs a contended 16 GB. Preflight-size the plans first; 256³ is the safe gate. *(Q-N1-4)*
- **R4 soliton κ provenance** — measure in-tool (relax 2 masses, read slope), literature value as labelled cross-check (D-021 `gamma_eff` discipline). *(Q-N1-2)*
- **R5 frame contract v2.0.0** — field textures replace SoA particles = MAJOR bump; a separate reviewed change, NOT done inside N1 (N1 is headless-first, names the bump as a precondition).

## The honest boundary (carried from the contract, PROPOSAL §2)
Many-body entanglement never unifies on a 3D grid — mean-field, forever. Single-particle interference exact; Bell/EPR faked or absent, declared. N1 = non-relativistic SP; no lapse/horizon/fusion (those are N2/N3/N4).

## Doctrine touchpoints
Invariant 1 (contract before source) · 3 (oracle named: N5 + analytic + v1 farm; SUSPECT on fail) · 4 (fixed-point deposit, no float atomics) · 5 (regimes emerge — ħ/m is a dial, **no per-cell regime bit**) · 8 (headless full citizen, exit 0/1/2) · 9 (HUD declares the boundary). Exit codes sacred. `[ARGUMENT-GRADE]` perf until `--bench` (ADR-007).
