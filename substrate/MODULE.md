# MODULE — substrate_nexus (v2 N0)

**Purpose.** The foundation stone of TINY UNIVERSE v2 "SUBSTRATE": a CPU fp64 spherically-symmetric Einstein–Klein–Gordon solver (polar-areal gauge, constrained evolution) — the standing oracle for the v2 GPU ladder (N1 field → N2 lapse → N3 horizon → N4 star), built **before any GPU substrate code** and needing no GPU (runs under any card contention).

**Contract:** `contracts/substrate_nexus.contract.md` v1.0.0.
**Units:** geometric G = c = 1 (the Choptuik phenomenon is dimensionless/universal — D-020).
**Oracle spine:** flat-space wave conservation · subcritical dispersal · supercritical horizon formation · the Type-II critical transition · massive-KG stability. Golden `13aa73e5` (N=800, r_max=24).

## Build (BUILD.md; no GPU)

```
cl /std:c++17 /EHsc /O2 /W4 substrate\substrate_nexus.cpp /Fe:build\substrate_nexus.exe
```

## Usage

```
substrate_nexus [--json] [--selftest] [--golden] [--only S1..S5] [--dev P] [--N n] [--rmax R]
```
`--golden` freezes/checks `goldens/substrate_nexus/golden.hash` (GOLDEN OK / NOT FROZEN / MISMATCH, tiny_nexus idiom). `--dev P` runs a single evolution at amplitude P and prints diagnostics (physics dev). `--json` prints the declared envelope (hash domain: params + per-test metrics + verdict; notes excluded).

## Formulation (polar-areal EKG, constrained evolution)

`ds² = −α²dt² + a²dr² + r²dΩ²`, real scalar φ, `V = ½μ²φ²`, first-order `Φ≡∂_rφ`, `Π≡(a/α)∂_tφ`.
- Evolution (method of lines, RK4, KO dissipation ε=0.5, origin parity ghosts, outgoing outer BC): `∂_tφ=(α/a)Π`, `∂_tΦ=∂_r[(α/a)Π]`, `∂_tΠ=(1/r²)∂_r[r²(α/a)Φ]−αaμ²φ`.
- Constraints (a, α by outward integration of the LOGS each RK substage — guarantees positivity through the stiff near-horizon region, D-021): `(ln a)_{,r}=(1−a²)/(2r)+2πr(Π²+Φ²)+4πr a²V`; `(ln α)_{,r}=(a²−1)/(2r)+2πr(Π²+Φ²)−4πr a²V`, α rescaled to `αa→1` at r_max.
- Mass aspect `m=(r/2)(1−1/a²)`; horizon at `2m/r→1` (polar slicing → lapse collapse `α(0)→0`).

## Battery (S1–S5 + S6 determinism)

| id | what | key result (golden config) |
|---|---|---|
| S1 | flat-space wave (ε pulse) | ADM mass conserved 7.9e-4; disperses; weak field |
| S2 | subcritical dispersal (0.6 p*) | no horizon (2m/r=0.23); lapse recovers 0.53; disperses |
| S3 | supercritical collapse (1.5 p*) | horizon forms 2m/r=0.98; M_BH=0.43; lapse collapses 3.9e-5 |
| S4 | **Choptuik Type-II transition** | p*=0.404 bracketed; M_BH 0.412→0.118 → 0 at p* (ratio 0.29); clean power law R²=0.998; effective exponent 0.24 (**grid-limited diagnostic; precise γ=0.374 deferred to AMR — D-021**) |
| S5 | massive-KG stability (μ²=1) | stable + mass conserved 3.3e-6; mass term correct (bound oscillaton deferred) |
| S6 | determinism | in-process S1 identical; out-of-process proven by GOLDEN OK on re-run |

## Known issues / deferred (D-021)

- **Precise Choptuik γ ≈ 0.374** needs adaptive mesh refinement (uniform grid caps the self-similar curvature; refining turns near-critical chaotic). Deferred to the AMR contract / possibly the GPU `horizon` (N3) stage. N0 proves the *transition*, not the exponent.
- **Bound oscillaton** (massive-field eigen-profile) deferred; S5 checks stability + conservation only.
- **clang++/g++ parity** build owed (same as nexus — MSVC-pinned golden for now).

---

# MODULE — field_nexus (v2 N1: the Schrödinger–Poisson weld)

**Purpose.** The first GPU substrate rung: M6's split-step ψ welded to M2's PM cuFFT-Poisson (`kGreen` verbatim) in one loop — quantum pressure vs self-gravity, on a 3D lattice (D-022).
**Contract:** `contracts/field.contract.md` v1.0.2 (FROZEN). **Oracle:** nexus N5 (freepacket/sho3d exact) + analytic SP soliton scale-covariance (3e-8) + v1 free-fall/two-body cross-checks + structural echo.
**Goldens (6/6, GPU, behind the harness preflight):** `field_freepacket` `03dd3a3b` · `field_sho3d` `dfbc6185` · `field_soliton` `d163d765` (the weld; ~7 min) · `field_echoF` `433ddcc8` · `field_cloudF` `2308ea49` · `field_mergerF` `a09dda6a`.

```
nvcc -O3 -arch=sm_89 -Xcompiler "/O2" -o build\field_nexus.exe substrate\field_nexus.cu cufft.lib
```

**Deferred (declared):** galaxyF rotation curve (needs quantized vortices + 512-su box — Q-N1-1).

# MODULE — lapse_nexus (v2 N2: the clock)

**Purpose.** The substrate grows a clock: per-cell lapse α(x)=√(1+2Φ/c²) + a declared, hashed proper-time field τ(x)=∫α dt (c goes LIVE). Exact Schwarzschild gravitational redshift; the PM weld makes the substrate's own gravity dilate time (D-023).
**Contract:** `contracts/lapse.contract.md` v1.0.0 (FROZEN). **Oracle:** analytic z = 1/√(1−r_s/r)−1 (α-err 5.7e-6) + the N1/M2 PM solver (Newtonian-depth well, 3.6%).
**Goldens (2/2, GPU):** `lapse_redshift` `e2c75be5` · `lapse_redshiftPM` `3dddb950`. `--selftest` = flatlapse (Φ=0 ⇒ α≡1 exact).

```
nvcc -O3 -arch=sm_89 -Xcompiler "/O2" -o build\lapse_nexus.exe substrate\lapse_nexus.cu cufft.lib
```

**Honest boundary:** temporal metric only — bent paths are N3 `curve`.

# MODULE — curve_nexus (v2 N3: geometry curves)

**Purpose.** Geodesics through the substrate's weak-field metric ds²=−(1+2Φ/c²)c²dt²+(1−2Φ/c²)dl²: light bends at exact GR 4GM/bc² (the 1919 factor of 2 decomposed — lapse half + space half), orbits precess at 6πGM/(c²a(1−e²)), and the Shapiro delay lands — with N2, ALL FOUR classical tests of GR pass on the substrate (D-024).
**Contract:** `contracts/curve.contract.md` v1.1.0 (FROZEN). **Oracle:** analytic GR closed forms (deflect 0.3–0.75% · precess 0.52% · shapiro 0.33%). CPU fp64 — no GPU, runs under any contention.
**Goldens (3/3):** `curve_deflect` `4e6c33ca` · `curve_precess` `67272705` · `curve_shapiro` `20bfd4d2`. `--selftest` = flat GM=0 (straight ray + closed orbit).

```
cl /std:c++17 /EHsc /O2 /W4 /nologo substrate\curve_nexus.cpp /Fe:build\curve_nexus.exe /Fo:build\curve_nexus.obj
```

**Honest boundary:** static weak-field geodesic oracle — a dynamical GPU metric field with back-reaction is N4+.

# MODULE — fluidcss_nexus (Axis-C crown: radiation-fluid CSS, Stage-A)

**Purpose.** The Evans–Coleman radiation-fluid CSS critical background in deterministic C++ fp64 (RK4 center→sonic shoot on the central density oi): oi\*=0.3750013 (exact 3/8), sonic point (3/2, 2/√3, 3/4, −1/√3) to ~1e-5, Misner–Sharp 2m/r=1/3 — everything EMERGES, nothing tuned (D-027).
**Contract:** `contracts/fluidcss_nexus.contract.md` v0.9.1 (Stage-A golden FROZEN; Stage-B scoped — do not change gates). **Oracle:** HKA gr-qc/9607010 §IV exact values + invariants A=1+(2−γ)ω (7.7e-6), N·eˣ const (8.6e-8); grid-convergent 3.8e-6.
**Golden:** `fluidcss_stageA` `b4f4e463` (`--stageA --golden`).

```
cl /std:c++17 /EHsc /O2 /W4 /nologo substrate\fluidcss_nexus.cpp /Fe:build\fluidcss_nexus.exe /Fo:build\fluidcss_nexus.obj
```

**Honest wall (Stage-B):** β=1/κ₀=0.35580192 NOT yet measured, none faked (D-016/D-021). The perturbation operator is PROVEN correct from first principles (D-031; `tournament/gamma/phase4/nr_rederive.py`); the sonic→center shoot surfaces only the gauge mode κ=1 — the wall is the eigenvalue-EXTRACTION method. Next: HKA's Lyapunov time-evolution (§V.G) → κ → β → Stage-B port → v1.0.0. The β research thread lives in `tournament/gamma/phase4/` (`RESULTS_hka_beta.md`).
