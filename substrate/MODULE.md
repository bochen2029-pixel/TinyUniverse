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
