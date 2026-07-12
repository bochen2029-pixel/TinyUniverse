# MODULE — curve (v2 N3: geometry curves)

**Purpose.** Give the substrate a **spatial** metric. N2 curved time (the lapse → redshift); N3 curves space. Geodesics through the weak-field metric `ds² = −(1+2Φ/c²)c²dt² + (1−2Φ/c²)dl²` bend light and precess orbits at **exact GR** — and the famous 1919 factor of 2 is **decomposed** into the N2 (time) half and the N3 (space) half.

**Contract.** [`contracts/curve.contract.md`](../contracts/curve.contract.md) v1.0.0.

**Tool.** `substrate/curve_nexus.cpp` — single-file **CPU fp64** (no GPU; runs under any card contention, like N0 `substrate_nexus`). Geodesic ODE integration (RK4): the eikonal ray for light deflection, the relativistic orbit equation for precession.

**Invariants touched.** 1 (contract-first) · 3 (named oracle: analytic GR) · 8 (headless envelope, exit codes) · 10 (fp64 = truth).

**Oracle.** Analytic GR: light deflection **4GM/bc²** (Eddington 1919) + the 2GM/bc² lapse-half; perihelion precession **6πGM/(c²a(1−e²))** (Einstein 1915).

**Build.**
```
cl /std:c++17 /EHsc /O2 /W4 substrate\curve_nexus.cpp /Fe:build\curve_nexus.exe
```
(vcvars64 first — BUILD.md CPU path.)

**Run (headless).**
```
build\curve_nexus.exe --scenario deflect --json     # light bending, the factor-2 decomposition
build\curve_nexus.exe --scenario precess --json     # perihelion precession
build\curve_nexus.exe --scenario NAME --golden       # vs goldens/curve_NAME/golden.hash
build\curve_nexus.exe --selftest                     # flat GM=0: straight ray, closed orbit
```

**Scenarios / goldens (seed 20260711).**

| scenario | proves | measured | golden |
|---|---|---|---|
| **deflect** | light bends at exact GR 4GM/bc²; **space curvature doubles** the lapse-only bending | b=100/200: full ×1.008/1.003 of 4GM/bc²; lapse ×1.004/1.001 of 2GM/bc²; **full/lapse = 2.008/2.004** (GR: exactly 2) | `4e6c33ca` |
| **precess** | orbits precess at the GR rate | a=1000, e=0.5: 0.72°/orbit vs 6πGM/(c²a(1−e²)), rel **5.2e-3** (weak-field, leading 1PN) | `67272705` |

`--selftest` (flat, GM=0): straight ray (deflection 0), closed orbit (precession 1e-11). PASS.

**Determinism.** fp64 RK4, fixed step counts → byte-identical declared JSON → blake2b golden. No RNG, no GPU, no atomics. Both goldens two-passed (GOLDEN OK on independent re-run). Runs with the CPU oracles in the harness (GPU-independent, no preflight).

**Honest boundary (printed in the contract).** N3 is the **geodesic/curvature oracle** — it proves the static weak-field metric (sourced by the substrate's Poisson Φ; g_tt lattice-validated by N2 `redshiftPM`, 3.6%) bends light and precesses orbits at exact GR. A **live GPU metric field a²(x) with dynamical back-reaction** (the energy density sourcing the metric each step) is the harder continuation toward N4/horizon — named, not faked. Weak-field/1PN only; strong-field is N4+.

**Next rung.** N4 `star` — fusion closure + radiation + the Ratchet lattice (the hydrogen-ball sentence), or the horizon/BSSN direction that may reach the precise Choptuik γ.
