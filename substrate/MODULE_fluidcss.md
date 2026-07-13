# MODULE — fluidcss (crown Stage A: radiation-fluid CSS background)

**Purpose.** Bank the **Evans–Coleman critical continuously-self-similar (CSS) radiation-fluid collapse background** — the fluid-β crown's **Stage A** — as a deterministic golden. This is the foundation stone of the critical-exponent measurement (β = 1/Re κ₀ = 0.35580192); Stage A is the background, Stage B (the perturbation eigenvalue that yields β) remains an **honest wall**.

**Contract.** [`contracts/fluidcss_nexus.contract.md`](../contracts/fluidcss_nexus.contract.md) v0.9.1.

**Tool.** `substrate/fluidcss_nexus.cpp` — single-file **CPU fp64** (no GPU). Deterministic C++ port of the overnight run's verified `hka_ec.py` (Hara–Koike–Adachi gr-qc/9607010; Evans–Coleman gr-qc/9402041). Fixed-step RK4 shoot from a regular center (ingoing sound cone, V<0) to the sonic point, tuning the central density `oi` (bisection) so the sonic velocity reaches V=−c_s. γ=4/3, c_s=1/√3, G=c=1.

**Invariants touched.** 1 · 3 (oracle: the analytic Evans–Coleman critical values) · 8 · 10 (fp64 = truth).

**Oracle.** The Evans–Coleman critical CSS solution's exact values (HKA 4.7–4.9): oi*=3/8, sonic (A0,N0,ω0,V0)=(3/2, 2/√3, 3/4, −1/√3), 2m/r=1/3, invariants N=N∞e⁻ˣ and A=1+⅔ω.

**Build.** `cl /std:c++17 /EHsc /O2 /W4 substrate\fluidcss_nexus.cpp /Fe:build\fluidcss_nexus.exe` (vcvars64; BUILD.md CPU path).

**Run.**
```
build\fluidcss_nexus.exe --stageA --json      # the Stage-A background object
build\fluidcss_nexus.exe --stageA --golden     # vs goldens/fluidcss_stageA/golden.hash
build\fluidcss_nexus.exe --selftest            # exact-sonic-state identities
```

**Golden (seed-free, deterministic).**

| golden | measured (EMERGES — nothing tuned) | gate |
|---|---|---|
| **fluidcss_stageA** `b4f4e463` | oi*=**0.3750013** (3/8) · sonic (1.50000, 1.15470, 0.75000, −0.57735) = (3/2,2/√3,3/4,−1/√3) to ~1e-5 · 2m/r=**0.333335** (1/3) · invariants max resid 7.7e-6 / 8.6e-8 · grid-converge 3.8e-6 | all under 2e-3 / 1e-4 / 5e-4 |

`--selftest`: the exact sonic state satisfies Dson=0 (1e-16), A=1+2ω/3 (exact), A_of match (exact).

**Determinism.** fp64 RK4 + fixed-iteration bisection, fixed grids → byte-identical declared JSON → blake2b golden. No RNG/GPU. Two-passed. Runs with the CPU oracles (GPU-independent).

**Honest boundary (D-027).** This banks **Stage A only**. **β = 1/Re κ₀ = 0.35580192 remains a WALL** — the Stage-B perturbation eigenvalue (the ∂ₛ-coupling of the perturbed fluid equations) fails its gauge-mode exactness gate (`tournament/gamma/phase4/RESULTS_hka_beta.md`); **β is NOT reported** (D-016/D-021 — a faked eigenvalue would poison the oracle farm). Supersedes the walled v0.9.0 4D reduction (which lacked a regular center on the ingoing sound cone and never reproduced oi*/the sonic point). Resume for β: re-transcribe HKA (5.5)–(5.10) from the primary TeX or re-derive the ∂ₛ-coupling, then the armed eigenvalue solver reads off the non-gauge zero.
