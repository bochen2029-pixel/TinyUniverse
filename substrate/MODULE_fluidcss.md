# MODULE — fluidcss (crown, COMPLETE: radiation-fluid CSS critical exponent β)

**Purpose.** Measure the **critical exponent of radiation-fluid gravitational collapse** — β = 1/Re κ₀ — from the true **Evans–Coleman** critical CSS background (Stage A) and its relevant perturbation eigenvalue (Stage B). **COMPLETE (D-032, 2026-07-19): κ₀ = 2.810577211 → β = 0.355798800** (literature 2.8105525488 / 0.35580192; |Δβ| = 3.1e-6 ≪ the 4e-3 G-ANCHOR gate). Two goldens: `fluidcss_stageA` `27af7920` (background) + `fluidcss_stageB` `9f8587fd` (β), both two-passed.

**Contract.** [`contracts/fluidcss_nexus.contract.md`](../contracts/fluidcss_nexus.contract.md) v1.0.0 (FROZEN — the original gates fire unchanged).

**Tool.** `substrate/fluidcss_nexus.cpp` — single-file **CPU fp64** (no GPU), γ=4/3, G=c=1. **Stage A** = HKA's own construction (gr-qc/9607010 §IV): sonic values 4.7–4.9 at free V₀; sonic→center RK4 shoot; deterministic bisection on the case1(A<1)/case2(second-sonic) transition → V₀\*=0.1124394014, N̄'(sonic)=−0.355699 (the paper's fingerprint), exactly one V-zero, center-relaunch closure 4.2e-6. **Stage B** = the §V shoot: sonic Frobenius modes (ker R + identity + gauge, order-18) → sonic→center eigen-ODE on the stored background → κ-bisections with TWO analytic controls (sonic-gauge 0.355740 vs analytic 0.355699; origin-gauge 0.999954 vs exact 1). **History (D-032):** the v0.9.x "Stage A" was the collapsing flat **Friedmann** solution mislabeled as EC (its sonic criterion V=−c_s IS the Friedmann point of the sonic line; measured V=−√(1−1/A) to 1.9e-10); the superseded golden `b4f4e463` is documented in `goldens/fluidcss_stageA/NOTE.md`; the Friedmann build is retained as the `--friedmann` control face. Python cross-validation at higher precision: `tournament/gamma/phase4/nr_ec2.py` (background) + `nr_lyap.py` (Lyapunov spectrum, κ=2.8105526) + `nr_shoot_ec.py` (shoot, κ=2.810552374).

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
