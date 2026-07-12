# MODULE — lapse (v2 N2: the clock)

**Purpose.** Give the substrate a clock. Turn N1's Newtonian potential Φ into a per-cell **lapse** α(x) = √(1 + 2Φ/c²) and a declared, hashed **proper-time field** τ(x) = ∫α dt. Demonstrates gravitational time dilation / redshift *on the substrate* — the bridge from M2's Newtonian potential to M4's relativistic time.

**Contract.** [`contracts/lapse.contract.md`](../contracts/lapse.contract.md) v1.0.0.

**Tool.** `substrate/lapse_nexus.cu` — single-file CUDA, `_nexus` family (self-contained blake2b golden hasher; envelope face; exit 0/1/2, GPU preflight 3). Reuses N1's PM cuFFT-Poisson (`kGreen`) **verbatim** for the `redshiftPM` weld.

**Invariants touched.** 1 (contract-first) · 3 (named oracle: analytic Schwarzschild) · 4 (fixed-point deposit, no float atomics — the PM path) · 6 (no fast-math) · 8 (headless envelope, exit codes) · 10 (fp32 core, fp64 host observables).

**Oracle.** Analytic **Schwarzschild gravitational redshift** z = 1/√(1 − r_s/r) − 1 (exact for the point-mass lapse, r_s = 2GM/c²). The `redshiftPM` weld additionally cross-checks the frozen N1/M2 PM-Poisson solver (well depth A = GM).

**Build.**
```
nvcc -O3 -arch=sm_89 -Xcompiler "/O2" -o build\lapse_nexus.exe substrate\lapse_nexus.cu cufft.lib
```
(vcvars64 first — see `BUILD.md`.)

**Run (headless).**
```
build\lapse_nexus.exe --scenario redshift   --json      # exact Schwarzschild time dilation
build\lapse_nexus.exe --scenario redshiftPM --json      # the substrate weld (PM well -> lapse)
build\lapse_nexus.exe --scenario NAME --golden          # vs goldens/lapse_NAME/golden.hash
build\lapse_nexus.exe --selftest                        # flatlapse: alpha=1, tau=N*dt exact
```

**Scenarios / goldens (seed 20260711, sm_89-pinned).**

| scenario | grid | proves | measured | golden |
|---|---|---|---|---|
| **redshift** | 128³ | the lapse+proper-time machinery reproduces **exact Schwarzschild** gravitational time dilation, weak→strong field | a clock at r≈2r_s ticks 40% slow (z=0.398); max \|α_meas/α_ana−1\| = **5.7e-6** (gate 1e-3); τ-integrator 6.6e-6 (gate 5e-5) | `e2c75be5` |
| **redshiftPM** | 128³ | **the weld** — the substrate's own PM-Poisson gravity, through the lapse, gives a well of Newtonian depth ⇒ correct gravitational redshift | fitted well depth A/GM = **0.9643** (rel **3.57e-2**, gate 5e-2 — the PM discretization floor, cf. N1 `poissontest` 3.6%); implied z(r=8)=0.066 | `3dddb950` |

`--selftest` (flatlapse): Φ=0 → α≡1 (max\|α−1\|=0), τ=N·dt (6.6e-7). PASS.

**Determinism.** Declared state = the τ field grid bytes → `state_b2b` → blake2b of the declared JSON = the golden. Two independent runs → identical bytes (verified: both goldens reproduce GOLDEN OK on fresh runs). Fixed-point deposit (Invariant 4); cuFFT bit-stable at fixed plan + sm_89.

**Dials.** c=20 su/s (**LIVE** — sets the lapse), G=2×10⁻³, m=1, dt=1/240, ħ=0.5 (inert). r_s = 2GM/c² = M·10⁻⁵ su.

**Honest boundary (printed in the contract).** N2 models only the **temporal** metric (g_tt): exact gravitational redshift/time-dilation, but **no light bending, no orbit precession, no Shapiro delay** (that is N3 `curve`). Static spacetime, no back-reaction; Schwarzschild-scalar (no frame-dragging). The lapse's exactness for a point mass is the g_tt = weak-field coincidence.

**Known issues / open questions.** Q-N2-1 (redshiftPM periodic-image background sets the 3.6% floor — the fit's C absorbs the constant part) · Q-N2-2 (fp32 τ accumulation floor ~7e-6, gated at 5e-5) · Q-N2-3 (innermost exact probe r≈2 su = 2r_s; closer to r_s the 1/√ nonlinearity × cell sampling would break the tight gate — not probed). All measured, none faked (D-016/D-021).

**Next rung.** N3 `curve` — the lattice metric back-reacts on energy density (light bending, precession); the N2 lapse becomes one component of a dynamical metric. The fluid-CSS Stage-A background plugs in as the static-collapse oracle.
