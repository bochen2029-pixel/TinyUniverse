# CONTRACT — lapse (v2 N2: the clock) · v1.0.0 · status: BUILDING 2026-07-12

**Purpose.** The second rung of the v2 SUBSTRATE ladder (`ROADMAP.md` §3, N2): give the substrate a **clock**. N0 gave the deterministic shell + a critical threshold; N1 `field` gave one self-gravitating field ψ and its Newtonian potential Φ. N2 turns that potential into **proper time**: a per-cell **lapse** α(x) sets how fast each cell's clock runs, and a declared, hashed **proper-time field** τ(x) integrates it. This is the bridge from M2's Newtonian potential to M4's relativistic time — *inside* the substrate:

```
α(x)      = √( 1 + 2Φ(x)/c² )              (the lapse: static-observer clock rate)
dτ(x)/dt  = α(x)                            (per-cell proper time; τ = ∫ α dt)
∇²Φ       = 4πG( ρ − ρ̄ )                    (Poisson — the SAME PM solve as N1/M2, verbatim)
```

**The physics claim (the emergence):** *clocks run slow in gravity wells, on the substrate.* A signal emitted at radius r from a lattice mass and observed far away is redshifted by **z = 1/√(1 − r_s/r) − 1** (r_s = 2GM/c²). This is Schwarzschild gravitational time dilation, produced by the substrate's own gravity — not a hard-coded formula.

**Why the oracle is EXACT, not merely weak-field.** For a point mass Φ = −GM/r, the lapse `α = √(1 + 2Φ/c²) = √(1 − 2GM/rc²) = √(1 − r_s/r)` is *algebraically identical* to the exact Schwarzschild metric's static-observer lapse (g_tt = −(1 − r_s/r)c²). So the redshift oracle z(r) = 1/√(1 − r_s/r) − 1 holds **exactly** for all r > r_s, spanning the weak field into the strongly-dilated regime — even though N2 models only the **temporal** part of the metric (see Honest boundary). The scale dials make it visible: with c = 20, G = 2×10⁻³, r_s = 2GM/c² = M·10⁻⁵ su, so a mass of ~10⁵ gives r_s ≈ 1 su and clocks at a few su tick tens of percent slow.

Single-file CUDA oracle tool `substrate/lapse_nexus.cu`, self-contained in the `_nexus` family idiom (blake2b golden hasher lifted verbatim; envelope face `--scenario X --json|--golden|--selftest`; exit 0/1/2, GPU preflight 3). It reuses N1's PM-Poisson (`kGreen`) **verbatim** — one solver, now a third caller (D-005 spirit: don't fork what is frozen).

## Units & dials (v1 dials — c is now LIVE)

N2 inherits the frozen nexus v0 dials and introduces **no new fundamental constants**. Unlike N1 (where c was declared inert — SP is non-relativistic), **N2 makes c live**: it is the constant that converts the Newtonian potential into a clock rate.

| dial | value | role in N2 |
|---|---|---|
| c | **20 su/s** | **LIVE** — sets the lapse α = √(1 + 2Φ/c²); r_s = 2GM/c². The whole rung is about c re-entering. |
| G | 2×10⁻³ su³/(sm·s²) | sets the well depth Φ and hence r_s = 2GM/c² (M·10⁻⁵ su) |
| m (m_p) | 1 sm | the mass quantum for the PM deposit (redshiftPM source) |
| ħ | 0.5 sm·su²/s | **inert in N2** (no ψ evolution; the source is a static density). Declared, unused. |
| dt | 1/240 s | fixed timestep; τ integrates α over N ticks (the per-cell clock) |
| L_box | 512 su | the periodic box for the PM Poisson solve (redshiftPM); analytic scenarios pick their own box |

Not geometric units (unlike N0's G=c=1): N2 lives in the v1 dial system so the lapse, r_s, and redshift are the *same numbers* the v1 relativity module (M4 `einstein`) and cosmos clocks used — the oracle farm N2 must stay consistent with.

## Declared engine (the lapse + proper-time integrator)

**Grid.** Cubic lattice `LN³` over a cubic box side `gL`; `dx = gL/LN`. Real fields Φ, α, τ (fp32). For the PM scenario, the N1 Poisson working set (fixed-point ρ grid + R2C/C2R plans) is allocated; analytic scenarios need only Φ/α/τ.

**Φ source — two declared modes (per scenario):**
- **ANALYTIC** (`redshift`, `flatlapse`): Φ filled directly on the grid. Point mass: Φ(x) = −G·M/√(r² + ε²), ε a tiny declared core-softening (ε ≪ smallest probed radius; ε = 0.01 su — affects only the unprobed center cell, not the gate). `flatlapse`: Φ ≡ 0.
- **PM POISSON** (`redshiftPM`): a resolved mass on the lattice → the N1 weld path **verbatim**: `kZeroFix → kPsiDeposit(|ψ|²→ρ) → kFixToReal → R2C → kGreen(−4πG/k², k=0 zeroed) → C2R` → Φ (true, M2 convention — the same solver N1's soliton/cloudF/mergerF use). The source ψ = √ρ is a Gaussian blob scaled to total mass M.

**Lapse.** `kLapse`: α[i] = √( max(1 + 2·Φ[i]/c², α_floor²) ), α_floor a small positive clamp (1e-4) so a cell that reaches the horizon (1 + 2Φ/c² ≤ 0) yields a finite tiny α rather than NaN — declared; the physics gates probe only r > r_s where the argument is safely positive. Φ < 0 in wells ⇒ α < 1 ⇒ the clock runs slow. Flat space Φ = 0 ⇒ α ≡ 1.

**Proper-time integration.** `kAccumTau`: τ[i] += α[i]·dt, iterated over N ticks (the declared step count). Because N2 is a **static** spacetime (no velocity field, static mass), the integral of the constant lapse is τ(x) = α(x)·N·dt — but it is computed as a genuine per-cell per-tick accumulation, so (a) it is the literal "per-cell proper-time field driving the clock" the ladder calls for, deterministic over steps, and (b) it is forward-compatible with N3, where a dynamical metric makes τ genuinely path-dependent. **τ is the new declared field** (the golden hashes it).

**Determinism.** Every kernel is a per-cell map (fill Φ, lapse, τ-accumulate) plus the N1 Poisson path (fixed-point deposit — Invariant 4, no float atomics; cuFFT bit-stable at fixed plan + sm_89). No scatter, nothing to race. **Golden = (scenario, dials, seed, grid, steps) → blake2b-256 of the declared JSON** (which carries `state_b2b` = blake2b of the τ field grid). sm_89-pinned. Two independent runs ⇒ identical declared bytes.

## Scenarios & gates (seed 20260711; the golden units)

| scenario | grid / box / Φ-source | what runs | gate | oracle |
|---|---|---|---|---|
| **redshift** | 128³, box 64, ANALYTIC Φ = −GM/r (M = 10⁵ → r_s = 1 su, ε = 0.01) | fill Φ → α → integrate τ over N = 480 ticks; sample the shell-mean clock rate α(r) = τ(r)/(N·dt) at radii r ∈ {2, 4, 8, 16} (r_s/r = 0.5 … 0.0625, weak→strong); form redshift ratios vs the outermost shell | **max over probed radii** of `|(1+z)_meas / (1+z)_analytic − 1| < 1e-3`, where (1+z)_meas = α(r_out)/α(r) and (1+z)_analytic = √((1−r_s/r_out)/(1−r_s/r)); + α(flat cell)→ within [floor,1]; NaN-free | **analytic Schwarzschild** z = 1/√(1−r_s/r) − 1 — EXACT (see Purpose). The lapse+proper-time machinery vs the closed form. |
| **redshiftPM** | 128³, box 64, PM POISSON of a Gaussian blob (σ = 4 su resolved; M = 10⁵) | N1 weld path → Φ_PM → α → τ; invert the lapse to Φ_meas(r) = ½c²(α²−1); radial-fit Φ_meas = −A/r + C over the far field r ∈ [8,24] (r ≫ σ, r ≪ box/2 — C absorbs the periodic-image background) | **`|A/(G·M) − 1| < 0.05`** (the substrate's clock field encodes a gravitational well of Newtonian depth GM ⇒ the redshift magnitude is right; A = GM ⇔ z(r)=1/√(1−2A/rc²)−1 correct) + implied z at r = 8 reported + mass-conservation + NaN-free | **the substrate weld** — PM Poisson (N1/M2) + the analytic −GM/r depth (cf. N1 `--poissontest`, 3.6%). The emergence claim: substrate gravity → correct gravitational redshift. |

**`--selftest` (flatlapse):** Φ ≡ 0 → α ≡ 1 (max|α−1| < 1e-6) and τ = N·dt exactly (|τ/(N·dt) − 1| < 1e-6). The no-gravity sanity + a determinism smoke. Mirrors N1's flat-field selftest.

**The gate is on the derived clock physics, not on IC bytes.** `redshift` is exact against the closed form (the machinery oracle, like N1's freepacket/sho3d vs N5); `redshiftPM` is the honest substrate weld with a measured lattice tolerance (like N1's soliton/poissontest). No scenario is gated against an unmeasured number (D-016/D-021).

## Declared state & the frame contract

- **Declared state = the τ field grid bytes** (`float[LN³]`), hashed blake2b-256 into `state_b2b`, canonical byte order (x-fastest, matching M2/M6/N1). The declared JSON (tool, seed, params, result{state_b2b, observables}, gates, verdict) is hashed → the golden.
- **Frame contract:** N2 adds the **lapse α and proper-time τ as declared field textures** to the substrate state (they already exist as named buffers in `frame.contract.md` v1 — `tau` per-particle — now promoted to fields). This is additive to the N1 v2.0.0 field-state; **N2 names, but does not itself execute, the frame-contract note** (headless-first, like N1). The renderer reading τ to drive animation rate is a non-declared visualization clause (ARCHITECTURE §7).

## The headless face (envelope-conformant)

```
lapse_nexus.exe --scenario redshift|redshiftPM [--grid N] [--seed N=20260711] --json   # declared JSON envelope
lapse_nexus.exe --scenario NAME --golden                                                # frozen params vs goldens/lapse_NAME/golden.hash
lapse_nexus.exe --selftest                                                              # flatlapse: alpha=1, tau=N*dt
```

Declared body (canonical order): `seed, params{scenario, grid, steps, c, G, m, hbar, dt, L_box, gL, M, r_s}, result{state_b2b, <z_r, z_analytic, redshift_rel | A_fit, A_expect, A_rel, z_at8>, tau_ref, alpha_min, nan_free}, gates{primary, norm, nan}, verdict`. Floats fmt6/fmt9. **Exit 0 pass · 1 declared gate fired · 2 error · 3 GPU contended** — never conflated (Invariant 8).

## Oracle pedigree (Invariant 3 — no silent fallback)

- **Exact analytic:** `redshift` — Schwarzschild gravitational redshift z = 1/√(1−r_s/r) − 1, exact for the point-mass lapse. Oracle-gate failure ⇒ run-state SUSPECT.
- **The substrate weld:** `redshiftPM` — the PM-Poisson potential's Newtonian depth (A = GM) read through the lapse; cross-checks the frozen N1/M2 solver (a free structural guard that fusing the lapse didn't perturb the potential).
- **The N0/N1 tie:** N2 does not gate against N0 (EKG, geometric units) — the tie lands at N3 `curve`/`horizon` where the metric back-reacts. N2's oracle is the analytic Schwarzschild lapse + the N1 solver. Stated so the pedigree is honest about which stone carries which rung.

## The honest boundary (permanent, printed)

1. **N2 models only the TEMPORAL part of the metric.** The lapse α = √(1 + 2Φ/c²) gives exact gravitational time dilation / redshift (g_tt), but N2 does **not** model the spatial curvature (the (1 − 2Φ/c²) on dx²): **no light bending, no orbit precession, no Shapiro delay** — those are N3 `curve`. A photon's *frequency* shift is right; its *path* is not yet bent.
2. **Static spacetime, no back-reaction.** The mass/potential is static; the lapse does not feed back into the field dynamics, and the field's energy does not source the metric beyond Newtonian Φ. Dynamical geometry is N3+.
3. **No frame-dragging, no rotation (Kerr) — scalar redshift only.** The v1 `gargantua` Kerr view remains the rendering-side BH; N2 is the substrate's own gravitational clock, Schwarzschild-static.
4. **Weak-field Poisson source.** Φ is the Newtonian potential (∇²Φ = 4πGρ); the lapse's exactness for a *point* mass is the g_tt coincidence, not a full strong-field solve. Honest for the redshift observable; the strong-field metric is N3/N4.

## Build runbook (design → code → golden → harness → register)

1. **Contract (this file).** ← *you are here.* (Operator directive "execute against ROADMAP immediately" = authorization to proceed; Invariant 1 satisfied by writing the contract first.)
2. **MODULE.md** — `substrate/MODULE_lapse.md`: purpose, this contract, invariants touched (1,3,4,8), oracle (analytic Schwarzschild + N1 solver), build command, known issues.
3. **Implement** `substrate/lapse_nexus.cu` — lift from `field_nexus.cu`: blake2b, Grid+gridAlloc, kGreen + the Poisson path, fmt6/9, gpuPreflight, the envelope/golden/main scaffold. Add: `kFillPointPhi`, `kLapse`, `kAccumTau`, radial shell-mean + fit host readouts.
4. **Build** — `nvcc -O3 -arch=sm_89 -Xcompiler "/O2" -o build\lapse_nexus.exe substrate\lapse_nexus.cu cufft.lib` (BUILD.md; no fast-math — Invariant 6).
5. **Golden freeze** — `goldens/lapse_redshift/golden.hash`, `goldens/lapse_redshiftPM/golden.hash` (run `--golden` → `GOLDEN NOT FROZEN <hash>` exit 2 → write the hash → rerun `GOLDEN OK`). sm_89-pinned.
6. **Harness row** — add a `LAPSE_GOLDENS` block + the build line to `harness/verify.py` (behind the GPU preflight).
7. **Register** — ARCHITECTURE §11 v2 N2 row, TASKLIST, DECISIONS (a D-0xx findings entry), RUN_STATE. Spec never lags code.

## Open questions (carried to build)

- **Q-N2-1 · redshiftPM periodic background.** The mean-subtracted PM Poisson gives Φ_PM = −GM/r + C(r) with a slowly-varying periodic-image background C. The fit `Φ_meas = −A/r + C_const` absorbs a *constant* C; residual C(r) curvature over r ∈ [8,24] is the tolerance floor. *Deciding measurement:* the fitted A/GM rel error — gate at 5% (poissontest's regime), **measure, don't pre-tune** (D-018). If the periodic curvature dominates, widen the far-field window or shrink the source, declared.
- **Q-N2-2 · fp32 proper-time accumulation.** τ += α·dt over N ticks in fp32 accumulates rounding; the ratio τ(r)/τ(r_ref) is fp32-clean to ~1e-6 (well under the redshift gate). Confirmed cheap at build; if the redshift ratio noise floor rises, read the clock rate from α directly (τ stays the declared field).
- **Q-N2-3 · strong-field probe.** `redshift` probes to r = 2 su (r_s/r = 0.5, z ≈ 0.41). How close to r_s can the shell-mean stay exact before the 1/√ nonlinearity × shell-binning breaks the 1e-3 gate? *Deciding measurement:* the per-radius rel error at build; keep the innermost probe where the gate holds, declared (not a faked-tight inner radius).

## Changelog

- v1.0.0 (2026-07-12) — initial contract + build (operator directive: execute against ROADMAP §3 N2 immediately). The lapse α = √(1+2Φ/c²) and the declared proper-time field τ = ∫α dt; c goes live. Two scenarios: `redshift` (exact Schwarzschild gravitational time dilation — the machinery oracle) + `redshiftPM` (the substrate weld — PM-Poisson well of Newtonian depth GM through the lapse). Reuses N1 `kGreen` verbatim. Honest boundary: temporal metric only (no light-bending/precession — that's N3), static, Schwarzschild-scalar. Tool: `substrate/lapse_nexus.cu`.
