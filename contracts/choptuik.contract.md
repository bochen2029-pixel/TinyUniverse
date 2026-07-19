# choptuik.contract.md — `choptuik_nexus` (direct Choptuik mass-scaling, CPU oracle)

**Version 1.0.0 — FROZEN** (freezes with its goldens; semver from here).
**Approved:** operator directive 2026-07-19 ("proceed to next at your best
recommendations") on the recorded recommendation to golden-gate the scaling
measurement (`tournament/gamma/dss/RESULTS_dss.md`, session-6 close).

## Purpose

Golden-gate the **direct classic measurement of the Choptuik critical exponent**:
supercritical mass scaling M_BH ∝ (p−p\*)^γ on the polar-areal massless-scalar
evolver, so that "γ = 0.37 ± 0.02 (uniform grid)" is a claim a stranger runs cold.
This is the **evolver-side redundant-recovery observable** of the γ crown
(`similarity_nexus.contract.md` v0.1.0 targets the eigenvalue route; its ±0.001
precision is AMR-gated per the measured D-021 triple — this module does NOT claim it).

## Module

`substrate/choptuik_nexus.cpp` — single file, C++17, CPU fp64, stdlib only, no GPU
(runs under any card contention). Faithful port of the research evolver
`tournament/gamma/dss/nr_evolve.py::Ev` (session-6 state, run-4 observables):

- staggered grid r_j=(j+½)dr, N cells, r_max=60; polar-areal metric by outward
  trapezoid log-integration, a(0)=1, outer normalization α·a|_{rmax}=1;
- Φ_t = ∂_r(wΠ) (plain gradient, even ghost), Π_t = ∂_r(wΦ)+2wΦ/r (odd ghost,
  outgoing-ish outer row), w = α/a; KO dissipation ε=0.5; RK4; CFL dt=0.4·dr/max(w);
- Gaussian data φ = p·exp(−((r−12)/2)²) at rest; φ(0,t) EVOLVED (∂_tφ|₀=(α/a)Π|₀);
- freeze detect: max 2m/r > 0.90 ∨ α(0) < 0.02; dispersal: t>40 ∧ inner-third quiet;
- mass observables per run: **M_frz** (freeze-peak) and first-crossing masses
  **M70/M65** at 2m/r = 0.70/0.65 (below the measured freeze floor 0.74, uniform
  across p) with **sub-cell parabolic peak interpolation** (kills the dr staircase).

Critical amplitudes (from the research bisections, reproducible via `--bisect`):
p\*(N=800) = 0.03751655962597 · p\*(N=1600) = 0.03732817692976 ·
p\*(N=3200) = 0.03739102496155509.

## Faces

- `--selftest` — closed-form checks (< 30 s): flat-space metric a≡α≡1 at machine;
  vacuum run disperses; parabolic-peak interpolator exact on a synthetic parabola.
- `--scaling [--golden|--json]` — **THE GOLDEN**: frozen config = the campaign run-4
  config EXACTLY (N=1600, p\*(1600), 26 points log-spaced δp ∈ [1e-4, 1e-2], cut
  rH ≥ 8·dr ⇒ 24 kept); declared = the (δp, M70, M65, M_frz) table + plain-slope
  γ[M70] (%.9e formatting). The 2-decade window is REQUIRED: it holds ≈1.3
  fine-structure periods and partially averages the wiggle — any sub-period window
  gives a wiggle-phase-biased slope (measured: [1e-3, 1e-2] alone → 0.55).
  **Frozen: γ[M70] = 0.3406859239, hash `86c68cf9…` (~13 s).**
- `--cross [--golden|--json]` — **THE ORACLE GATE**: N=1600 at 4 anchor δp; declared
  = masses + max relative deviation vs the committed Python research values
  (`tournament/gamma/dss/gamma_scaling.npy`, run 4, in git):

  | idx | δp | M70 | M65 | M_frz |
  |---|---|---|---|---|
  | 13 | 0.0010964781961431851 | 0.16446089316636592 | 0.15673331402621063 | 0.1758678147376308 |
  | 17 | 0.0022908676527677745 | 0.20794151004571293 | 0.19959739389584485 | 0.23191592433047195 |
  | 21 | 0.0047863009232263802 | 0.31810886500533275 | 0.34128180133591657 | 0.33960666294128899 |
  | 25 | 0.01                  | 0.55267430178984989 | 0.55294965651816241 | 0.54268339619480543 |

- `--ladder N lo hi npts pstar` — general research face (ungated).
- `--bisect N [iters]` — reproduce p\* by bisection (ungated; ~minutes–hours by N).

## Gates (all must PASS for exit 0)

- **G-DET** — declared-state hash byte-identical across two cold runs (two-pass at
  freeze; regression thereafter via `goldens/choptuik_scaling/golden.hash` +
  `goldens/choptuik_cross/golden.hash`).
- **G-GAMMA** — golden-config plain-slope γ[M70] ∈ [0.30, 0.45]. The band is wide
  BY DESIGN: even the 2-decade window (≈1.3 wiggle periods) carries residual
  wiggle-phase bias ~±0.05 (measured, session 6). The exact value is frozen in the
  hash; the band only guards gross regression.
- **G-CUT** — exactly 24 of 26 rows survive the rH ≥ 8·dr resolvability cut (the
  two smallest-δp rows sit below it — measured).
- **G-ORACLE** (cross face) — max relative deviation of the 12 cross masses vs the
  Python table < **1e-6**. **Measured at port time: 0.000e+00 — BIT-EXACT** across
  all 4 anchors × 3 observables (~3000 RK4 steps × 1600 nodes each).
  **Frozen: hash `0e04f941…` (~2 s).**
- **G-FATE** — all golden/cross ladder runs collapse (`fate=bh`); the selftest
  vacuum run disperses.

## Oracle

(a) The committed Python research pipeline (`nr_evolve.py` + `gamma_scaling.npy`,
cross-validated against the literature in RESULTS_dss.md) — implementation oracle;
(b) the literature exponent γ = 0.374 ± 0.001 (Choptuik 1993; Gundlach reviews) —
physics band, entered as the wide G-GAMMA gate honestly (uniform-grid wiggle-phase
systematics ±0.05 documented).

## Honest boundary (declared limits — the module names what it cannot do)

- Uniform grid: the fine-structure period Δ is NOT measurable (window < 1 period in
  ln δp — the D-021 measured triple); the module DETECTS the wiggle indirectly only
  through the frozen residuals.
- γ to ±0.001 (the crown target) requires AMR (N4). This module's claim class is
  **γ = 0.37 ± 0.02** (N=3200 research campaign, RESULTS_dss.md) with the golden
  locking the N=1600 pipeline, not the headline precision.
- Goldens are hardware-pinned (libm; BUILD.md rule). Re-baseline = operator-signed
  note in goldens/.

## Determinism

Single-threaded; no atomics; no wall-clock in the declared path; fixed config in the
golden faces; adaptive dt is a pure function of state. Declared strings via fixed
%.9e formatting. Visuals: none (headless tool).
