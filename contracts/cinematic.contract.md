# cinematic.contract.md — `cinematic` (R1: the CINEMATIC post-chain on the interop path)

**Version 1.0.0 — FROZEN** (approved 2026-07-19, operator: "okay try supernova first
… proceed" [scene pick + green-light]; frozen with golden `cinematic_r1` `4962558c…`,
two-passed same day; `--shot` frames delivered for the judged half of G-CHECKLIST).
**Amendment at freeze (measured, two steps):** the §4 exposure meter for astro
content is **energy-weighted** (key = ΣE·log₂L/ΣE over lit bins) — the percentile
2–98% meter was blind to the flash twice (bin-0 sky domination, then
percentile-robustness against small bright regions, EV moved 0.03); a camera is
blinded by a point flash because the flash carries the frame's ENERGY. The
percentile form remains right for non-astro scenes; the meter choice is per-scene
and declared. Basis: `CINEMATIC.md` (the rendering law — binds from R1
by design) · R0 `interop` v1.0.0 (the proven presentation path, golden `4ba7fbcb`) ·
the M1 canvas stack in `app/tinyuniverse.cu` (CINEMATIC 10/10, the LIFT SOURCE) ·
GARGANTUA `C:\blackhole\files\blackhole.cu` (the reference implementation).

## Purpose

R1 moves the **entire CINEMATIC §1 pipeline** onto the R0 presentation path:

```
HDR test scene (blackbody starfield + hot emitters, deterministic)
  → RGBA16F linear accumulation (no clamping)
  → mip-chain threshold-free bloom (§3)
  → log-luminance auto-exposure, EV-based, temporally smoothed in TICKS (§4)
  → AgX tone map (ACES-Narkowicz = GARGANTUA-parity mode)
  → triangular dither → single sRGB encode → 8-bit pack (fixes the R0 R/B swizzle)
  → the R0 shared buffer → swapchain
```

The scene is synthetic BUT physical (blackbody temperatures, 10–10⁴× median-luminance
emitters — §2): rich enough that bloom and exposure demonstrably WORK, simple enough
to be byte-deterministic. R2 replaces the scene with the substrate; the pipeline and
its golden survive unchanged.

## Module

`render/cinematic.cu` — single file, CUDA + Vulkan (the R0 plumbing pattern repeated;
factoring a shared header waits for the rule of three at R2). Fixed 1280×720.

**Lifted verbatim where proven (D-005 spirit; sources named in MODULE.md):**
GARGANTUA `blackbody()` (Tanner–Helland fit incl. 2.2 linearization) · the app's M1
mip-bloom and AgX/ACES kernels · the D-013 bitmap-font HUD kernel · the damped
inertial orbit constants (§5: stiffness 8–12 s⁻¹, ζ 0.9–1.0, wheel ×0.90).

## Faces

- **windowed (default)** — the user-facing view: damped inertial orbit camera (§5),
  `C` toggles cinematic/physical with HUD-enumerated deviations (§6), typed HUD
  meters (§7: physics-cyan vs render-amber), manual EV +/− keys.
- **`--headless --frames N [--golden|--json]`** — fixed camera, deterministic
  tick-driven scene + exposure adaptation; declared = blake2b of the FINAL 8-bit
  frames {0, N/2, N−1} + the pipeline config + the gate booleans.
- **`--shot PATH.png`** — one frame to disk for the operator's visual review (the
  aesthetic half of the gate is operator-judged, as at M1).
- **`--selftest`** — closed-form KATs, < 10 s (see gates).

## Gates

- **G-DET** — two-pass byte-identical declared hash (`goldens/cinematic/golden.hash`).
- **G-KAT** (selftest) — AgX at pinned sample values · blackbody(3200/6500/11000 K)
  vs the Tanner–Helland reference triples · sRGB encode round-trip · triangular-dither
  mean/variance bounds · mip-chain energy conservation to ≤ 1% on a flat field.
- **G-RANGE** — the headless scene's in-frame luminance range ≥ 10³ (bloom and
  exposure are exercised, not decorative), asserted in the declared string.
- **G-EXPOSURE** — the auto-exposure settles (declared EV at frame N−1 within a
  pinned band) and manual override works (selftest step).
- **G-VALIDATION** — Khronos layers ON headless: 0 errors, 0 warnings (R0 idiom).
- **G-CHECKLIST** — the §7 boxes, split honestly:
  *mechanical* (golden-pinned): RGBA16F end-to-end · single sRGB encode · tone map
  after bloom · threshold-free mip bloom · EV exposure · blackbody color · dither
  before quantization; *judged* (operator sign-off on `--shot` output, recorded in
  TASKLIST like M1): overall look, camera feel, HUD typing.
  AA box: satisfied at **SSAA-dev level** (2× supersample of the point splats) —
  TAA/DLSS explicitly deferred to the DLSS milestone (declared, not skipped).

## Honest boundary

Synthetic scene (the substrate arrives at R2 — this contract's scene face is then
swapped, pipeline + goldens retained). Fixed resolution; no HDR swapchain; no TAA;
no LUT grading at v1 (vignette/grain within §1 bounds allowed). The windowed face
inherits R0's FIFO vsync; pacing non-declared (frame contract).

## Determinism

Headless face only: tick-parameterized scene and adaptation (NO wall-clock), fixed
frame count, single queue, house blake2b. Windowed face non-declared.

## Open questions — RESOLVED at approval

- Q-R1-1 → **the SUPERNOVA scene** (operator's pick): ~2500 counter-hashed blackbody
  stars (T ∈ [3200, 11000] K cool-weighted, lognormal-ish L) as the backdrop; star 0
  is the progenitor — light curve quiescent → ×10⁵ flash (hot blue ~15000 K) →
  cooling power-law decay, all pure functions of the tick, looping. The flash IS the
  G-RANGE source and the exposure-adaptation demo. (The other brainstormed presets —
  inspiral, Choptuik-breathing, ladder — join later as windowed-face additions,
  non-declared, no contract churn.)
- Q-R1-2 → minimal-AgX polynomial (portable, KAT-able); ACES-Narkowicz as the
  GARGANTUA-parity mode on the T key.
- Q-R1-3 → the golden pins ONLY the final packed 8-bit frames (+ declared EV values);
  fp intermediates stay unfrozen so look-tuning doesn't churn goldens.
- Amendment (AA box): satisfied at R1 by **analytic band-limited splats** (gaussian
  footprint σ ≥ 1.2 px — the shimmer the checklist bans comes from sub-pixel point
  sprites; band-limited splats cannot alias). SSAA/TAA/DLSS remain the upgrade path,
  declared not skipped.
- Amendment (precision): the linear HDR intermediates are fp32 (`float4`) — strictly
  exceeds the RGBA16F floor of §7 box 1; declared in MODULE.md.
