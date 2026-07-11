# CINEMATIC.md — the rendering law

**No user-facing visual ships unless every box in §7 is checked.** This file exists because agent-written renderers default to "clear, draw, 8-bit swapchain" and look like programmer art. The gap between that and three.js/GARGANTUA is this checklist — ~2–3k lines of pipeline, not a platform. Reference implementation of most of it: `C:\blackhole\files\blackhole.cu` (GARGANTUA).

## 1 · The pipeline (exact order; deviations are DECISIONS.md entries)

```
sim emission (physically-driven HDR color)
  → accumulation buffer  RGBA16F, linear light, no clamping
  → [TAA now | DLSS when integrated]
  → bloom                mip-chain, threshold-free (§3)
  → auto-exposure        log-luminance histogram (§4)
  → tone map             AgX (default) | ACES-Narkowicz (GARGANTUA-parity mode)
  → grade                optional LUT · vignette ≤ 0.15 · grain ≤ 0.5%
  → dither               triangular, BEFORE quantization (mandatory — kills banding)
  → sRGB encode          the ONLY place gamma exists
```

Linear light everywhere internally. Never write raw linear to the swapchain. Never tone map before bloom.

## 2 · Color is physics

- **Blackbody:** Tanner-Helland fit (lift `blackbody()` from GARGANTUA verbatim, incl. the 2.2 linearization).
- **Redshift/blueshift = the g factor**, never a hue rotation: shift T_obs = g·T_em then look up blackbody; intensity scales g⁴ in physical mode (GARGANTUA `shade_disk` pattern).
- Emissive particles sit at 10–10⁴× scene median luminance — real dynamic range is what makes bloom and exposure *do* something.
- Star palette: T ∈ [3200 K, 11000 K] weighted toward cool (GARGANTUA `sky()` pattern).

## 3 · Bloom (mip-chain, CoD-style — upgrade over GARGANTUA's 13-tap)

Downsample 5–6 mips with 13-tap filter, upsample with 3×3 tent + additive blend; final mix 0.04–0.08 scene-referred. **No brightness threshold** (energy-conserving; the threshold look is 2010). GARGANTUA's half-res Gaussian (σ=3.4, bright-pass −0.70) is the compat fallback for the M1 first light only.

## 4 · Exposure

Physical camera: exposure in EV, manual override always available (+/− keys, GARGANTUA parity ×1.15 steps). Auto mode: log₂-luminance histogram, meter 2%–98% clip, temporal adaptation τ ≈ 0.5–1.5 s. For extreme-range views (cosmology timelapse), the Ward-Larson histogram-EQ mode (Buddhabrot v4 `k_eq_accum` pattern) is the sanctioned alternative — declared on the HUD when active.

## 5 · Camera

Damped inertial orbit: critically-damped spring (stiffness 8–12 s⁻¹, damping ratio 0.9–1.0) on azimuth/elevation/log-distance. Wheel zoom exponential ×0.90 per notch, cursor-anchored where meaningful (Buddhabrot). Inclination clamp 1.5°–178.5° (GARGANTUA). FOV presets cycle 45/58/70/85°. View changes crossfade from a reprojection of the previous frame — never flash black (Buddhabrot).

## 6 · The toggle (mandatory)

Every renderer ships **cinematic** and **physical** modes on one key (`C`, GARGANTUA convention). Physical = honest (g⁴ beaming, true dynamic range, no artistic remaps). Cinematic = beautiful, and its deviations from physical are *enumerable on the HUD* (e.g. "beaming softened g⁴→g^1.4; redshift blended 35%"). Educational honesty is a feature, not a mode.

## 7 · The gate (all boxes or it doesn't ship)

- [ ] RGBA16F linear HDR end-to-end; single sRGB encode at the end
- [ ] Tone map: AgX or ACES, after bloom, before grade
- [ ] Mip-chain threshold-free bloom (or declared compat fallback)
- [ ] Exposure: EV-based, manual + auto, temporally smoothed
- [ ] Blackbody + g-factor color (no arbitrary palettes on physical quantities)
- [ ] TAA or DLSS (or SSAA at dev quality); no shimmering point sprites
- [ ] Triangular dither before quantization
- [ ] Damped inertial camera; no raw-input snapping; crossfade on view change
- [ ] Cinematic/physical toggle with HUD-enumerated deviations
- [ ] HUD meters typed: physics meters (τ, dS/dτ, conservation drift) visually distinct from render meters (FPS, exposure)

## 8 · References

GARGANTUA post chain: `blackhole.cu` `kBright`/`kBlur`/`kComposite`/`aces()` (lines ~360–424) · three.js `ACESFilmicToneMapping` + `UnrealBloomPass` (MIT — portable line-for-line) · AgX: Blender 4.x default view transform · Buddhabrot v4 display pipeline (`k_eq_accum → k_tonemap → denoise → composite`).
