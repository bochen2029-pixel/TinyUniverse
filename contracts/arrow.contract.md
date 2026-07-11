# CONTRACT — arrow (M3 thermodynamics + inscription) · v1.0.0 · status: FROZEN 2026-07-11 (operator directive "go M3")

**Purpose.** The arrow of time as declared, gated behavior: coarse-grained entropy meters over the reversible microdynamics, a Loschmidt echo, and the **Ratchet inscription engine** (THE UNFINISHED MIRROR §6; closed form verified in-repo by nexus N6). Extends `app/tinyuniverse.cu` (newton contract still holds).

## Declared machinery

- **Entropy meter:** S = S_x + S_v. S_x from a 32³ position histogram (box/16-su cells), S_v from a 32³ velocity histogram (±8 su/s, clamped); counts via integer atomics (order-invariant); S = −Σ (n/N)ln(n/N) computed host-side fp64 from downloaded counts. Coarse-graining scale is declared, not tuned. HUD and declared output show S and its *change* (the ledger rule: gates on direction/Δ, never on absolute calibration).
- **Ratchet engine:** per-record redundancy n does a biased walk, one event per tick: down w.p. p/(p+(1−p)ρ), else up; absorb at 0 (unwritten) or R+30 (survived); p = 0.3, ρ = 0.6 (λ = p/((1−p)ρ) = 0.714 — supercritical persistence), the N6/ORRERY constants. RNG counter-keyed f(seed, record, tick) — replay-exact. A record is **RECORDED** (regime bit 0x40, reserved in frame.contract v1.0.0, activated here — no bump) once n ≥ 16 (unwrite probability λ¹⁶ ≈ 0.46%: past the point of no refund).
- **Detector:** a slab region that *writes*: particles transiting it gain redundancy deterministically (+1/tick in-slab, cap 64); outside, the record evolves by the ratchet walk. Recorded particles render tinted (visuals non-declared; the bit is declared).
- **SOLV_NONE:** scenarios may declare no gravity (detector/ratchet) — pure drift.

## Scenarios (golden params frozen; seed 20260711)

| scenario | N | solver | init | steps | gates |
|---|---|---|---|---|---|
| **merger** | 10⁶ | pm | two clouds R=60 at (∓80,0,∓10), approach ±1.5 su/s, sub-virial σ, species-colored (4500 K vs 9500 K) | 12000 | S_end − S₀ > 0.3 nat · monotone fraction (ΔS ≥ −0.01 nat per 250-tick sample; fluctuation-scale tolerance, D-015 — collisionless infall transiently compresses) ≥ 0.75 · de_rel < 0.08 |
| **echo** | 10⁶ | pm | merger IC; at tick 6000 momenta flip p → −p (Kahan buffers zeroed, declared) | 12000 | rise = S_mid − S₀ > 0.25 nat · **S_end − S₀ < 0.35·rise** (Loschmidt: the arrow is statistical, not microscopic) |
| **ratchet** | 10³ inert + 3×10⁶ records | none | records: 10⁶ each at R∈{1,4,16} | 4000 | unwrite fraction vs min(1,λ)^R: rel err < 1% (R1) · < 1% (R4) · < 5% (R16, MC-limited) — the in-sim engine reproduces nexus N6 |
| **detector** | 2×10⁵ | none | cold stream v = +2 su/s through slab x∈[0,8] | 4000 | recorded/crossed > 0.95 · crossed > 5×10⁴ |

Declared result adds (per scenario as applicable): `s0, s_mid, s_end, mono_frac, unwrite_r1/4/16, n_crossed, n_recorded`. Record aggregates are declared numbers; the record buffer itself is sim-internal (event-log semantics per frame.contract).

## Gate (M3, from TASKLIST)

dS/dτ ≥ 0 emerges from reversible microdynamics (merger green) · the echo *reverses* the rise (arrow = statistics, proven by its own exception) · in-sim Ratchet matches the closed form · detector writes records before anything collapses (M6 prep). Plus: `harness/verify.py` automates nexus + all scenario goldens.

## Changelog
- v1.0.0 (2026-07-11) — initial freeze under operator directive. P³M/spatial-hash deferred again (collisionless violent relaxation needs no collisions — declared); they arrive with close-encounter physics (M4/M5).
