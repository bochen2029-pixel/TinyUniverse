# CONTRACT — gargantua (M5 compact objects) · v1.0.0 · status: FROZEN 2026-07-11 (operator directive "go M5")

**Purpose.** Black holes as first-class entities: unscripted formation from collapse, Paczyński–Wiita near-field dynamics (the strong-field oracle einstein v1.0.1 deferred to), Hawking evaporation with an exact energy ledger, and gravitational lensing in the render.

## Declared machinery

- **BH entities** (≤ 8, device array, 1-thread step kernel — fully deterministic): mass M, position, momentum. BHs feel the PM field (trilinear gather at their position) + each other (P–W pairwise); BHs are **not** deposited into the PM grid — their field acts on particles directly (no double counting).
- **P–W force on particles/photons:** g = −GM/(r − r_s)² r̂ per BH, r_s = 2GM/c²; photons take the factor-2 kick (einstein contract); Φ_PW enters τ via the same channel. COMPACT bit (0x08) set within 10 r_s.
- **Absorption:** r < 1.2·r_s → particle DEAD (new regime bit **0x100**, additive to frame contract; parked at +10⁶ su, excluded from deposit/splat/meters — its energy/momentum move to the BH ledger). Deltas accumulate in fixed-point u64 slots during tick T, apply to the BH at T+1 (order-invariant).
- **Formation (unscripted):** every 24 ticks, argmax over PM mass cells (packed-u64 atomicMax — deterministic tiebreak); if M_cell ≥ **M_FORM = 10⁵** (r_s = CELL/4 = 1 su — the dials-audit N_BH crossover, nexus N1) and no BH within 8 su: spawn seed BH at cell center, then absorption at r < CELL for one formation tick binds the core.
- **Hawking:** dM/dt = −K_H/M², K_H = ħc⁴/(15360π G²) ≈ 4.145×10⁵ (dial-derived; lifetime M³/3K_H — M=250 pops in ~12.5 s, M=10⁵ lives 8×10⁸ s). Integrated by the exact cube step M³ ← M³ − 3K_H·dt (deterministic). At M ≤ M_POP = 50 the BH pops; the remainder joins the radiated ledger, so **radiated_total = (M₀ + absorbed − 0)·c² exactly**. No recoil, no emitted quanta in v1 (declared; the render shows T_H-driven glow instead).
- **Lensing (visuals, non-declared):** screen-space point-lens warp per BH before bloom — β = θ(1 − θ_E²/θ²), θ_E = √(4GM/(c²D)), photon-shadow disc at 2.6 r_s/D. Declared approximation; the full GARGANTUA per-pixel Kerr geodesic mode + the OptiX-vs-compute spike (D-007 → ORRERY `lens`) are **deferred to the render/art milestone** — logged, not forgotten (D-017).
- **State hash:** scenarios with BHs append the BH array + ledgers to the declared state; BH-free scenarios' hash domain is untouched (the 12 existing goldens stand).

## Scenarios (seed 20260711)

| scenario | init | steps | gates |
|---|---|---|---|
| **collapse** | cold dense cloud (10⁶, R=40, σ=0.15 v_vir), no seeded BH | 12000 | n_bh ≥ 1 formed unscripted · M_bh grows ≥ 1.5×10⁵ · ledger: absorbed mass = Σ dead-particle mass exactly |
| **isco** | BH M=10⁵ (r_s=1, Newtonian-P–W ISCO=3) + 8 test orbits at r = {2.4, 2.55, 2.7} (plunge set) and {4.5, 5.0, 5.5, 6.5, 8.0} (survive set), SR-circular init, −10⁻³u inward kick | 6000 | all 3 inner absorbed · all 5 outer alive · (declared + **measured**: the SR-inertia+P–W marginally-stable radius lies in (2.7, 4.5) — an r=3.3 orbit plunges; nexus N4 pins the Newtonian-dynamics limit at exactly 3.0. D-017.) |
| **hawking** | BH M=250 + 32 dust at r=5 | 4000 | pop tick vs analytic (M₀³−M_POP³)/(3K_H·dt) within 0.5% · |radiated − M₀c²|/M₀c² < 1e-6 (fp32 BH-mass round-trip; measured 4.3e-8 — D-017) · n_bh 1→0 |

**Ladder scope (honest):** the M5 gate's "cloud → star → remnant → BH" compresses to **cloud → BH** — the ignition/fusion rung needs thermal physics not yet in scope (backlog with M6+). Declared here, not fudged.

## Gate (M5)

Collapse forms a BH unscripted · ISCO plunge/survive split matches P–W · a small BH evaporates on the analytic clock with an exact ledger · 15/15 goldens green.

## Changelog
- v1.0.0 (2026-07-11) — initial freeze under operator directive.
