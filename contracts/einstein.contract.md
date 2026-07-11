# CONTRACT — einstein (M4 relativistic layer) · v1.0.0 · status: FROZEN 2026-07-11 (operator directive "go M4")

**Purpose.** Relativity always-on for every particle (ARCHITECTURE §6: no mode switches — the M2 Newtonian integrator is retired as the approximation it was). **Supersedes the dynamics clauses of newton.contract v1.0.0; all prior scenario goldens are re-frozen under this contract (operator-signed note in goldens/).**

## The declared model (weak-field hybrid; exact GR arrives with M5's compact objects)

- **Momentum:** state stores p = γmv (frame contract unchanged). Velocity recovered cancellation-free: v = p/√(m² + |p|²/c²) — no 1/√(1−β²) anywhere; nexus N7's rapidity discipline satisfied by construction. Massless: v = c·p̂.
- **Massive kick:** dp = m·[g + (1/c²)((v² + 4Φ)g − 4(g·v)v)]·dt — the harmonic-gauge 1PN field term generalized to (Φ, g) from the PM grid / pairwise sums. kGather/kDirect now deliver Φ in acc.w.
- **Photon kick:** dp = (2|p|/c)·g·dt — the weak-field factor-2 bending, declared approximation.
- **Proper time:** dτ = dt·√(max(1 + (2Φ − v²)/c², 0)) — the nexus N8 form with the potential term live (Φ tick-lagged for PM, declared). Photons: the clamp gives dτ = 0 — **light does not age.**
- **Energy meter:** KE = p²c²/(√(m²c⁴ + p²c²) + mc²) (cancellation-free (γ−1)mc²); photons |p|c; PE unchanged (massive only, declared).
- **Precession oracle:** SR inertia alone contributes π·GM/(a(1−e²)c²) per orbit (1/6 of GR); the 1PN field term contributes 6π (nexus N3). The combined model's leading-order expectation is **Δφ = 7π·GM/(a(1−e²)c²)** — declared as *this model's* analytic truth, with the 1/6 overshoot vs pure GR stated openly (retired when M5's pseudo-potential/Kerr oracles take over near compact objects). Cross-terms O(ε²), ε = 5×10⁻³ here.
- Tiny solver: same model in fp64 (u = p/m internal). 2.5PN drag and t_emit rendering: deferred to M5/M7, declared.

## New scenarios (seed 20260711)

| scenario | N | solver | init | steps | gates |
|---|---|---|---|---|---|
| **keprel** | 2 | tiny | m={10⁴,1}, a=10, e=0.6 (ε = GM/ac² = 5e-3), apoapsis, CM frame | 320000 (≈30 orbits) | measured precession (LRL angle, batch-sampled + unwrapped) vs 7πGM/(a(1−e²)c²) rel err < 5% · bounded r < 30 |
| **clocks** | 3 | tiny | central 10⁴; clock A circular r=2 (v≈0.16c, Φ/c²=−0.025), clock B r=100 | 24000 (100 s) | |τ_A/τ_B − analytic| < 2e-3 (analytic = √(1+(2Φ_A−v_A²)/c²)/√(…_B), orbit-averaged) |
| **photons** | 65 | direct | central 10⁴ + 64 photons (m=0, |p|=1) from x=−100, impact b = 2+0.25i | 3000 | mean rel err of deflection angle vs 4GM/(bc²) over b ∈ [6.5, 17.75] < 5% · all photons |v| = c (by construction, asserted) |

**Golden supersession:** kepler/threebody/cloud/galaxy/merger/echo/ratchet/detector re-frozen (dynamics changed for every massive particle; gates unchanged and must still pass — the echo must still reverse exactly under relativistic KDK).

## Gate (M4)

keprel precession green vs the declared 7π oracle · clocks desync green · photons bend at the GR factor · all 12 goldens green via harness.

## Changelog
- **v1.0.1 (2026-07-11) — D-016 erratum: the 1PN field term is WITHDRAWN.** Measured combined precession was 6.41 π-units vs the 7 π-units the linear-superposition argument claimed; the derivation could not be settled under session constraints, so per the RAYFORMER protocol the term is retired rather than shipped unproven (Q-006 holds the open derivation). The declared model is now: **SR inertia (v = p/√(m²+p²/c²)) + Newtonian field + factor-2 photon bending + full proper time** — every piece with an exact oracle. keprel gates against the exact Sommerfeld relativistic-Kepler precession (2π(1/Λ−1), Λ = √(1−(μ/cl)²)) — measured agreement **0.50%**. Also D-016: precession measurement uses a least-squares slope over 320 dense LRL samples (endpoint sampling aliases the intra-orbit oscillation — measured 6.7% artifact, eliminated). Strong-field precession (6π-grade) arrives with M5's pseudo-potential/Kerr oracles. Results: keprel 0.50% · clocks τ-ratio to 6.3e-4 · photons GR bending to 0.83%.
- v1.0.0 (2026-07-11) — initial freeze; 7π combined-model claim (retired above).
