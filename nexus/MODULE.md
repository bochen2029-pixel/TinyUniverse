# MODULE — nexus (tiny_nexus v1.0.0)

**Purpose:** M0 composition proof and the standing fp64 oracle. Proves the four regimes (quantum, classical, relativistic, compact) compose without contradiction under the v0 dials, and freezes the dial defaults. Every later CUDA module gates against nexus values at its declared tolerance.

**Contract:** `contracts/nexus.contract.md` v1.0.0 (operator-approved 2026-07-11; errata **D-011**). Schema: `contracts/nexus.schema.json`.
**Oracle:** the battery *is* analytic — closed-form Kepler/1PN/P–W/QM/Ratchet/SR results (Invariant 3's root of trust).
**Status:** DONE v1.0.0 · golden `ad64f810` · full battery 24.1 s (< 30 s NFR) · N11 in-process ×2 byte-identical + 3× out-of-process `GOLDEN OK`.

## Build

```
cmd /c '"C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat" >nul 2>&1 && cl /std:c++17 /EHsc /O2 /W4 /nologo nexus\tiny_nexus.cpp /Fe:build\tiny_nexus.exe /Fo:build\tiny_nexus.obj'
```

Run from the repo root (`--golden` resolves `goldens/nexus/` from CWD). MSVC 2022 is the golden platform.

## Frozen measurements (runs/nexus_v1.0.0_freeze.json, seed 20260711)

| test | key result |
|---|---|
| N2 kepler | dE/E 7.4e-10 · dL/L 5.8e-13 · peri drift 6.7e-10 rad/orbit (100 orbits, 6.2e7 steps) |
| N3 pn1 | precession 0.150714 vs analytic 0.147262 rad/orbit — 2.3% (gate 5%; residual = higher-order PN) |
| N4 pw | ISCO 3.000000008 vs 3.0 (2.7e-9) · capture from 5.5GM/c² at t=2.62 s |
| N5 qm | free-packet σ err 5.8e-15 · SHO E₀ err 1.3e-13 rel |
| N6 ratchet | rel err R1 2.1e-5 · R4 3.2e-4 · R16 4.5e-3 (gate 1%; MC 1σ at R16 ≈ 0.47%) |
| N7 rapidity | naive-path max err 1.0 (γ=10⁸ diverges) vs cosh path 1.6e-15 — ratio 6.4e14 · γ(0.99999999)=7071.07 |
| N8 tau | conventions lock: 0 · 0 · 1.1e-16 |
| N9 t_emit | static rate exactly 0.8 = 1−v/c · orbiting mean rate 0.79762 (0.3%) |
| N10 compose | all drifts ≤ 2.4e-13 · masks binary=0x12 photon=0x24 orbiter=0x1C stable over 10⁴ ticks |

## Internal design notes

- **Test parameters** (implementer-chosen per contract): N2 m₁=10⁴, m₂=1, a=150, e=0.6, KDK at dial dt; N3 test particle, M=10⁴, a=10, e=0.6, harmonic-gauge 1PN accel, RK4, 30 perihelion passages, precession from LRL-vector angle at perihelia; N4 M=10⁵ (r_s=1 su), ISCO via ternary search on L(r)=√(GMr³)/(r−r_s), capture via −10⁻³·v_φ inward kick; N5 grids 1024/64 su (free, k₀=2) and 512/32 su (SHO ω=1, imaginary-time τ=0.002 × 20000); N9 v_obs=4 (β=0.2), D₀=100, r_orb=2; N10 subsystems dynamically independent by declaration — composition = shared clock/typing/bookkeeping (cross-gravity is the CUDA sim's job).
- **N6 RNG:** mt19937_64 per contract; each 64-bit draw consumed as two 32-bit integer-threshold Bernoullis (quantization bias 2⁻³², five orders below the 1% gate) — this keeps the double-battery run inside the 30 s NFR. Walk capped at R+30 (survival truncation bias λ³⁰ ≈ 4e-5 relative).
- **N7 ladder** extends to γ=10⁸ where β=tanh(acosh γ) rounds to exactly 1.0 in fp64 and the naive path diverges — the cancellation demonstration is deterministic, not rounding-luck-dependent (D-011).
- **Canonical serialization:** `%.9e`, −0 normalized, `notes` excluded from the hash domain. Regime mask constants (N10): QUANTUM 0x1 · CLASSICAL 0x2 · REL 0x4 · COMPACT 0x8 · BOUND 0x10 · MASSLESS 0x20; REL threshold v > 0.05c; COMPACT within 10 r_s. The frame contract freezes the canonical table at M1.

## Known limitations (honest scope, v1.0.0)

- **BLAKE2b-256 is self-contained**, not yet byte-verified against liborrery's envelope; the M1 lift swaps it in — if bytes differ, the golden is superseded under an operator-signed `goldens/nexus/NOTE.md` entry (expected, planned).
- **clang++/g++ parity unverified** on this machine (no non-MSVC toolchain installed) — owed; tracked in TASKLIST M1 chores. Gates are tolerance-based and expected to pass; the golden hash remains MSVC-pinned regardless.
- N6's R=16 gate carries ~2σ headroom at the frozen seed (measured 0.45% vs 1% gate) — deterministic once frozen, documented here for anyone retuning seeds.
- Sims prove structure, never metaphysics: nexus verifies the Ratchet *closed form*; what counts as a "record" in-sim is M3's contract.
