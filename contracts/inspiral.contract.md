# CONTRACT — inspiral (v1 polish: 2.5PN gravitational-wave inspiral) · v1.0.0 · status: FROZEN 2026-07-12 (CLOSED — 2/2 goldens green: circular `2eba79de` · eccentric `4578d3ac`; D-025)

**Purpose.** Resolve the long-owed **2.5PN inspiral** chore (TASKLIST M4/M5 "deferred → 2.5PN inspiral drag → binaries inspiral"; operator queue #3, v1 polish). A relativistic binary radiates **gravitational waves**, losing energy and angular momentum, so its orbit **shrinks and circularizes**, spiralling to merger. This is the leading-order (2.5PN) radiation-reaction effect — the **Peters (1964)** result — and it is the physics behind every LIGO chirp. A clean CPU fp64 oracle, gated against exact GR (analytic Peters formulae), in the `_nexus` family idiom.

**The physics (Peters 1964, orbit-averaged 2.5PN radiation reaction):**

```
da/dt = −(64/5)(G³/c⁵)·(m₁m₂M/a³)·(1−e²)^(−7/2)·(1 + (73/24)e² + (37/96)e⁴)
de/dt = −(304/15)(G³/c⁵)·(m₁m₂M/a⁴)·e·(1−e²)^(−5/2)·(1 + (121/304)e²)
```
(M = m₁+m₂.) Two exact consequences are the oracles:
- **Circular merger time:** `T_c = (5/256)·(c⁵/G³)·a₀⁴/(m₁m₂M)` — the time for a circular binary (e=0) to spiral from a₀ to a=0.
- **The circularization curve:** integrating da/de gives `a(e) = c₀·e^(12/19)/(1−e²)·[1 + (121/304)e²]^(870/2299)` — the binary sheds eccentricity as it inspirals (Peters' famous a(e) relation; the orbit is nearly circular by merger).

Single-file CPU fp64 oracle `nexus/inspiral_nexus.cpp` — no GPU (RK4 ODE integration; deterministic, runs under any card contention). Envelope face (`--scenario X --json|--golden|--selftest`), exit 0/1/2, blake2b golden (the CPU-oracle idiom — the measured numbers are the declared state).

## Units & dials (v1 dials)

| dial | value | role |
|---|---|---|
| c | 20 su/s | radiation reaction ∝ 1/c⁵ (leading 2.5PN); r_s = 2GM/c² |
| G | 2×10⁻³ | with the masses sets the coupling K ≡ G³m₁m₂M/c⁵ |
| m₁, m₂ | 10⁴ each (M=2×10⁴) | the binary; r_s = 2GM/c² = 0.2 su, so a₀ = 10 su gives a/r_s = 50 — the weak-field PN regime where the leading 2.5PN result is the answer |

## Declared engine

**Integrators (fp64 RK4).**
- **Circular inspiral:** integrate `a(t)` via da/dt (e=0) from a₀ forward in t until a drops below a_end = 1 su (a/r_s = 5); T_merge = that time. (Analytically a⁴ is linear in t → a(t) = (a₀⁴ − (256/5)Kt)^{1/4}; RK4 must recover it.) The GW chirp f_GW ∝ (T_c − t)^{−3/8} follows.
- **Eccentric circularization:** integrate `a(e)` via `da/de = (da/dt)/(de/dt)` (RK4 in e as the independent variable) from e₀ = 0.7 down to e checkpoints; compare to the closed-form Peters a(e). This tests the coupled radiation-reaction dynamics, not just the rate.

**Determinism.** Pure fp64 RK4, fixed step counts → byte-identical declared JSON → blake2b golden. No RNG, no GPU, no atomics. CPU (sm-independent).

## Scenarios & gates (the golden units)

| scenario | what runs | gate | oracle |
|---|---|---|---|
| **circular** | m₁=m₂=10⁴, a₀=10, e=0; integrate a(t) → T_merge (to a_end=1) | **`|T_merge/(T_c·(1−(a_end/a₀)⁴)) − 1| < 1e-3`** | **analytic Peters** T_c = (5/256)c⁵a₀⁴/(G³m₁m₂M). The 2.5PN inspiral rate + the chirp. |
| **eccentric** | m₁=m₂=10⁴, a₀=10, e₀=0.7; integrate a(e) → a at e ∈ {0.6,0.5,0.4,0.3,0.2} | **`max |a_num(e)/a_Peters(e) − 1| < 1e-3`** (the orbit circularizes along the Peters curve) | **analytic Peters** a(e) = a₀·g(e)/g(e₀), g(e)=e^(12/19)/(1−e²)·[1+(121/304)e²]^(870/2299). |

**`--selftest`:** K=0 (radiation off) → da/dt = de/dt = 0 → a, e constant (drift < 1e-12). The no-radiation sanity + determinism smoke.

Both gates are exact against closed-form GR (Peters 1964); the 1e-3 tolerance is RK4/step-count residual (measured, not pre-tuned — D-018). No scenario is gated against an unmeasured number.

## The headless face (envelope-conformant)

```
inspiral_nexus.exe --scenario circular|eccentric [--seed N=20260711] --json    # declared JSON
inspiral_nexus.exe --scenario NAME --golden                                    # vs goldens/inspiral_NAME/golden.hash
inspiral_nexus.exe --selftest                                                  # K=0: a,e constant
```

Declared body (canonical order): `seed, params{scenario, c, G, m1, m2, a0, e0, a_end, steps}, result{<T_merge, T_c, T_rel, chirp_exp | a_e checkpoints, a_rel_max>, nan_free}, gates{primary, nan}, verdict`. fmt6/fmt9. Exit 0/1/2. Build with MSVC `cl` (BUILD.md CPU path).

## Oracle pedigree (Invariant 3)

- **Exact analytic GR:** Peters (1964) — the circular merger time T_c and the eccentric circularization curve a(e). Oracle-gate failure ⇒ run-state SUSPECT.
- **v1 tie:** this is the 2.5PN radiation reaction the app's `einstein`/`gargantua` modules deferred (M4/M5 chores). The oracle validates the physics; app integration (a drag term in the relativistic KDK so `merger`/binary scenarios chirp) is a separate later step, named here.

## The honest boundary (permanent, printed)

1. **Orbit-averaged (secular) 2.5PN — the Peters result.** The oracle integrates the orbit-averaged energy/angular-momentum loss (Peters 1964), the leading radiation reaction. It does **not** integrate the instantaneous 2.5PN equations of motion (the oscillating reactive force) — that is a possible extension; the secular inspiral is the standard, unambiguous result and what governs the observable chirp.
2. **Leading order (2.5PN), weak-field.** Valid for a ≫ r_s (the oracle runs a/r_s ≥ 5); higher PN orders and the plunge near merger are out of scope. The formulae are used down to a_end where the PN expansion is still controlled.
3. **Point masses, no spin.** No spin-orbit/spin-spin, no tidal effects — the leading quadrupole radiation of two point masses.

## Build runbook

1. **Contract (this file).** ← *here.*
2. **`nexus/inspiral_nexus.cpp`** — blake2b + fmt + the Peters ODE integrators + envelope/golden/main. CPU fp64.
3. **Build** — `cl /std:c++17 /EHsc /O2 /W4 nexus\inspiral_nexus.cpp /Fe:build\inspiral_nexus.exe`.
4. **Golden** — `goldens/inspiral_circular/golden.hash`, `goldens/inspiral_eccentric/golden.hash`.
5. **Harness** — an `INSPIRAL_GOLDENS` block runs with the CPU oracles (GPU-independent, no preflight).
6. **Register** — TASKLIST (close the 2.5PN chore), DECISIONS (D-025), RUN_STATE, `nexus/MODULE_inspiral.md`. (Not an ARCHITECTURE §11 module row — it's a v1-polish oracle, noted in the chore closure.)

## Changelog

- v1.0.0 (2026-07-12) — initial contract + build (operator queue #3, v1 polish). The 2.5PN gravitational-wave inspiral (Peters 1964): circular merger time T_c and eccentric circularization a(e), vs exact GR. CPU fp64 oracle `nexus/inspiral_nexus.cpp`. Honest boundary: orbit-averaged (secular) leading-order 2.5PN, weak-field, point masses.
