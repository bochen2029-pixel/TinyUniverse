# MODULE — precession (v1 polish: Q-006 resolution)

**Purpose.** Resolve **Q-006** (the D-016 withdrawal): was the combined SR-inertia + 1PN-field precession really 6.41π, contradicting the 7π superposition? The three fp64 isolation experiments settle it: **the superposition (π + 6π = 7π) is correct**; the app's 6.41π was a **normalization artifact** (nominal vs the orbit's actual, force-distorted semi-latus rectum).

**Contract.** [`contracts/precession.contract.md`](../contracts/precession.contract.md) v1.0.0.

**Tool.** `nexus/precession_nexus.cpp` — single-file **CPU fp64** (no GPU). RK4 orbit integration with pluggable force (Newtonian / harmonic-1PN) + inertia (Newtonian / SR); precession by perihelion tracking + least-squares slope (D-016's method).

**Invariants touched.** 1 · 3 (oracle: exact Sommerfeld + analytic GR) · 8 · 10.

**Oracle.** Exact Sommerfeld relativistic-Kepler precession (sr); analytic GR 6πGM/c²p (pn1); the superposition 7π (combined).

**Build.** `cl /std:c++17 /EHsc /O2 /W4 nexus\precession_nexus.cpp /Fe:build\precession_nexus.exe` (vcvars64; BUILD.md CPU path).

**Run.**
```
build\precession_nexus.exe --scenario sr|pn1|combined --json
build\precession_nexus.exe --scenario NAME --golden
build\precession_nexus.exe --selftest        # Newtonian orbit closes (no precession)
```

**Scenarios / goldens (seed 20260711; measured with the orbit's ACTUAL semi-latus rectum).**

| scenario | force + inertia | measured (actual p) | golden |
|---|---|---|---|
| **sr** | Newtonian + SR | **1.00π** (exact Sommerfeld, 9e-9) — relativistic Kepler, what the app ships | `a0e180df` |
| **pn1** | 1PN + Newtonian | **6.03π** = the correct full-GR precession (6πGM/c²p) | `db0818f2` |
| **combined** | 1PN + SR | **6.95π ≈ 7π** — the superposition CONFIRMED; nominal-p gives 6.45π = the app's 6.41π artifact | `f9df648f` |

`--selftest`: Newtonian force + Newtonian inertia → a closed Kepler ellipse, no precession (1e-10). PASS.

**The resolution (D-026).** π + 6π = 7π **holds** (combined = 6.95π with the actual p, to ~1%). The app's 6.41π came from dividing by the nominal (undistorted) p=6.40 instead of the orbit's actual p=6.90 (the forces expand the orbit ~8%): 6.95 × 6.40/6.90 = 6.45 ≈ 6.41π. So the superposition *reasoning* was right; the *measurement normalization* was the flaw. **The resolution flipped once p_actual was measured** — the initial contract draft had it backwards (RAYFORMER earned its keep). Separately, the combined model's 7π still over-counts the correct full-GR 6π (pn1 = N3 `curve`), so retiring it (D-016) stands.

**Determinism.** fp64 RK4, fixed steps → byte-identical declared JSON → blake2b golden. Both... three goldens two-passed. Runs with the CPU oracles (GPU-independent).
