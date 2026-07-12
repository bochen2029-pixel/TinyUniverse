# MODULE — inspiral (v1 polish: 2.5PN gravitational-wave inspiral)

**Purpose.** Resolve the long-owed **2.5PN inspiral** chore. A relativistic binary radiates gravitational waves, losing energy and angular momentum, so its orbit **shrinks and circularizes** — the **Peters (1964)** result, the physics behind every LIGO chirp. Validated to machine precision against exact GR.

**Contract.** [`contracts/inspiral.contract.md`](../contracts/inspiral.contract.md) v1.0.0.

**Tool.** `nexus/inspiral_nexus.cpp` — single-file **CPU fp64** (no GPU; RK4 ODE integration, runs under any card contention). `_nexus` family idiom (blake2b golden, envelope face, exit 0/1/2).

**Invariants touched.** 1 (contract-first) · 3 (named oracle: analytic Peters GR) · 8 (headless envelope) · 10 (fp64 = truth).

**Oracle.** Analytic GR — **Peters (1964)**: circular merger time T_c = (5/256)(c⁵/G³)a₀⁴/(m₁m₂M); eccentric circularization a(e) = a₀·g(e)/g(e₀), g(e)=e^(12/19)/(1−e²)[1+(121/304)e²]^(870/2299).

**Build.**
```
cl /std:c++17 /EHsc /O2 /W4 nexus\inspiral_nexus.cpp /Fe:build\inspiral_nexus.exe
```
(vcvars64 first — BUILD.md CPU path.)

**Run (headless).**
```
build\inspiral_nexus.exe --scenario circular  --json     # merger time vs Peters T_c
build\inspiral_nexus.exe --scenario eccentric --json     # a(e) circularization curve
build\inspiral_nexus.exe --scenario NAME --golden        # vs goldens/inspiral_NAME/golden.hash
build\inspiral_nexus.exe --selftest                      # K=0: no radiation, a unchanged
```

**Scenarios / goldens (seed 20260711).**

| scenario | proves | measured | golden |
|---|---|---|---|
| **circular** | a binary inspirals to merger by GW emission at the Peters rate | m₁=m₂=10⁴, a₀=10: T_merge vs Peters T_c, rel **1.3e-13** | `2eba79de` |
| **eccentric** | the orbit **circularizes** along the Peters a(e) curve as it inspirals | e₀=0.7→0.2: a(e) vs closed form, max rel **5.2e-11**; a shrinks 10.0→2.26 | `4578d3ac` |

`--selftest`: K=0 (radiation off) → a unchanged (< 1e-12); g(e) monotone. PASS.

**Determinism.** fp64 RK4, fixed step counts → byte-identical declared JSON → blake2b golden. No RNG/GPU/atomics. Both goldens two-passed (GOLDEN OK on independent re-run). Runs with the CPU oracles in the harness (GPU-independent).

**Honest boundary (printed in the contract).** Orbit-averaged (secular) **leading-order 2.5PN**, weak-field (a/r_s=50), **point masses, no spin**. It does not integrate the instantaneous 2.5PN equations of motion, and the **app integration** (a radiation-reaction drag term in the app's relativistic KDK so the `merger`/binary scenarios chirp) is named future work, not done here — this module is the physics *oracle*.

**Relation to N3 `curve`.** Together they cover the two GR orbital effects: `curve` does the **conservative** precession (6πGM/c²p, D-024); `inspiral` does the **dissipative** radiation reaction (the chirp, D-025).
