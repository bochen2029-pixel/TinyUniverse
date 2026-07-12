# TASKLIST.md — the build plan

Status legend: ☐ planned · ◐ in progress · ☑ done (golden frozen) · ✖ killed (DECISIONS entry).

## M0 · `nexus` — the composition proof ☑ **CLOSED 2026-07-11**
- ☑ Contract drafted + operator-approved (v1.0.0; errata D-011 applied at implementation)
- ☑ `nexus/tiny_nexus.cpp` implemented (C++17, fp64, single file, MSVC-built)
- ☑ Battery N1–N11 green in 24.1 s; golden frozen `ad64f810`; N11 + 3× out-of-process byte-identical
- ☑ Dial defaults v0 frozen; nexus is the standing oracle
- ☐ *(owed, moved to M1 chores)* clang++/g++ parity build — no non-MSVC toolchain on this machine
- **Gate: MET** — all four regimes stepped in one scene, zero contradictions (N10 masks 0x12/0x24/0x1C stable, drifts ≤ 2.4e-13).

## M1 · `canvas` — first light ☑ **CLOSED 2026-07-11** (P1 presentation owed)
- ☑ frame.contract.md v1.0.0 frozen (Q-001 → god-hand; regime hex table; D-012 P0/P1 clause)
- ☑ liborrery lifted verbatim (`core/lib/`, ORRERY commit d56c4c7, per-file sha256 in LIFT.md); BLAKE2b byte-compat verified (`harness/hash_compat.cpp` green) — **nexus golden stands**
- ☑ CINEMATIC stack in CUDA: float4 HDR → 13-tap mip bloom → log-lum auto-exposure → AgX + ACES parity → astro stretch (cinematic-only, HUD-declared) → triangular dither
- ☑ Two-stream frame loop (sim low-prio ping-pong publish, present high-prio); damped spring camera; CUDA bitmap HUD (D-013; ImGui deferred to P1)
- ☑ Presentation P0: thin GL blit (D-012 — Vulkan SDK absent, measured)
- ☐ *(owed → backlog)* P1 Vulkan swapchain + external-memory presentation (gate: SDK install + DLSS milestone)
- ☐ *(carried)* clang++/g++ parity build of nexus
- **Gate: MET** — 1M particles 1080p: **499 fps avg / 226 min** (SSAA 2×: 178/152); CINEMATIC §7 **10/10** (evidence: `app/MODULE.md`, `runs/firstlight_*.png`).

## M2 · `newton` — the classical tier ☑ **CLOSED 2026-07-11**
- ☑ PM gravity: 128³ CIC (fixed-point uint64 deposit — Invariant 4) → cuFFT Poisson → FD force grid → gather; periodic box = the torus, natively
- ☑ KDK leapfrog + Kahan drift; tiny solver (N≤32) fp64-internal (D-014); direct solver ≤4096
- ☑ Conservation meters via `fixed_atomic_add` (device) / fp64 (tiny); HUD drift readouts
- ☑ Scenario contracts + envelope face (liborrery `full_envelope`/`golden_check`): kepler · threebody · cloud · galaxy — **all gates green, goldens frozen + GOLDEN OK on re-run** (kepler de<1e-6 @10⁶ ticks; figure-8 bounded 29.5 periods; cloud 0.46%; galaxy 1.9%, p_drift 2e-6)
- ☐ *(deferred → M3)* P³M short-range correction (PM cell softening is the declared v1 resolution); spatial hash (arrives with short-range forces); harness/verify.py
- **Gate: MET** — 1M gravitating @1080p: **347 fps avg / 159 min**; 10⁶-tick drift evidence in goldens (kepler/threebody run 10⁶ ticks as their golden params).

## M3 · `arrow` — thermodynamics + inscription ☑ **CLOSED 2026-07-11**
- ☑ Entropy meter (32³ position + velocity histograms, integer-atomic counts, host fp64 S; HUD shows S and ΔS)
- ☑ `merger` golden `34a2db77`: S rises 0.90 nat through violent relaxation, mono 0.77 ≥ 0.75 (D-015)
- ☑ `echo` golden `2f02d94f`: **Loschmidt reversal exact to 6 decimals** — S returns 14.084043 → 14.570027 → 14.084043
- ☑ `ratchet` golden `ccf4a3f8`: in-sim engine vs closed form at 1.5e-5/0.23%/0.33% (R=1/4/16) — nexus N6 parity in production
- ☑ `detector` golden `a0c31f74`: 200k/200k transiting particles inscribed + RECORDED (0x40/0x80 latches live; recorded particles render blue-tinted)
- ☑ `harness/verify.py`: **9/9 GREEN in 55 s** (carried chore closed)
- ☐ *(deferred again, declared)* P³M/spatial hash — arrives with close-encounter physics (M4/M5)
- **Gate: MET** — the arrow emerges from reversible microdynamics, reverses on cue, and records ratchet exactly as the theory's closed form demands.

## M4 · `einstein` — the relativistic layer ☑ **CLOSED 2026-07-11**
- ☑ SR inertia always-on for every particle: v = p/√(m²+p²/c²) (cancellation-free; v<c structural — N7 discipline by construction); photons v = c·p̂ with factor-2 bending; full proper time (photons don't age)
- ☑ `keprel` golden `f985e473`: precession vs **exact Sommerfeld** formula — 0.50% (LS-fit over 320 LRL samples, D-016)
- ☑ `clocks` golden `330c86a7`: τ_A/τ_B across a potential well to 6.3e-4 of analytic
- ☑ `photons` golden `c4c565de`: GR factor-2 deflection to 0.83% over 46 impact parameters
- ☑ All 8 prior goldens superseded + re-frozen (D-016, operator-signed) — **12/12 GREEN in 66 s**
- ✖ 1PN field term — withdrawn after measurement (6.41π vs claimed 7π; Q-006); strong-field precession → M5 pseudo-potential oracles
- ☐ *(deferred → M5/M7, declared)* 2.5PN inspiral drag · Kepler-at-t_emit rendering
- **Gate: MET** — relativistic kepler at the exact analytic rate; clocks desync as GR demands; light bends at the GR factor and does not age.

## M5 · `gargantua` — compact objects ☑ **CLOSED 2026-07-12**
- ☑ BH entities (≤8, device-resident, deterministic): P–W force + Φ channel, absorption → DEAD latch (0x100) + fixed-point ledgers, unscripted formation from the PM density argmax (M_FORM = 10⁵ = the N1 N_BH crossover)
- ☑ `collapse` golden `5bcb5f58`: a 10⁶-particle cloud forms a BH **unscripted** (peak cell 182,720), grows past 1.5×10⁵, absorbed-mass ledger exact
- ☑ `isco` golden `5801ed2f`: inner 3 orbits plunge, outer 5 survive — **measured**: the SR+P–W marginally-stable radius lies in (3.3, 4.5); nexus N4 pins the Newtonian limit at 3.0 (D-017)
- ☑ `hawking` golden `6bd3faeb`: M=250 BH evaporates on the analytic clock (pop tick 2991 = expected 2991, exact) with the energy ledger closed to 4.3e-8 (gate 1e-6, fp32 floor — D-017)
- ☑ Screen-space BH lensing (θ_E warp + shadow disc; the hole's Hawking glow lenses into an Einstein ring — `runs/gargantua_isco_lens.png`) + T_H-driven glow (small BHs visibly run hot as they shrink)
- ☑ **15/15 goldens GREEN in 75 s** (12 pre-M5 goldens unaffected by the DEAD/BH plumbing — verified)
- ☐ *(deferred → render/art milestone, D-017)* GARGANTUA per-pixel Kerr geodesic view · OptiX-vs-compute spike (→ ORRERY `lens` report) · 2.5PN inspiral
- ☐ *(scope note, contract)* the "star" rung (ignition) needs thermal physics — ladder shipped as cloud → BH; fusion backlogged with M6+
- **Gate: MET** — formation unscripted, ISCO split measured and declared, evaporation exact, harness green.

## M6 · `planck` — quantum bubbles ☑ **CLOSED 2026-07-12**
- ☑ ψ engine: 2D 256² split-step cuFFT (real + imaginary time), absorbing edges, wall/barrier/well potentials, host-fp64 observables, counter-keyed inverse-CDF collapse sampling
- ☑ `doubleslit` golden `47a67d66`: **fringes from single collapses — contrast 0.83 (ψ) / 0.82 (4096 dots); the which-way detector kills it to 0.052** (16×); fringe frequency at the analytic k* within near-field tolerance (D-018)
- ☑ `tunneling` golden `f1e7a061`: T = 0.038624 vs same-grid fp64 oracle 0.038623 — **1e-6 agreement** (D-018: oracle isolates implementation, not resolution)
- ☑ `shoq` golden `fa2e009e`: E₀ = 0.500000 = ħω exactly; σ_ground = 0.499995
- ☑ Q-004 resolved (contract): ~5 MB per 64³ bubble → hundreds coexist with the full universe on 16 GB
- ☑ **18/18 goldens GREEN** — reproduced bit-exact on a 100%-contended GPU (D-018 finding)
- ☐ *(deferred → M7, declared)* roaming 3D bubbles in the live universe (spawn-on-isolation, live Ratchet-collapse binding, bubbleVis); live batched-FFT budget probe
- **Gate: MET** — interference emerges one particle at a time and dies under observation; every quantum number oracle-checked.

## M7 · `cosmos` — the tiny planet ☑ **CLOSED 2026-07-11**
- ☑ Contract frozen (`contracts/cosmos.contract.md` v1.0.0); frame contract → v1.1.0 (canonical-range clause + live bubble handle)
- ☑ 3-torus live: canonical coordinates `[−256,256)` in `kDriftK` (exact ±512, Kahan-safe), minimum-image in direct/BH kernels, BH positions wrap. PM field periodic since M2 — coordinates now agree with it
- ☑ `circumnav` golden `3f18f02c`: photon laps the universe **twice** and returns to 3.1e-5 su; 0.5c particle laps once (1.5e-5 su); SR clock τ=44.34 to 5.5e-5; photon τ **exactly 0** (D-019 fp32 leak caught + fixed)
- ☑ `expand` golden `ce448f2b`: EdS comoving cosmology on the PM machinery (Green ×1/a, drift ×1/a², Zel'dovich lattice) — **linear growth D ∝ a measured**: slope 0.988, amplitude ratio 0.7% vs a_end 1.927, H₀=0.011173
- ☑ `bubbles` golden `78b753f1`: roaming 3D quantum bubbles (64³ batched split-step) — spawn on PM-isolation (8/8), free-packet σ to **0.57%** of the analytic law, collapse by Ratchet intrusion (4 collapsed / 4 alive), RECORDED particles stay classical
- ☑ Light-history ring (retarded-time sampling, "see your own past"), 27-image torus splat, **stereographic little-planet projection**, bubbleVis ψ-glow — all non-declared, CINEMATIC-gated
- ☑ 5 goldens superseded (D-019: cloud/galaxy/merger/collapse seam-crossing ejecta + photons τ-fix), physics gates re-verified pass; 13 byte-identical incl. the echo; **21/21 GREEN**
- ☑ Perf: 1M windowed galaxy with 27-image splat + history **225 fps avg / 58 min**; CINEMATIC §7 10/10 (little-planet, lensed, bubbles shots in `runs/cosmos_*.png`)
- ☐ *(deferred → backlog)* Φ-coupled bubbles (V=0 free evolution v1); moving bubble window; a(t) little-planet timelapse video; light-history for BH/bubble glows
- **Gate: MET** — the universe is a globe, wrapped and lensed and self-visible; cosmology's linear growth is a frozen number; quantum bubbles roam and collapse; every regime golden a stranger runs cold.

## v2 SUBSTRATE track (D-020; proposal `docs/PROPOSAL_2026-07-12_v2_substrate.md`)

### N0 · `substrate_nexus` — the substrate oracle ☑ **CLOSED 2026-07-11**
- ☑ Contract frozen (`contracts/substrate_nexus.contract.md` v1.0.0); geometric G=c=1 (D-020)
- ☑ `substrate/substrate_nexus.cpp` — CPU fp64 spherical Einstein–Klein–Gordon, polar-areal constrained evolution (RK4, log-metric integration for positivity, KO dissipation, origin parity, outgoing BC). No GPU.
- ☑ Battery S1–S5 + determinism green in ~20 s; golden `13aa73e5` frozen (N=800, r_max=24); GOLDEN OK out-of-process
- ☑ S1 flat-space mass conservation 7.9e-4 · S2 subcritical dispersal · S3 supercritical horizon (2m/r→0.98, lapse→3.9e-5, M_BH=0.43) · **S4 Choptuik Type-II transition: p*=0.404 bracketed, M_BH 0.412→0.118 → 0 at p* (ratio 0.29), clean power law R²=0.998** · S5 massive-KG stable + conserved 3.3e-6
- ☑ Harness: CPU oracles (nexus + substrate_nexus) run regardless of GPU contention; **22/22 GREEN**
- ☐ *(deferred → AMR contract, D-021)* **precise Choptuik γ ≈ 0.374** — uniform grid caps the self-similar curvature (measured); N0 proves the transition, not the exponent. Also: bound oscillaton eigen-profile; clang/g++ parity.
- **Gate: MET (honest scope)** — the spherical EKG evolution is trustworthy (four green oracles) and the Choptuik *phenomenon* (Type-II, arbitrarily small BHs) is demonstrated and frozen, before any GPU substrate code. Precise γ named as future work (D-021, RAYFORMER discipline).

### N1 · `field` — the Schrödinger–Poisson weld ☑ **CLOSED 2026-07-12** (merged to master `8ace261`)
- ☑ Contract `contracts/field.contract.md` v1.0.2; tool `substrate/field_nexus.cu` (M6 split-step ψ welded to M2 PM cuFFT-Poisson `kGreen` verbatim, one loop)
- ☑ 6/6 goldens GREEN: freepacket `03dd3a3b` + sho3d `dfbc6185` (exact vs nexus N5) · **soliton `d163d765`** (the weld — SP self-similar r_c·M scale-covariant to 3e-8, two-passed) · echoF `433ddcc8` (reversal receipt) · cloudF `2308ea49` (collapse@t_ff) · mergerF `a09dda6a` (two-body attraction)
- ☐ *(deferred → Q-N1-1, D-022)* galaxyF (v1 rotation curve — ψ needs quantized vortices + 512-su box; cloudF is the built irrotational cross-check)
- **Gate: MET** — the PM+ψ gravity weld across three dynamical regimes + both Madelung limits (D-022).

### N2 · `lapse` — the clock ☑ **CLOSED 2026-07-12**
- ☑ Contract `contracts/lapse.contract.md` v1.0.0; tool `substrate/lapse_nexus.cu` (lapse α=√(1+2Φ/c²) + declared proper-time field τ=∫α dt; reuses N1 `kGreen` verbatim; c goes LIVE)
- ☑ `lapse_redshift` `e2c75be5`: **exact Schwarzschild gravitational time dilation** — a clock at r≈2r_s ticks 40% slow (z=0.398), α-err **5.7e-6** vs 1/√(1−r_s/r)−1; τ-integrator 6.6e-6
- ☑ `lapse_redshiftPM` `3dddb950`: **the substrate weld** — the PM-Poisson well through the lapse has Newtonian depth A/GM=0.9643 (3.6%, the PM floor) → correct gravitational redshift from the substrate's own gravity
- ☑ `--selftest` (flatlapse): Φ=0 → α≡1 exact, τ=N·dt to 6.6e-7. Harness rows added (behind the GPU preflight); both goldens two-passed cold.
- ☐ *(deferred → N3 `curve`, D-023)* the spatial metric — light bending, orbit precession, Shapiro delay (N2 is temporal-metric-only: exact redshift, no bent paths)
- **Gate: MET** — clocks run slow in gravity wells, on the substrate; the redshift is exact Schwarzschild and the substrate's own gravity produces it (D-023).

### N3 · `curve` — geometry curves ☑ **CLOSED 2026-07-12**
- ☑ Contract `contracts/curve.contract.md` v1.0.0; tool `substrate/curve_nexus.cpp` (CPU fp64 geodesic oracle — no GPU, runs under any contention, like N0)
- ☑ `curve_deflect` `4e6c33ca`: **light bends at exact GR 4GM/bc²** (0.3–0.75%) — the 1919 **factor of 2 decomposed**: the lapse (time, N2) gives 2GM/bc², the spatial curvature (N3) doubles it (full/lapse = 2.004–2.008 vs exactly 2)
- ☑ `curve_precess` `67272705`: **perihelion precession** 0.72°/orbit vs exact GR 6πGM/(c²a(1−e²)) — **0.52%** (weak-field, leading 1PN)
- ☑ `curve_shapiro` `20bfd4d2`: **Shapiro time delay** (excess light-travel time past the mass) vs (4GM/c³)asinh(X₀/b) — **0.33%**. With N2's redshift, this completes **all FOUR classical tests of GR** from the substrate metric (redshift · light-bending · precession · Shapiro)
- ☑ `--selftest` (flat GM=0): straight ray + closed orbit. Harness rows added to the CPU-independent oracles (no GPU preflight); both goldens two-passed.
- ☐ *(deferred → N4+, D-024)* a live GPU metric field a²(x) with dynamical back-reaction (energy density sourcing the metric each step); strong-field / higher-PN
- **Gate: MET** — geometry curves: light bends at the exact GR angle (space curvature is the difference over N2), orbits precess at the GR rate (D-024).

### N4 (GPU, future — gated against N0 + the 21 v1 goldens)
- ☐ N4 `star` (fusion closure + radiation + Ratchet lattice — the hydrogen-ball sentence)

## Backlog (measured adoption only; each needs an honest baseline)
- ☐ DLSS 4.5 via Streamline · ☐ RTXDI/ReSTIR emissive-particle lighting · ☐ NRD · ☐ Slang port of shared math (D-008) · ☐ UE 5.8 TextureShare showcase shell · ☐ MCP surface for scenario driving · ☐ df64 zoom ladder (D-010 gate)
