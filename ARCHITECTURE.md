# TINY UNIVERSE — Architecture

### A scale-compressed, contract-bounded, deterministic universe: quantum → Newtonian → relativistic → black holes, emergent within one particle substrate, on one GPU.

**Version:** 0.1.0 (spec-first; no module source yet) · Founded 2026-07-11 · Owner: Bo Chen · Built with Claude (Fable 5).

---

## 1 · System Intent

In nature the physics regimes are separated by ~60 orders of magnitude because ħ is tiny, c is huge, G is feeble. TINY UNIVERSE retunes the three constants so every crossover lands between 1 and ~10⁶ particles: **one substrate, one rule set, all regimes on one screen.** It is simultaneously a real-time windowed experience (beautiful under `CINEMATIC.md`) and a headless instrument (deterministic, golden-gated, oracle-checked) — the same binary, two faces.

Success = a windowed universe whose every physics claim is backed by a green oracle gate a stranger can reproduce cold, and whose every shipped frame passes the CINEMATIC checklist.

## 2 · The Split (the top-level contract)

```
  THE CLIENT  (app/: window, input, renderer)      THE INSTRUMENT  (core/: libtinyverse)
  ──────────────────────────────────────────       ─────────────────────────────────────
  presents, interacts, beautifies                  simulates, measures, conserves
  reads GPU state ────────────────►  FRAME CONTRACT (SoA buffers + fields + sync semantics)
  writes input events ────────────►  INPUT TRACE (recorded; replay = determinism)
  visuals/timing = NON-DECLARED                    declared state = hashable, golden-gated
```

The renderer depends on the core **only** through the frame contract (`contracts/frame.contract.md`). Break that seam and the renderer becomes physics — and physics becomes un-goldenable. The seam is also what makes the renderer swappable (Vulkan today; UE 5.8 TextureShare showcase shell someday; D-002).

**Sibling, not child, of ORRERY** (`C:\ORRERY`): same doctrine, different clock (interactive 240 Hz ticks vs batch campaigns). TU never imports ORRERY source except the **verbatim liborrery lift** (D-005). ORRERY is referenced, never edited, from here.

## 3 · The Doctrine

> **The spec is the product. The contract is sacred. The golden is load-bearing. The code is ephemeral.**

Adopted whole from ORRERY (`C:\ORRERY\ARCHITECTURE.md` §3), with one adaptation for interactivity: determinism is defined over **(dials, seed, input trace, step count)** — the game-industry replay solution welded to the golden discipline. See I-2.

## 4 · Glossary (canonical; forbidden synonyms in parens)

| term | definition |
|---|---|
| **Dial** | A fundamental constant of this universe (ħ, c, G, m_p, dt, L_box), fixed in `dials.json`, semver'd with the frame contract. (not: "setting", "config" loosely) |
| **Regime** | An emergent physics behavior band (QUANTUM, CLASSICAL, RELATIVISTIC, COMPACT/BH), carried per particle as a **derived-only** bitmask (liborrery `regime.h` pattern) — computed from state, never stored authoritative. |
| **Frame contract** | The sacred seam: SoA particle buffers + field textures + sync semantics the renderer may read. `contracts/frame.contract.md`. |
| **Scenario** | A named, seeded initial condition + dial set, runnable headless (`--scenario doubleslit --steps N --seed S --json`) or windowed. Scenarios are the golden units. |
| **Input trace** | A recorded, timestamped (in ticks) stream of user interventions. Replay of a trace is bit-exact. Empty trace = pure simulation. |
| **Golden** | Frozen `(dials, seed, trace, steps) → blake2b-256 of declared state` in `goldens/`. Hardware-pinned sm_89. |
| **Oracle** | The independent reference a CUDA module is judged against: analytic result or `tiny_nexus` (CPU, fp64). Named per module (I-3). |
| **Quantum bubble** | A local comoving grid holding one particle's ψ, evolved by split-step cuFFT; created on isolation, collapsed by inscription. |
| **Inscription / Ratchet** | The decoherence mechanic (D-004): a quantum state collapses when redundantly recorded into environment DOF; unwrite probability `min(1, p/((1−p)ρ))^R`, crossover at `(1−p)ρ = p`. Parameters anchored to ORRERY `ratchet` receipts. |
| **Proper time τ** | Per-particle clock, `dτ = dt·√(1 − v²/c² − 2Φ/c²)` (weak field). Drives all visible animation rates. |
| **Declared output** | What goldens hash: particle state, conserved quantities, event log. **Non-declared:** frames, wall-clock, FPS, bloom. |
| **CINEMATIC gate** | The `CINEMATIC.md` checklist. A visual that fails it does not ship. |

## 5 · Invariants (system physics — always hold)

1. **Every module has a contract** before source exists; scenarios have schemas; breaking change = MAJOR bump + migration note.
2. **Determinism:** fixed timestep (dt is a dial, decoupled from render rate); `(dials, seed, trace, steps) ⇒` byte-identical declared state. All randomness counter-based: `f(seed, id, tick)` — stateless, replay-safe (liborrery `rng.cuh`). Wall-clock, frame pacing, and every pixel are non-declared.
3. **Every physics module names its oracle** (analytic or `tiny_nexus` fp64) in its MODULE.md and §11. Oracle-gate failure sets run-state **SUSPECT** — never a silent fallback. (RAYFORMER ADR-007 lesson.)
4. **Conservation gates are declared:** energy, momentum, angular momentum drift bounds per scenario, computed with deterministic reductions only (liborrery `reduce.cuh` shapes; no float atomics in declared sums).
5. **Regimes emerge — no hard switches.** The only representation handoffs (ψ ↔ point particle) are physics-triggered (Ratchet redundancy threshold crossed), declared in the event log, and deterministic.
6. **Unbiased acceleration** (Buddhabrot v4 rule): every performance path (mixed precision, importance sampling, graphs, batching) ships with a proof or paired-oracle test that declared expectations are unchanged within declared tolerance — or it does not ship. `--use_fast_math` banned in declared paths.
7. **liborrery verbatim rule:** lifted files (`rng.cuh`, `reduce.cuh`, envelope, `ckpt.h`, `regime.h`) are byte-identical to `C:\ORRERY\lib` at the pinned commit; any divergence is a named fork in DECISIONS.md.
8. **The headless face is a full citizen:** every scenario runs windowless with the universal envelope (JSON + exit 0/1/2, `--selftest`, `--golden`). Exit 1 = declared gate fired (real negative result); exit 2 = error. Never conflated.
9. **CINEMATIC.md gates every shipped visual**; the cinematic/physical toggle is mandatory, and physical mode's deviations from cinematic mode are enumerable on the HUD.
10. **fp32/fp64 policy:** the CUDA core runs fp32 (df64 ladder if deep zoom demands it, per Buddhabrot); `tiny_nexus` is fp64 and defines truth; per-module tolerance vs oracle is declared in the contract. Host≠device libm pinned separately (liborrery KAT lesson — never assert host==device transcendentals).
11. **Frozen external data:** any external dataset ships with sha256 + provenance; goldens cover the hash. (None expected in v1.)

## 6 · The Physics (what the dials buy)

**Dials — v0 strawman (frozen for M0 review in `contracts/nexus.contract.md`; tuned only via nexus N1):**

| dial | v0 value | consequence |
|---|---|---|
| m_p (particle mass) | 1 sm | one massive species in v1, + photons (Q-002) |
| c | 20 su/s | orbital speeds hit 0.05–0.3 c; r_s = 2GM/c² is visible; light lag perceptible across the box |
| ħ | 0.5 sm·su²/s | λ_dB(v=1) = 0.5 su — a slow lone particle is visibly wavy; Hawking T ∝ ħ makes small BHs pop |
| G | 2×10⁻³ su³/(sm·s²) | r_s(10⁵ particles) ≈ 1 su; clouds of 10⁴ self-bind at v/c ≈ 0.07 |
| dt (tick) | 1/240 s | 4 ticks per 60 Hz frame; fixed, never adaptive globally (hierarchical sub-stepping is per-particle and deterministic) |
| L_box | 512 su | 3-torus period (M7) |

**The regime ladder (mechanism → algorithm → oracle):**

| tier | mechanism | algorithm | oracle |
|---|---|---|---|
| QUANTUM (M6) | lone/isolated particles carry ψ | quantum bubbles: local comoving 64³–128³ grid, split-step cuFFT; spawn on isolation, collapse on inscription | free-packet σ(t), SHO E₀ = ħω/2 (nexus N5) |
| CLASSICAL (M2) | decohered bulk | SoA leapfrog (KDK, symplectic); spatial hash short-range; PM gravity (cuFFT Poisson) + P³M near-field | two-body Kepler (nexus N2); conservation gates |
| ARROW (M3) | thermodynamics + decoherence | temperature/pressure/entropy meters (deterministic reductions); **Ratchet inscription**: records written into environment DOF, redundancy R per bubble, collapse past `(1−p)ρ = p` | ORRERY `ratchet` closed form (nexus N6); entropy meter shows **dS/dτ only** — flow, not level |
| RELATIVISTIC (M4) | always-on corrections, reachable because c is small | p = γmv integration (rapidity form, γ = cosh ω — no cancellation); per-particle τ; 1PN precession from Φ field; photons as c-speed particles; 2.5PN drag → binaries inspiral | precession Δφ = 6πGM/(a(1−e²)c²) (N3); τ vs exact SR/weak-field (N8); t_emit retro-phase (N9) |
| COMPACT/BH (M5) | gravity wins: R < r_s(M) | BH entity absorbing boundary; **Paczyński–Wiita** Φ = −GM/(r − r_s) for disk dynamics (ISCO at 6GM/c²); Hawking T ∝ ħ/M — small BHs visibly evaporate; interior/near view = GARGANTUA Kerr renderer with dynamic mass | P–W ISCO (N4); analytic Kerr (GARGANTUA's own validated integrator) |

**The mass ladder (emergent progression):** gas cloud → collapse → ignition (density/temperature threshold) → star → exhaustion → degeneracy-pressure remnant (ħ is big: white-dwarf physics at toy scale) → over the limit → black hole. "Keep adding particles and matter itself changes phase."

**Topology & cosmology (M7):** space is a 3-torus (L_box); light wraps — with small c and a light-history ring buffer you see your own system's past. Zoomed out, render via stereographic **little-planet projection** (the product's visual signature). Cosmology mode: scale factor a(t), hot-dense start → structure → heat death as timelapse.

**Framing note (positioning, not physics):** the Ratchet mechanic and per-seat clocks are the measured mathematics from ORRERY's receipts; TINY UNIVERSE renders mechanisms and proves structure. It is not a proof of any theory, and takes no metaphysical position (Q-003 holds the branding question).

## 7 · The Frame Contract (v0 sketch — DRAFT; frozen at M1)

Vulkan-allocated, CUDA-imported external memory (`cudaImportExternalMemory` + external semaphores; zero copy). Sim ticks (fixed dt, low-priority stream) and render frames (high-priority stream) are decoupled — **presentation never waits on accumulation** (Buddhabrot two-stream pattern; no device-wide syncs in the frame loop).

| buffer | type | notes |
|---|---|---|
| pos | float3 SoA | su; df64 ladder only if zoom demands (D-010 gate) |
| mom | float3 SoA | relativistic momentum p = γmv |
| mass / temp | float | temp is derived for HUD, non-authoritative |
| tau | float | per-particle proper time — drives all visible animation |
| regime | u32 | derived-only bitmask (`regime.h` pattern) |
| bubble | u32 | ψ-handle (0 = none); bubbles own local grids |
| Φ field | 3D texture | PM potential — reused for τ, 1PN, lensing |
| event log | append buffer | inscriptions, collapses, ignitions, captures — declared output |

Full schema, ownership, and semaphore semantics land in `contracts/frame.contract.md` v1.0.0 at M1.

## 8 · Rendering ladder (each rung measured before adoption — ADR-007 protocol)

1. **M1:** Vulkan interop canvas + the GARGANTUA post chain ported per `CINEMATIC.md` (HDR → mip bloom → AgX/ACES) — gate: 1M inert particles, gorgeous, 60 Hz.
2. **M5:** geodesic/lensed rendering — OptiX via function-table dispatch against the driver DLL, compile-out fallback (Buddhabrot recipe). This **is** ORRERY's parked `lens`/RT spike: RT cores vs CUDA-compute baseline, honestly measured; result flows back to ORRERY either way.
3. **Backlog, measured:** DLSS 4.5 (Streamline), RTXDI/ReSTIR (emissive particles as real lights — the accretion-disk look), NRD, Slang single-source math (D-008), UE 5.8 TextureShare showcase shell (display-only; never the home).

## 9 · Verification model (compile-as-proof)

`harness/` CI: compile every module → run every `--selftest` → run every golden → red/green. `tiny_nexus` (M0) is the fp64 oracle spine: the CUDA core is later gated against its values at declared tolerances. Two-pass (fresh cold-context verifier) mandatory for any publicly citable claim.

## 10 · Tech decisions

ADRs in `DECISIONS.md`. Summary: CUDA 13.1 · sm_89 (fat binary sm_89+sm_90 SASS, compute_120 PTX) · MSVC 2022 · fast-math banned · GLFW/GLAD/ImGui/stb via FetchContent · Vulkan default (Q-005) · OptiX stubs recipe · one binary, two faces · liborrery verbatim lift · Ratchet decoherence · Slang at M4+ · single-file-per-module where possible, CMake preset for multi-file.

## 11 · Module inventory (build order in `TASKLIST.md`)

| module | lang | what it does | seeded from | oracle | status |
|---|---|---|---|---|---|
| **nexus** | C++17 fp64 | the composition proof: all four regimes step together without contradiction; freezes dial defaults | `astra_nexus.cpp` pattern | analytic battery N1–N11 (contract) | **CLOSED (M0) — golden `ad64f810`; the standing v1 oracle** |
| canvas | C++/CUDA/Vulkan | interop skeleton + CINEMATIC post stack | GARGANTUA post chain, Buddhabrot loop | visual gate + perf gate | **CLOSED (M1) — 1M @ 499 fps, CINEMATIC 10/10 (P1 Vulkan owed, D-012)** |
| newton | CUDA | classical tier: hash, PM/P³M, leapfrog | DSA nbody kernels | nexus N2 + conservation | **CLOSED (M2) — 4 goldens (kepler <1e-6 @ 10⁶ ticks)** |
| arrow | CUDA | thermo meters + Ratchet inscription | ORRERY ratchet | nexus N6; DP universality ref | **CLOSED (M3) — 4 goldens (Loschmidt echo exact to 6 decimals)** |
| einstein | CUDA | γ/τ/1PN/photons/2.5PN | astra_nexus formulas | nexus N3/N7/N8/N9 | **CLOSED (M4) — 3 goldens (exact Sommerfeld 0.5%)** |
| gargantua | CUDA(/OptiX) | BH entities, P–W disks, Hawking, lensed view | blackhole.cu (lift) | nexus N4 + analytic Kerr | **CLOSED (M5) — 3 goldens (unscripted collapse; Hawking on the analytic clock)** |
| planck | CUDA (cuFFT) | quantum bubbles + inscription collapse | new | nexus N5 | **CLOSED (M6) — 3 goldens (which-way detection kills fringes 16×)** |
| cosmos | CUDA | torus wrap, little-planet, light history, a(t), roaming bubbles | new | circumnav/expand/bubbles + wrap consistency | **CLOSED (M7) — 21/21 green** |
| **substrate_nexus** | C++17 fp64 | **v2 N0**: spherical Einstein–Klein–Gordon oracle; demonstrates Choptuik Type-II critical collapse (BH mass → 0 at p*) | new (NR literature; G=c=1) | S1–S5 battery (flat/sub/super/transition/massive) | **CLOSED (v2 N0) — golden `13aa73e5`; precise γ deferred to AMR (D-021)** |
| **field_nexus** | CUDA (cuFFT) | **v2 N1**: Schrödinger–Poisson on a 3D lattice — M6 split-step ψ **welded** to M2 PM cuFFT-Poisson (one loop, one field). The PM+ψ gravity weld demonstrated across 3 dynamical regimes (soliton / collapse / two-body) | M6 `kQ3PhaseK` + M2 `kGreen` (verbatim) | N5 (freepacket/sho3d exact) + analytic SP soliton + free-fall/two-body + structural echo | **CLOSED (v2 N1) — 6/6 green: freepacket `03dd3a3b`, sho3d `dfbc6185`, soliton `d163d765` (weld, scale_rel 3e-8, two-passed), echoF `433ddcc8`, cloudF `2308ea49` (collapse@t_ff), mergerF `a09dda6a` (2-body attraction); galaxyF (rotation curve) deferred Q-N1-1 (D-022)** |
| **lapse_nexus** | CUDA (cuFFT) | **v2 N2**: the clock — lapse α(x) = √(1+2Φ/c²) from the substrate potential + a declared proper-time field τ = ∫α dt. Gravitational time dilation / redshift, on the substrate (bridges M2 Newtonian Φ → M4 relativistic time) | analytic Schwarzschild z = 1/√(1−r_s/r)−1 (exact) + the N1/M2 PM solver (weld) | **CLOSED (v2 N2) — 2/2 green: redshift `e2c75be5` (exact Schwarzschild; a clock at 2r_s ticks 40% slow; α-err 5.7e-6), redshiftPM `3dddb950` (the weld: PM well of Newtonian depth, 3.6%). Temporal metric only; light-bending/precession → N3 (D-023)** |
| **curve_nexus** | C++17 fp64 | **v2 N3**: geometry curves — geodesics through the substrate's weak-field metric bend light + precess orbits at exact GR. The 1919 factor of 2 DECOMPOSED: the lapse (time, N2) bends 2GM/bc², the spatial curvature (N3) adds the other half → 4GM/bc² | analytic GR: deflection 4GM/bc² (+ 2GM/bc² lapse-half), precession 6πGM/(c²a(1−e²)) | **CLOSED (v2 N3) — 3/3 green = ALL FOUR classical tests of GR from the substrate metric (with N2 redshift): deflect `4e6c33ca` (light bends 4GM/bc² to 0.3%; space doubles it, full/lapse=2.004), precess `67272705` (0.52%), shapiro `20bfd4d2` (time delay 0.33%). CPU fp64 geodesic oracle; static weak-field (dynamical GPU metric field → N4+) (D-024)** |
| **inspiral_nexus** | C++17 fp64 | **v1 polish**: 2.5PN gravitational-wave inspiral (Peters 1964) — binaries merge on the quadrupole clock; eccentric orbits circularize | Peters 1964 closed forms | analytic: circular merger time (1.3e-13) + a(e) circularization (5e-11) | **CLOSED — 2/2 green: circular `2eba79de` · eccentric `4578d3ac` (D-025)** |
| **precession_nexus** | C++17 fp64 | **v1 polish**: Q-006 resolution — SR / 1PN / combined precession isolation (the π+6π=7π superposition CONFIRMED; the app's 6.41π was a normalization artifact — nominal vs force-distorted semi-latus rectum) | Sommerfeld + GR closed forms | analytic: sr=1.00π · pn1=6.03π · combined=6.95π at the orbit's *actual* p | **CLOSED — 3/3 green: sr `a0e180df` · pn1 `db0818f2` · combined `f9df648f` (D-026; Q-006 RESOLVED)** |
| **fluidcss_nexus** | C++17 fp64 | **Axis-C crown, COMPLETE**: the radiation-fluid critical exponent — the TRUE Evans–Coleman CSS background (V₀\*=0.1124394, one V-zero, N̄'(sonic)=−0.355699) + the relevant perturbation eigenvalue | HKA gr-qc/9607010 §IV/§V; lit κ₀=2.8105525488, β=0.35580192 | two analytic gauge controls (sonic-gauge −N̄'_ss, origin-gauge 1) + G-ANCHOR/G-CONVERGE/G-UNIQUE | **CLOSED (v1.0.0, D-032) — 2/2 green: stageA `27af7920` · stageB `9f8587fd`; **β = 0.3557988 MEASURED** (|Δβ|=3.1e-6). The v0.9.x golden `b4f4e463` banked FRIEDMANN-as-EC — superseded (NOTE.md)** |
| **interop** | CUDA + Vulkan | **R0 (Axis B unparked, D-034)**: zero-copy presentation — shared `VkBuffer` (external memory) written by CUDA, presented by Vulkan, ordered by exported timeline semaphores; LUID-matched device; live window + gated headless face | D-002 Path A · D-012 P1 · CUDA interop samples pattern | validation layers (0 errors/warnings gated) + byte-identity round-trip + monotonic-counter sync | **CLOSED (v1.0.0) — 1/1 green: `interop_r0` `4ba7fbcb` (240 frames ~1.2 s); windowed face live. R1 `cinematic` next** |
| **cinematic** | CUDA + Vulkan | **R1 (the CINEMATIC law binds)**: full §1 post-chain on the R0 path — fp32 HDR → mip bloom → energy-weighted EV auto-exposure → AgX/ACES → astro-stretch+grain (cinematic mode) → triangular dither → single sRGB → BGRA. SUPERNOVA scene (blackbody starfield + ×10⁵ flash) | CINEMATIC.md · GARGANTUA blackbody · minimal-AgX · CoD mip bloom | KAT selftest + G-RANGE ≥ 2^10 + G-EXPOSURE (adaptation moves ≥ 1.5 EV) + validation-clean | **CLOSED (v1.0.0) — 1/1 green: `cinematic_r1` `4962558c`; flash/decay shots delivered for the judged §7 half. R2 (substrate wiring) next** |
| **choptuik_nexus** | C++17 fp64 | **γ crown, direct route**: Choptuik mass scaling M_BH ∝ (p−p\*)^γ on the polar-areal massless-scalar evolver (C++ port of the session-6 research pipeline; sub-cell crossing masses M70/M65 + freeze mass) — **campaign γ = 0.37 ± 0.02** (N=3200); the fine structure DETECTED (wiggle A≈0.15–0.20) | the committed Python table (`gamma_scaling.npy`) — port measured **BIT-EXACT** — + lit γ=0.374(1) | G-DET two-pass · G-ORACLE (< 1e-6; measured 0.0) · G-GAMMA [0.30,0.45] · G-CUT · G-FATE | **CLOSED (v1.0.0) — 2/2 green: scaling `86c68cf9` (γ[M70]=0.3406859239 frozen) · cross `0e04f941`. Δ + γ-to-±0.001 = AMR-gated (D-021 measured TRIPLE, session 6)** |

Shared infrastructure (not a module): `core/lib/` — the liborrery lift (D-005), KAT-selftested, verbatim.

**v2 SUBSTRATE track** (D-020; `docs/PROPOSAL_2026-07-12_v2_substrate.md`): the one-field rewrite (Schrödinger–Poisson → lapse → Einstein–Klein–Gordon → closures + Ratchet lattice). N0 `substrate_nexus` is its standing fp64 oracle, built CPU-first (no GPU); the GPU ladder is climbing — **N1 `field` (the SP gravity weld, 6/6), N2 `lapse` (the clock, 2/2), and N3 `curve` (geometry, 3/3 — with N2, ALL FOUR classical tests of GR pass on the substrate) are CLOSED** — with N4 `horizon`/`star` (GPU; the AMR/Choptuik crown per D-028) ahead, each gated against N0 and the 21 v1 goldens (the oracle farm).

## 12 · What this is NOT

- Not an ORRERY tool: it opens a window; ORRERY's invariants forbid that. Sibling, same doctrine.
- Not the science, and not a proof of any theory: TU renders mechanisms (measured, gated); interpretation lives elsewhere.
- Not full numerical relativity (weak-field + kinematics + pseudo-potentials — the right scope, stated) and not many-body QM (bubbles are single-particle ψ; entanglement beyond pairs is out of scope v1 and said so on the HUD in physical mode).
- Not a place where beauty and honesty trade off silently: the toggle makes the difference visible, always.

*The spec is the product. Freeze the dials, gate the goldens, and let the universe compound.*
