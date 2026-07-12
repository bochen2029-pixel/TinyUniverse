# TINY UNIVERSE — THE AMBITION & THE ROADMAP
*The maximal plan. Written 2026-07-12. Grounded in the doctrine (contract-first, deterministic-or-it-doesn't-ship, golden-gated, honest-oracle). This is the north star and the path; `ARCHITECTURE.md` is the spec, `TASKLIST.md` is the tactical queue, this is the arc that binds them.*

---

## 0. THE NORTH STAR

**One native Windows binary that is simultaneously (a) a cinematic, real-time, interactive universe you can fly through, and (b) a deterministic physics instrument where every regime — quantum → Newtonian → relativistic → black-hole → cosmological — *emerges* from a single substrate under retuned constants (ħ↑, c↓, G↑), and every claim is reproducible cold from a seed + input-trace.**

The scale trick is the whole soul: by retuning the constants so the quantum, gravitational, and relativistic scales collide at human-observable particle counts, phenomena that in the real universe are separated by 40 orders of magnitude become *visible in the same window at the same time*. You watch decoherence, gravitational collapse, horizon formation, Hawking glow, and cosmic expansion as one continuous, beautiful, honest thing.

**The definition of MAXIMAL done:** a stranger clones the repo, reads `BUILD.md`, and in one pass produces a windowed universe that (1) passes every physics oracle gate cold, (2) passes the `CINEMATIC.md` beauty gate, (3) runs interactively at 60 fps on an RTX 4070 Ti SUPER, (4) reproduces at least one published critical exponent (Choptuik β or γ) on the scale-compressed lattice, and (5) is byte-identical in its declared state across runs. Beauty AND honesty, both load-bearing, both provable.

*"The playable sibling of ORRERY."* ORRERY is the deterministic instrument; TINY UNIVERSE is the instrument you can also *inhabit*.

---

## 1. THE THREE AXES (and the thesis that binds them)

The project advances on three axes that must eventually converge into the single binary:

- **Axis A — SUBSTRATE DEPTH** (the physics). The v2 program: collapse the many hard-coded v1 modules into *one field, one clock, one lattice* from which every regime emerges rather than being separately programmed. This is the intellectual core.
- **Axis B — BEAUTY** (the renderer). The parked-but-inevitable Vulkan⇄CUDA / RTX real-time window that makes the substrate *visible* and *gorgeous*. Non-negotiable for the north star; currently the largest gap.
- **Axis C — RIGOR** (the crowns). Validating the whole methodology on analytic ground truth: the critical-collapse eigenvalues (fluid β=0.35580192, scalar γ=0.374), the ORRERY evidence-grade coupling, and the RAYFORMER honest-baseline discipline for any isomorphism/speedup claim.

**The thesis:** these are not three projects. A black hole in the final demo is (A) an emergent horizon in the substrate, (B) a path-traced Kerr-shadow with a photon ring, and (C) a configuration whose formation threshold reproduces the published critical exponent. When the same object satisfies all three gates at once, the universe is *real* — beautiful, interactive, and true. Everything below is in service of that single convergence.

---

## 2. HONEST CURRENT STATE (2026-07-12)

**Axis A — Substrate:**
- v1 physics oracles **M0–M7 COMPLETE**, 22 goldens green (`nexus, canvas, newton, arrow, einstein, gargantua, planck, cosmos`) — quantum bubbles, thermo/Ratchet decoherence, relativity always-on, Kerr black holes, quantum mechanics, torus cosmos. These are the *reference regimes* the substrate must reproduce.
- v2 substrate **N0 `substrate_nexus`** (Choptuik transition oracle) **DONE**; **N1 `field`** (Schrödinger–Poisson: M2 PM-Poisson ⊕ M6 split-step ψ fused into one deterministic CUDA loop) **DONE — 6/6 goldens, merged to `master` @ `0eefafa`**. The weld is demonstrated across three dynamical regimes (bound soliton / classical collapse / two-body merger) and both Madelung limits.

**Axis B — Renderer:** `app/MODULE.md` skeleton exists; the full Vulkan⇄CUDA interop + CINEMATIC + RTX real-time window is **PARKED**. This is the biggest single gap to the north star.

**Axis C — Rigor:**
- **fluid-β:** Stage A (the Evans–Coleman self-similar background) **LANDED at machine precision** (σᵢ*=3/8, sonic point exact, 2m/r=1/3). Stage B (the perturbation eigenvalue) **WALLED but massively advanced this session**: the operator is now *gauge-exact with correct indicial structure at both boundaries*, and a *validated compound-matrix Evans-function solver* exists — but the operator is missing the true ∂ₛ content of the fluid PDEs, so it carries no physical mode yet. Precisely diagnosed; no number faked.
- **scalar-γ:** not started (the harder DSS crown).
- **ORRERY / RAYFORMER discipline:** coupled and honored (D-021 gauge×discretization ruling; honest walls logged, never tuned).

**Reproducibility:** everything on GitHub (`github.com/bochen2029-pixel/TinyUniverse`), `master` @ `0eefafa`, tournament branch @ `e289922`.

---

## 3. AXIS A — THE SUBSTRATE LADDER (the detailed N-roadmap)

The v2 vision — **one field ψ, one clock (lapse α), one lattice (metric a²)** — realized as a rung-by-rung ladder. Each rung is contract-first, names its analytic oracle, freezes a golden, and makes ONE emergence claim. The deep bet: by the top of the ladder the *separate* v1 modules are re-derived as *limits* of a single running substrate.

| rung | name | the physics it adds | oracle (gate) | emergence claim |
|---|---|---|---|---|
| N0 | `substrate_nexus` ✅ | the deterministic substrate shell + Choptuik transition detector | analytic transition | "there is a critical threshold" |
| N1 | `field` ✅ | ψ on a 3D lattice, self-gravitating (Schrödinger–Poisson) | free-packet σ(t), SHO E₀, soliton scale-covariance, free-fall collapse | "one field carries quantum AND classical gravity" |
| **N2** | **`lapse`** | **the clock: proper time τ per cell set by the lapse α; redshift; time dilation** | Schwarzschild-lite redshift z=1/√(1−2GM/rc²)−1 near a lattice mass; clock-rate ratio | **"clocks run slow in gravity wells, on the substrate"** |
| N3 | `curve` | the lattice metric a² back-reacts on energy density (GR-lite constraint) | the fluid-CSS background (Stage A!) / TOV-lite equilibrium; Misner–Sharp mass | "geometry curves; light bends; orbits precess" |
| N4 | `horizon` | emergent trapped surfaces when collapse × curvature cross threshold | **Choptuik critical exponent (β fluid / γ scalar) — the crowns** | "collapse *makes* horizons; the threshold is universal" |
| N5 | `hawking` | semiclassical field modes near the emergent horizon → thermal spectrum | Hawking T_H(M) scaling; Bogoliubov/greybody | "black holes glow and evaporate" |
| N6 | `expanse` | the substrate on a 3-torus with scale factor a(t); structure from ψ fluctuations | Friedmann (retuned), linear growth factor, Jeans scale | "the universe expands; structure forms" |
| N7 | `unification` | ALL regimes coexisting in ONE running substrate | multi-region consistency: each sub-domain matches its own regime oracle simultaneously | **"it is all one thing"** |

**Design notes / the hard parts, per rung:**
- **N2 `lapse`** — cleanest next rung, fully in hand. The lapse `N ≡ α/(a·eˣ)` already appears in the fluid-CSS formalism; here it becomes a live declared field driving per-cell proper-time integration. Determinism: proper time is a declared, hashed field. Bridges Newtonian potential (M2) → relativistic time (M4) *inside* the substrate. **This is the recommended immediate build.**
- **N3 `curve`** — the fluidcss_nexus Stage-A work plugs in here as the static-collapse oracle. The ambition: a scale-compressed Einstein constraint (Hamiltonian + momentum) solved on the lattice each step, cheaply. Risk: full GR is expensive; the art is the minimal back-reaction that still bends light and precesses orbits honestly.
- **N4 `horizon`** — where Axis C pays off: the substrate's collapse threshold *should reproduce* the Choptuik exponent. Landing β and/or γ here turns "our black holes look right" into "our black holes form with the correct universal scaling." The apparent-horizon finder + mass-scaling measurement is the golden.
- **N5 `hawking`** — the moonshot rung: a scale-compressed black hole that *actually evaporates in real time*, with T_H tracking mass. Semiclassical Bogoliubov on the lattice.
- **N6 `expanse` / N7 `unification`** — the payoff: a single substrate where you pan from a quantum bubble to a forming horizon to cosmic web, each region passing its own oracle. N7's golden is the crown jewel of Axis A.

---

## 4. AXIS B — THE RENDERER (from parked spike to inhabited universe)

The substrate produces fields (|ψ|², ρ, lapse α, curvature). The renderer makes them a place you can be. Lift aggressively from the quarries (GARGANTUA Kerr renderer, Buddhabrot OptiX two-stream, the one-binary-two-faces pattern).

| rung | name | deliverable | lift from |
|---|---|---|---|
| R0 | `interop` | Vulkan⇄CUDA shared buffer; live substrate field → swapchain; a window showing the sim breathe | the parked M1 canvas spike; VK_KHR_external_memory |
| R1 | `cinematic` | the `CINEMATIC.md` checklist as real shaders: volumetric |ψ|², density-color, HDR/bloom/ACES tonemap, depth cueing | CINEMATIC.md gate; ASTRA-7 look |
| R2 | `lensing` | ray-marched gravitational lensing of emergent horizons driven by the SUBSTRATE metric (not a hardcoded Kerr) — shadow, photon ring, accretion glow | `C:\blackhole\blackhole.cu`, GARGANTUA Kerr renderer |
| R3 | `rtx` | RTX Kit: ReSTIR/RTXDI for millions of emissive ψ-lumps; path-traced GI; hardware-RT lensing; DLSS for 4K@60 | RTX Kit, OptiX recipe from Buddhabrot |
| R4 | `interactive` | the camera & controls: free-fly, little-planet projection, light-history (see your past light cone), deterministic time scrub (pause/rewind), constant-dials UI | M7 cosmos little-planet + light-history |
| R5 | `shell` | polish: Slang shaders, optional UE 5.8 TextureShare fallback renderer, presets, input-trace capture/replay, screenshot/video capture | UE-TextureShare, Slang; backlog per ADR |

**Invariant that keeps rigor intact:** visuals are **non-declared**. The renderer reads the substrate; it never writes the declared state. Frame pacing, resolution, RT quality — all swappable without touching a golden. The renderer boundary stays a clean seam (memory: swappable renderer boundary), so RTX Kit vs UE 5.8 is a measured adoption decision, never a lock-in.

---

## 5. AXIS C — RIGOR & THE CROWNS

The crowns validate that the whole approach *measures the right numbers*. They are the difference between "a pretty physics toy" and "a scale-compressed instrument whose emergent phenomena carry the correct universal constants."

- **fluid-β = 0.35580192** *(Stage A ✅, Stage B ~85%)*. Remaining: pin the true ∂ₛ coupling of the fluid perturbation PDEs. Two clean paths, both written up in `RESULTS_hka_beta.md` / `STAGE_B_UPDATE.md`:
  1. the first-principles `∇ₐTᵃᵇ=0` operator (`hka_pert_firstprinc.py`) completed with HKA's (4.4) 2×2 sonic recombination;
  2. the long-paper (gr-qc/9607010) perturbation coefficients transcribed from the actual PDF, resolving the flagged brace.
  Then the *already-validated* Evans-function solver (`hka_beta_evans.py`) reads off the non-gauge zero → β. **The method is done; only the operator's last coupling is open.**
- **scalar-γ = 0.374, Δ=3.4453** *(not started)*. Harder: a τ-periodic (discretely self-similar) hyperbolic BVP + Floquet analysis. The template that just cracked the fluid operator — *derive from the full PDEs, don't transcribe; validate with an amplification-immune solver* — is the playbook to reuse. Fetch G97 (gr-qc/9604019) appendices first.
- **ORRERY coupling** — every citable claim gets an `autotune`/`posit` receipt; internal convergence is weak evidence; two-pass anything public.
- **RAYFORMER discipline (ADR-007)** — the highest-ambition rigor: if any regime turns out to be *isomorphic* to another computation (a genuine surprise the substrate could reveal), or if any GPU path claims a speedup, it gets an honest baseline and a measurement before it is believed. This is where a *publishable* result could live.

---

## 6. THE CONVERGENCE DEMOS (the moonshots — the actual "wow")

These are the deliverables that prove the three axes have merged. Each is a single deterministic scenario (seed + trace) that is simultaneously beautiful and gate-passing.

1. **"THE COLLAPSE."** Fly alongside a quantum ψ-cloud. It collapses past the Choptuik threshold; a horizon blinks into existence (N4, correct critical scaling); you orbit it and watch the precession (N3); it Hawking-glows and evaporates (N5). One substrate, RTX-gorgeous, byte-reproducible. *The single most ambitious artifact — it exercises N1→N5 + R0→R4 + the β/γ crown at once.*
2. **"THE LADDER."** A continuous, unbroken zoom: Planck-scale ψ-bubble → atom-like bound state → Newtonian orbit → relativistic binary → black hole → galaxy → cosmological horizon. Every scale physics-honest; the camera never cuts. Proves the *emergence* thesis viscerally.
3. **"THE TELESCOPE FOR EMERGENCE."** A research mode: live sliders on (ħ, c, G). Watch classicality, gravity, and horizons appear and dissolve as you turn the dials; measure the critical exponents on the fly. The tiny universe as a genuine instrument for studying *how* the classical/relativistic/gravitational worlds emerge — a scientific contribution, not just a demo.

---

## 7. CROSS-CUTTING INVARIANTS (the doctrine that makes it trustworthy)

Non-negotiable, every rung, every axis:
- **Determinism in the declared path.** Fixed timestep; same (dials, seed, input-trace, steps) ⇒ byte-identical declared-state hash (blake2b). Counter-based RNG only (value = f(seed, id, tick)); no stateful randomness; no float atomics in declared reductions; no `--use_fast_math`.
- **Contract-first.** No module source before its reviewed `contracts/<name>.contract.md`. Semver from 1.0.0; breaking change ⇒ MAJOR bump + migration note.
- **Every physics module names its oracle.** Gate against analytic / CPU-double truth. Oracle-gate failure ⇒ run-state SUSPECT, never silent fallback. `tiny_nexus` is oracle-grade.
- **Regimes emerge; nothing hard-switches.** Representation handoffs (wavefunction ↔ point particle) are physics-triggered (Ratchet redundancy threshold), declared, deterministic.
- **Beauty is law.** Any user-facing visual passes the `CINEMATIC.md` checklist before it ships. "Looks fine" is not a gate; the checklist is.
- **Honesty over progress.** A real wall is reported as a wall, with its exact resume state — never tuned away. (This session's fluid-β Stage-B is the model: correct operator + validated method + precise diagnosis, zero faked numbers.)
- **Exit codes 0/1/2**, never conflated. liborrery lifted verbatim (D-005); any divergence is a named fork.

---

## 8. PHASED SEQUENCING (near → mid → long)

**PHASE I — "The substrate stands up" (near-term, weeks).** *Goal: the physics ladder reaches horizons.*
- N2 `lapse` (contract → source → golden). ← **immediate next build.**
- Consolidate fluid-β Stage-A into `fluidcss_nexus.cpp` + freeze a Stage-A golden (bank the landed win).
- N3 `curve` (the GR-lite back-reaction; Stage-A as static oracle).
- *Optional in parallel:* one focused shot at fluid-β Stage-B closure (the (4.4)-recombination).

**PHASE II — "It can be seen" (mid-term).** *Goal: the substrate is visible and beautiful.*
- R0 `interop` + R1 `cinematic` (unpark the renderer; live window; CINEMATIC gate green).
- N4 `horizon` (emergent black holes; apparent-horizon finder; land β and/or γ *here* against the collapse threshold).
- R2 `lensing` (substrate-driven black-hole rendering).

**PHASE III — "It is real" (long-term).** *Goal: gorgeous, interactive, validated.*
- R3 `rtx` + R4 `interactive` (real-time RTX; the playable camera).
- N5 `hawking`, N6 `expanse`.
- The **"Collapse"** convergence demo (Moonshot #1).

**PHASE IV — "It is a universe" (legendary).** *Goal: the north star.*
- N7 `unification` (the one substrate, all regimes).
- Moonshots #2 (Ladder) and #3 (Telescope).
- Both crowns landed; a publishable isomorphism/speedup result under RAYFORMER discipline.

---

## 9. RISK REGISTER (honest walls & mitigations)

| risk | axis | severity | mitigation |
|---|---|---|---|
| fluid-β / scalar-γ eigenvalue ∂ₛ content is genuinely hard to pin down | C | med | method is already validated (Evans); land it *at* N4 against the live collapse threshold, where the physics constrains the coupling independently of the paper transcription |
| ψ can't carry disk rotation without quantized vortices + huge box (Q-N1-1) | A | med | known; use vortex-carrying ψ or a hybrid particle handoff (Ratchet); documented, not faked |
| real-time RTX **and** honest physics at 60 fps may not co-fit on one GPU | A+B | high | renderer is non-declared & swappable; degrade visuals, never physics; DLSS + ReSTIR budget; async compute (physics on compute queue, render on graphics queue) |
| GR-lite back-reaction (N3) too expensive per step | A | med | minimal constraint (conformal-flat / CTT), not full BSSN; measure before believing (ADR-007) |
| scope: three axes, one person | all | high | the ladder is the spine — each rung ships value alone; the renderer can trail; crowns are parallelizable via ORRERY-armed subagents |
| determinism vs GPU nondeterminism (atomics, FP order) | all | med | liborrery reduction shapes only; the golden catches any leak immediately |

---

## 10. AMBITION TIERS (definition of done, escalating)

- **🥉 BRONZE — "The physics is real and unified."** N0–N7 headless, golden-gated. The substrate reproduces every v1 regime as a limit. No renderer required. *This alone is a serious achievement.*
- **🥈 SILVER — "You can see it, beautifully."** + R0–R2. The substrate renders live and passes the CINEMATIC gate; emergent black holes are lensed from the real metric.
- **🥇 GOLD — "The playable, validated universe."** + R3–R4 (RTX real-time interactivity) + at least one crown landed (β or γ). Fly through it at 60 fps; the black holes form with the correct universal scaling.
- **🏆 LEGENDARY — "The north star."** + both crowns + the "Collapse" moonshot + a publishable RAYFORMER-disciplined isomorphism/speedup + the research-instrument "Telescope" mode. A tiny universe that is beautiful, inhabitable, byte-honest, and a genuine contribution to how we understand emergence.

---

### THE ONE-LINE MANDATE
*Build one rung right. Freeze its golden. Make it beautiful. Prove its number. Let the universe compound — until you can fly through a quantum bubble as it collapses into a correctly-scaling black hole and watch it evaporate, all in one deterministic, gorgeous, honest window.*
