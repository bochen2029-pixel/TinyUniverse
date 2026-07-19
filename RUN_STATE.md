# RUN_STATE.md

## SESSION 4 (2026-07-19, autonomous) — **THE β WALL FALLS: β = 0.3558019 MEASURED (D-032)**

- **The three-session wall was the BACKGROUND.** The banked "Evans–Coleman" background (Stage-A golden `b4f4e463`, D-027) is actually the **collapsing flat Friedmann solution** — measured: V≡−√(1−1/A) to 1.9e-10; A=1+(2−γ)ω, oi\*=3/8, 2m/r=1/3 are all FRW radiation identities. `hka_ec` selected it by construction (sonic criterion "V=−c_s" = precisely the Friedmann point of the §IV V₀-parametrized sonic line). The gauge gate is **background-blind** — the true bottom of D-030's "necessary but not sufficient". The 0.35699 footnote was the true background's fingerprint (N̄'_EC(sonic)=−0.355699; paper typo), not a red herring — corrects D-031.
- **True EC background built** (`nr_ec2.py`, HKA's own construction): **V₀\* = 0.112439401388**, one V-zero ✓, center-relaunch closes to 6 digits ✓.
- **Lyapunov evolution built + validated** (`nr_lyap.py`: gauge control 1.000047 exact; fifth-equation formulation B; one-step-map spectrum). On the true EC background: **κ₀ = 2.8105526(3) (ref 2.8105525488) ⇒ β = 1/κ = 0.3558019 (ref 0.35580192)** — G-CONVERGE ✓ (N=200/300/400 + edge-standoff sweep) · G-UNIQUE ✓ (one relevant mode; gauge → 1) · **G-ANCHOR ✓ (|Δβ|≈1e-8 ≪ 4e-3)**.
- **DONE (same session): the full close.** (a) **Redundant recovery** — the shoot on the EC background (`nr_shoot_ec.py`): κ₀=**2.810552374**, β=0.355801945, PLUS both gauge controls land analytically (sonic-gauge 0.355698925 = −N̄'_ss; origin-gauge dip at exactly 1.00) — the "walled" shoot machinery was never broken, it was fed Friedmann. (b) **C++ port + goldens** — `fluidcss_nexus.cpp` **v1.0.0**: Stage-A = the true EC construction (V₀\*=0.1124394014 matches Python to 10 digits; fingerprint −0.355699; closure 4.2e-6) golden `27af7920` (supersedes `b4f4e463`, NOTE.md); **Stage-B = the eigenvalue: κ₀=2.810577211 → β=0.355798800, golden `9f8587fd`** — G-ANCHOR (|Δβ|=3.1e-6 ≪ 4e-3) + G-CONVERGE (1.4e-8) + G-UNIQUE (1 root in (1.1,3.6)) + both controls PASS; both goldens two-passed; harness row `fluidcss_stageB` added — **full suite re-run cold same session: ALL GREEN, 40 goldens, 686 s** (`runs/verify_40_2026-07-19.log`). Contract v1.0.0 FROZEN (original gates unchanged). **The fluid-β crown is CLOSED.**
- **NEXT: the scalar-DSS γ=0.374 crown — `contracts/similarity_nexus.contract.md` v0.1.0 DRAFTED (this session) — ⏸ OPERATOR REVIEW REQUIRED before any source** (contract-first hard rule). Grounded in Gundlach gr-qc/9604019 read verbatim (§2.2–2.3 DSS-as-eigenvalue-problem, §4.1 Floquet perturbations; source fetched to `tournament/gamma/GUN96_src/`, gitignored like HKA99_src): Stage-A = DSS background on the (ζ, τ-periodic) cylinder, Δ eigenvalue (lit 3.4453±0.0005); Stage-B = Floquet λ₁=−2.674 → γ=−1/λ₁=0.374±0.001; analytic controls = the λ=0 τ-translation gauge mode + the Δ/2-impostor check (the D-032 "which solution did you find" lesson institutionalized) + Δ-vs-N0-echo redundant recovery (D-028). Then per the operator queue: renderer (Axis B) / AMR-N4 / v1 polish tail.

## QC PASS (2026-07-19) — full-repo verification + reconciliation

- **Full harness cold: 39/39 GREEN in 1179 s** (`verify.py --build`; all 9 binaries compiled clean; contended card; log `runs/qc_harness_2026-07-19.log`). **Reality = claimed state.**
- **β-thread session-3 claims re-verified cold:** `nr_rederive` all 16 coeffs diff=0 · gauge gate 1e-10 · `nr_laurent` μ=1−2κ exact + Frobenius gate 9.8e-14 · `nr_sonic` vs integrated EC 1.5e-10 @ t=−0.02 · `nr_shoot` gauge positive control κ=0.99998705 (≈1 to 1.3e-5; the recorded 0.99999924 differs in last digits — env/tolerance-sensitive, same mode).
- **Reconciled/remediated this pass:** `RESULTS_hka_beta.md` ended on the FALSIFIED D-030 conclusion → **D-031 addendum appended** (operator PROVEN; wall = extraction method; next = Lyapunov §V.G). ARCHITECTURE §11 refreshed (8 stale rows CLOSED'd; inspiral/precession/fluidcss rows added). 5 contract headers BUILDING/DRAFT → FROZEN (curve title also v1.0.0→v1.1.0). MODULE.md sections added for field/lapse/curve/fluidcss (substrate) + inspiral/precession (nexus). README 22/22 → 39/39 + N1–N3/crown rows. Harness budget "< 5 min" → measured ~20 min (BUILD.md + harness/README, amended). Orphan `goldens/fluidcss` (the superseded v0.9.0 golden, gated by nothing) removed. 10 pre-gitignore `.pkl` caches (~51 MB) untracked from the index (files kept on disk; `phase4/.gitignore` *.pkl now effective). `nr_sonic`/`nr_laurent` `__main__` validation prints made stranger-safe (gate inside the series radius; out-of-radius rows labeled expected).
- **Push decision EXECUTED (operator directive 2026-07-19: "excise HKA99_src, then push public").** The 23 unpushed commits were rewritten (`git filter-branch --index-filter`, range `92ce270..master` only — pushed history untouched) to remove `tournament/gamma/phase4/HKA99_src/` (arXiv source of gr-qc/9607010, third-party ©), verified zero HKA99_src objects in the outgoing set, then pushed. **Old→new hash map: `docs/HASHMAP_2026-07-19.md`** — docs/DECISIONS/memory written before this date cite the OLD hashes; they resolve locally via the **LOCAL-ONLY branch `backup/pre-excise-2026-07-19` (NEVER push it — it still contains the paper source)**. The paper files stay on disk, now gitignored (`phase4/.gitignore`); strangers re-fetch from arXiv (command in the hashmap doc). (`bg_h1e4.pkl`, 40 MB, was already in pushed public history — moot.)
- **NEXT unchanged:** implement HKA's Lyapunov time-evolution (§V.G) per `docs/CONTINUATION_2026-07-13_session3.md` §5a → κ → β → C++ Stage-B port + golden.

**As of:** 2026-07-12 · **Milestone:** v2 **N3 `curve` — geometry curves: light bends + orbits precess at exact GR, on the substrate** (3/3 goldens GREEN — all four classical tests of GR, on master). **State:** v1 complete (M0–M7); v2 ladder N0 (oracle) + N1 `field` (SP weld, 6/6) + N2 `lapse` (the clock, 2/2) + **N3 `curve` (geometry, 3/3) all CLOSED on master**. N3 integrates geodesics through the substrate's weak-field metric ds²=−(1+2Φ/c²)c²dt²+(1−2Φ/c²)dl²: `curve_deflect` `4e6c33ca` (**light bends at exact GR 4GM/bc²** — the 1919 factor of 2 DECOMPOSED: lapse/time [N2] = one half, spatial curvature [N3] = the other; full/lapse=2.004) + `curve_precess` `67272705` (**perihelion precession** 0.72°/orbit, 0.52%) + `curve_shapiro` `20bfd4d2` (**Shapiro time delay** 0.33%). **With N2's redshift, the substrate now passes ALL FOUR classical tests of GR from its own metric** (redshift · light-bending · precession · Shapiro delay). CPU fp64 geodesic oracle (no GPU, like N0). Tool `substrate/curve_nexus.cpp`; contract v1.1.0; D-024. Honest boundary: static weak-field geodesic oracle — a dynamical GPU metric field with back-reaction is N4+. **Next = operator's call: N4 `horizon` per the ROADMAP is the Choptuik critical-exponent crown, which is the documented AMR research wall (D-021, and the overnight fluid-β/scalar-γ walls) — NOT a clean one-pass rung. The weak-field GR arc (N1 field → N2 lapse → N3 curve) is COMPLETE; genuine next options: the AMR contract (unlocks the crown), the renderer (Axis B), v1 polish, or the fluid-β/scalar-γ crown research.** — *(prior N2 milestone line, superseded:)* the clock: gravitational time dilation, on the substrate (2/2 goldens GREEN, on master). **State:** v1 complete (M0–M7); v2 N0 CLOSED; **N1 `field` CLOSED** (merged to master `8ace261`); **N2 `lapse` CLOSED** — the substrate now has a clock. N2 turns N1's Newtonian potential Φ into a per-cell lapse α(x)=√(1+2Φ/c²) and a declared, hashed proper-time field τ(x)=∫α dt (c goes LIVE). Goldens: `lapse_redshift` `e2c75be5` (**EXACT Schwarzschild redshift** — a clock at r≈2r_s ticks 40% slow, z=0.398, α-err 5.7e-6 vs 1/√(1−r_s/r)−1) + `lapse_redshiftPM` `3dddb950` (**the weld** — the substrate's own PM-Poisson gravity makes a well of Newtonian depth 3.6% → correct gravitational time dilation). Both two-passed cold. Tool `substrate/lapse_nexus.cu`; contract `lapse.contract.md` v1.0.0; D-023. Honest boundary: temporal metric only (redshift exact; light-bending/precession → N3). *Reconciliation note: RUN_STATE/TASKLIST had lagged git — N1 was merged to master at `8ace261`; both are now corrected to ground truth.*

## CURRENT DIRECTIVE (2026-07-12→13, this session)

**Axis C — the CSS×spectral crown.** Reconciled in **D-028** (AMR-for-N0 graveyarded → N4-GPU; γ route = eigenvalue γ=1/Re λ₀; D-021 amended to gauge×discretization). Warm-up first (the ruling): land the known **β=0.35580192** by closing `fluidcss_nexus` **Stage-B**.

**PROGRESS — SESSION 3 (2026-07-13): sonic Frobenius EXACT; β wall MOVED to the operator (D-030).**
- ✅ **The sonic analytic machinery is REBUILT + EXACT** (the §5a step-1 deliverable). `nr_sonic.py` = analytic background sonic Taylor series (order-1 solvability is *quadratic*, EC branch via the desingularized eigenvector; validated reference-free ODE-residual ~1e-15, matches integrated bg to 3.7e-11). `nr_laurent.py` = exact operator Laurent, **residue μ=1−2κ EXACTLY** (reconciles the −4.48/+10.5/−4.62 disagreement); 3 analytic modes solve the direct ODE to **9.8e-14** — the Frobenius GATE PASSES. Fixes the polyfit root-inaccuracy the prior session named.
- ✅ **The OPERATOR is PROVEN correct from first principles (D-031 — corrects D-030).** `nr_rederive.py` symbolically linearizes the primary EOM (`rflanl.tex` 4.1/4.2) keeping ∂ₛ→κ; all 16 fluid coeffs match §V (`hka_pert_hka99`) EXACTLY (diff=0); metric rows verified by hand. Background = true EC solution (N∝e^{−x}, A=1+(2−γ)ω exact, two methods). The 0.35699 gauge-mode footnote is a red herring (needs ω0<0). So D-030's "operator error" is FALSIFIED — operator+background+identity all correct.
- 🔴 **The wall is the eigenvalue-EXTRACTION method.** HKA's own sonic→center shoot (`nr_shoot.py`, verified: resolves the gauge mode to **κ=0.99999924**) surfaces ONLY the gauge mode κ=1 across κ∈[−1.5,12]+complex — NOT the physical κ=2.81. Since the operator is now proven correct, κ=2.81 must be extractable; the shoot's failure to surface it is a known-hard-numerics issue, not an operator issue.
- **NEXT: implement HKA's SECOND method — the Lyapunov time-evolution** (`rflanl.tex` §V.G, "very regular": M_s non-singular ⇒ NO sonic singularity in the s-evolution `∂ₛ(ω̄,V)=M_s,fl⁻¹[Gmat_fl·Ψ − M_x,fl ∂ₓ(ω̄,V)]`, metric from spatial constraints). Power-iterate → dominant growth rate = relevant κ=2.81. Then C++ port + β golden + `fluidcss_nexus` v1.0.0. **β NOT measured, none faked (D-016/D-021).**

## OPERATOR WORK QUEUE (priority order; #1 corrected per D-028)

1. **#3 · v1 polish** ← **IN PROGRESS.** a(t) little-planet timelapse video · GARGANTUA Kerr geodesic art pass · **2.5PN gravitational-wave inspiral** ✅ **DONE** (`nexus/inspiral_nexus.cpp`, Peters 1964 — circular merger time to 1.3e-13, eccentric a(e) circularization to 5e-11; goldens `2eba79de`/`4578d3ac`; D-025) · Q-006 ✅ **DONE** (`nexus/precession_nexus.cpp`, D-026 — the **7π superposition is CORRECT**; the app's 6.41π was a normalization artifact, nominal vs force-distorted p; goldens `a0e180df`/`db0818f2`/`f9df648f`). **#3 status:** 2.5PN ✅ + Q-006 ✅ DONE; **a(t) little-planet timelapse = CINEMATIC WALL** (the `expand` box is dim/near-uniform, ~2× growth over a=1→1.93 too subtle to read — physics-honest but not beautiful; frames `runs/expand_ts*.png`). Operator chose **option 2: the galaxy little-planet orbit** as the beauty win (`runs/galaxy_orbit.gif` — 24-frame 360° az orbit, delivered; `--nohud` render flag added to the app). **#3 remaining: GARGANTUA Kerr art pass.**
2. **#4 · crown research (the CSS×spectral γ, D-028)** ← **ACTIVE.** Stage-A ✅ BANKED (`substrate/fluidcss_nexus.cpp` v0.9.1, golden `fluidcss_stageA` `b4f4e463`: Evans–Coleman background — oi*=3/8, sonic point + 2m/r=1/3 EXACT; D-027). **NOW: close Stage-B** (β=0.35580192 warm-up — fix the perturbation κ-coupling K so the gauge-mode gate passes, then the validated QR shoot reads β; discard gauge modes κ≈0.357/1). **THEN:** scalar-γ = 0.374 (DSS Floquet, Gundlach gr-qc/9604019) → new `similarity_nexus` contract.
3. **#2 · renderer (Axis B)** — R0 `interop` Vulkan⇄CUDA live window (NOTE: Vulkan SDK not installed — measured gate).
4. **#1 · AMR contract → N4-GPU (reinstated, D-028)** — AMR-for-N0 graveyarded on the ≤5-min budget; its recipe (`tournament/gamma/phase1/amr-berger-oliger.md`, fully specified) is the spec for the future **N4 `horizon` GPU** stage where the 2^L subcycling is affordable. Not N0's γ route.

---

**Prior (master):** v2 **N0 `substrate_nexus` CLOSED** — the substrate ladder's foundation stone. v1 complete (M0–M7); v2 track begun (D-020). **22/22 goldens GREEN** (~121 s): nexus + substrate_nexus (CPU fp64, no GPU) + 20 GPU. Harness runs CPU oracles under any card contention.

## What N0 established (v2 foundation stone; D-020/D-021)

- **A spherical Einstein–Klein–Gordon solver** (`substrate/substrate_nexus.cpp`, CPU fp64, polar-areal constrained evolution, G=c=1) — the standing oracle for the v2 GPU ladder, built **before any GPU code** and needing no card. Golden `13aa73e5` (N=800, r_max=24, ~20 s).
- **Four rigorous oracles** a stranger runs cold: flat-space mass conservation (7.9e-4) · subcritical dispersal · supercritical horizon formation (2m/r→0.98, lapse→3.9e-5) · massive-KG stability+conservation (3.3e-6).
- **The Choptuik phenomenon, honestly**: S4 demonstrates the **Type-II critical transition** — a bracketed critical amplitude p*=0.404 where the black-hole mass falls continuously to zero (M_BH 0.412→0.118, ratio 0.29; arbitrarily small black holes), clean far-field power law R²=0.998. **The precise universal exponent γ ≈ 0.374 is DEFERRED to the AMR contract (D-021)**: a uniform grid caps the self-similar curvature (measured), refining turns near-critical chaotic — Choptuik needed AMR. Per RAYFORMER/D-016, N0 ships what it proves (the transition) and names what it cannot (the exponent) rather than fake a number that would poison the oracle farm.

## Next — the v2 GPU ladder, or v1 polish (operator's call)

- **v2 continues:** N1 `field` (Schrödinger–Poisson 256³–512³ GPU, the two v1 engines fused) → N2 `lapse` → N3 `horizon` (BSSN+scalar; may reach precise γ) → N4 `star` (the hydrogen-ball sentence). Each gated against N0 + the 21 v1 goldens (the oracle farm). Contract-first per module.
- **v1 polish backlog:** a(t) little-planet cosmology timelapse video · Φ-coupled bubbles · P1 Vulkan/ImGui · GARGANTUA Kerr geodesic art pass · 2.5PN · Q-006.
- **AMR contract** (would let N0/N3 reach the precise Choptuik γ) — named future work.

## Chores carried

clang/g++ parity (nexus + substrate_nexus) · P1 Vulkan · ImGui · TAA · art pass (Kerr view + OptiX → ORRERY `lens`) · cufft64 dll packaging · P³M/spatial hash · Q-006 · 2.5PN · bound oscillaton (N0 S5 deferred).

## Standing context

v1: the regime ladder + the world, complete (M0–M7, 21 goldens). v2: the substrate rewrite begun — N0 the fp64 oracle spine (22nd golden). Every physics claim a stranger reproduces cold from BUILD.md; every claim the tool can't back is named, not faked (D-016/D-021). GPU shared (preflight exits 3 <2 GB; CPU oracles run regardless). Repo docs authoritative over agent memories.

## What M7 established

- **Topology is live and exact.** Space is a 3-torus in coordinates, not just in the force field: `kDriftK` wraps positions into `[−256,256)` (the ±512 shift is a power-of-two, bit-reversible — the Loschmidt echo stayed byte-identical through it, D-019), minimum-image separations everywhere. `circumnav` `3f18f02c`: a photon laps the universe twice and comes home to 3.1e-5 su; photons are now **structurally ageless** (τ≡0) after an fp32 clamp-leak fix the gate caught.
- **Cosmology is a frozen number.** `expand` `ce448f2b`: Einstein–de Sitter comoving expansion on the existing PM solver (Green ×1/a, drift ×1/a², Zel'dovich growing mode) reproduces **linear growth D ∝ a** — slope 0.988, amplitude ratio 0.7% at a_end=1.927.
- **Quantum went roaming.** `bubbles` `78b753f1`: 64³ batched split-step bubbles spawn on PM-isolation (8/8), spread at the free-packet σ(t) to **0.57%**, and collapse by M3 Ratchet intrusion (4 collapsed / 4 alive). RECORDED particles stay classical.
- **The namesake shot shipped.** Stereographic little-planet projection + 27-image torus splat + light-history retarded-time sampling + bubbleVis ψ-glow, all non-declared, CINEMATIC §7 10/10. `runs/cosmos_littleplanet.png` (the globe), `runs/cosmos_planet_lensed.png` (a BH's Einstein ring on the little planet), `runs/cosmos_bubbles.png`. Perf: **225 fps avg / 58 min** windowed 1M with the 27-image splat + history.

## Next — operator's call (both doors open, one directive away)

- **v2 SUBSTRATE** (`docs/PROPOSAL_2026-07-12_v2_substrate.md`, operator-review pending): the one-field rewrite. First act on "go substrate": `contracts/substrate_nexus.contract.md` + build N0 (CPU fp64 spherical EKG, **Choptuik γ** golden — no GPU needed).
- **v1 polish backlog:** a(t) little-planet cosmology timelapse video · Φ-coupled bubbles · P1 Vulkan/ImGui · GARGANTUA per-pixel Kerr geodesic art pass · 2.5PN · Q-006.

## Chores carried (from M0–M6)

P1 Vulkan (SDK) · ImGui · TAA · clang/g++ nexus parity · art pass (Kerr geodesic view + OptiX spike → ORRERY `lens`) · cufft64 dll packaging · P³M/spatial hash · Q-006 · 2.5PN.

## Standing context

The regime ladder is COMPLETE (beauty → gravity → arrow → relativity → black holes → quantum) **and the world is complete** (wrapped, expanding, self-visible, roaming-quantum). Every physics claim a stranger reproduces cold from BUILD.md. 21 goldens, 8 contracts, ~105 s harness. GPU is shared (preflight exits 3 if <2 GB free; timings on a contended card not citable). Repo docs authoritative over agent memories.

## What M6 established

- **The ψ engine** (2D 256² split-step cuFFT, real + imaginary time, absorbing edges, host-fp64 observables) and **collapse as counter-keyed sampling** — the Tonomura experiment computed once, sampled 4096 times.
- **The measurement problem is a golden** (`doubleslit` `47a67d66`): fringe contrast 0.83 unobserved; the which-way detector (M3 record semantics) collapses it to 0.052. Interference emerges one particle at a time and dies when watched.
- **Oracle exactness**: tunneling T vs same-grid fp64 oracle to 1e-6 (`f1e7a061`); SHO ground state E₀ = ħω to 6 digits, σ exact (`fa2e009e`). D-018 findings: observable-vs-evolution coefficient bug caught by σ-exactness; oracle-isolates-implementation principle; near-field fringe shift is real physics.
- Q-004 resolved: ~5 MB per 64³ bubble — hundreds coexist with the full universe.
- Goldens now 18 (new: doubleslit `47a67d66` · tunneling `f1e7a061` · shoq `fa2e009e`). Contracts: + planck v1.0.0. App v0.4.x: 17 scenarios.

## Next task — M7 `cosmos` (the tiny planet)

Contract first (`contracts/cosmos.contract.md`), then:
1. **3-torus wrap live** (forces already periodic via PM — positions/light wrap next); light-history ring buffer ("see your own past" — ASTRA-7 Kepler-at-t_emit machinery, nexus N9-proven).
2. **Stereographic little-planet projection** — the product's namesake visual.
3. Scale factor a(t) cosmology mode (hot start → structure → heat death timelapse).
4. Roaming 3D bubbles integration (deferred from M6): spawn-on-isolation in the live universe, live Ratchet-collapse, bubbleVis.

**M7 gate:** the signature screenshot — the universe as a globe, wrapped and lensed, self-visible; cosmology timelapse runs; goldens + harness green.

## Chores carried

P1 Vulkan (SDK) · ImGui · TAA · clang/g++ nexus parity · art pass (Kerr geodesic view + OptiX spike → ORRERY `lens`; doubleslit fringe display saturation; collapse-scenario framing) · cufft64 dll packaging · P³M/spatial hash · Q-006 · 2.5PN · live bubble-budget probe.

## Standing context

The regime ladder is COMPLETE: beauty → gravity → thermodynamics/inscription → relativity → black holes → **quantum mechanics** — every rung golden-gated, every physics claim a stranger can reproduce cold from BUILD.md. 18 goldens, 7 contracts, ~199 s harness. GPU is shared (preflight in verify.py; timings on a contended card are not citable). Repo docs authoritative over agent memories.
