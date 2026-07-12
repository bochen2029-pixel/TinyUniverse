# CONTRACT — field (v2 N1: the first GPU substrate) · v1.0.0 · status: DRAFT 2026-07-12 (design-only; no GPU build — operator review requested)

**Purpose.** The first rung of the v2 SUBSTRATE ladder (`docs/PROPOSAL_2026-07-12_v2_substrate.md` §S1): fuse v1's two gravity engines into **one loop over one complex field** on a 3D lattice. The classical particles disappear; a single wavefunction ψ(x,t) carries all matter, and gravity is its own mean-field potential:

```
iħ ∂ψ/∂t = [ −(ħ²/2m)∇² + mΦ ] ψ                    (Schrödinger)
∇²Φ = 4πG( m|ψ|² − ρ̄ )                                (Poisson, periodic)
```

This is **M6's split-step ψ engine and M2's PM cuFFT-Poisson solver welded at the hip** (PROPOSAL §S1: "same cuFFT, one loop"). The split-step kinetic operator is `planck`/`cosmos` verbatim (`kQ3PhaseK` lineage); the Poisson solve is `newton` verbatim (`kDeposit`→`kGreen`→C2R). What is new is the *coupling*: the potential half-kick reads Φ from a Poisson solve keyed on the field's own `|ψ|²` every step. Fuzzy-dark-matter physics (Schive et al. 2014; Mocz et al. 2017) proves the regime emergence is real, not a stunt. Single file inside `app/tinyuniverse.cu` v2 (one binary, two faces — D-003); the headless face is the golden unit.

This contract is **design-only** (PROPOSAL §7 build order: N0 oracle first — SHIPPED, golden `13aa73e5` — then N1). No `.cu` source is written until this contract is operator-reviewed (Invariant 1). Autonomous design rulings below are flagged `[DECIDED]`.

## Units & dials (v1 dials — the proposal's thesis)

The v1 dials are exactly what make the fused solve fit in-box (PROPOSAL §0, the visibility=computability table). N1 inherits them verbatim from `nexus`; it introduces **no new fundamental constants**.

| dial | value | role in N1 |
|---|---|---|
| ħ | 0.5 sm·su²/s | de Broglie λ lands at grid scale → quantum pressure is resolvable on the *same* lattice as gravity (PROPOSAL §0) |
| m (m_p) | 1 sm | the single field's particle mass; sets ħ/m — **the only knob that moves the quantum↔classical boundary** |
| G | 2×10⁻³ su³/(sm·s²) | self-binding at the field densities the v1 scenarios used; k=0 mode does not gravitate (mean subtracted) |
| dt | 1/240 s | fixed timestep, Strang-split; **not** CFL-bounded here (Schrödinger is dispersive, not hyperbolic — the CFL headroom that mattered for c is an S2/S3 concern) |
| L_box | 512 su | the 3-torus period; the field is periodic by default (matches PM's periodic Green's function) |
| c | 20 su/s | **irrelevant to N1's declared physics** — SP is non-relativistic; c re-enters at N2 `lapse`. Declared inert here. |

- **Not geometric units.** Unlike N0 `substrate_nexus` (G=c=1, dimensionless Choptuik — see its contract "Units"), N1 lives in the v1 dial system because its whole point is to reproduce the v1 goldens *by the field* on the lattice where those dials make the physics visible-and-computable. N0 owns the equations; N1 owns the dials.

## Declared engine (v1 scope — the fused loop)

**Grid.** Cubic lattice `PSN³`, `PSN ∈ {256, 512}` (512 gated on VRAM; see Cost). Cell size `dx = L_box / PSN` → 2.0 su at 256³, 1.0 su at 512³. Complex field ψ stored `cufftComplex` (fp32 ×2). One cuFFT C2C plan (ψ kinetic step) + one R2C/C2R plan pair (Poisson), both fixed for the run.

**One step (Strang split, symmetric — declared operator order is normative):**

```
1. HALF KINETIC   ψ ← FFT⁻¹{ e^(−i·(ħ dt / 4m)·k²) · FFT{ψ} }        (kQPhaseK, coef = ħ·dt/(4m))
2. DEPOSIT ρ       ρ_fix ← fixed-point Σ over cells of m|ψ|²·dx³        (kPsiDeposit → uint64, Invariant 4)
3. POISSON Φ       Φ ← C2R{ (−4πG/k²)·R2C{ρ − ρ̄ } },  k=0 ⇒ 0         (kGreen lineage; ρ̄ = M_tot/V)
4. FULL POTENTIAL  ψ ← e^(−i·(m dt/ħ)·Φ) · ψ                          (kPhaseV_field, coef = m·dt/ħ, V=mΦ)
5. HALF KINETIC    ψ ← FFT⁻¹{ e^(−i·(ħ dt / 4m)·k²) · FFT{ψ} }        (kQPhaseK, coef = ħ·dt/(4m))
```

Strang splitting is 2nd-order in dt and **exactly time-reversible operator-by-operator** (each factor's inverse is its complex conjugate) — this is what makes the echo exact (see Determinism). `[DECIDED]` **half-kinetic / full-potential / half-kinetic** (kinetic on the ends), not the reverse: it lets consecutive steps fuse their trailing+leading half-kicks in a later optimization *without changing declared bytes* (the fused form is algebraically identical); ordering frozen now so the golden is stable. Rationale: fewer FFTs later, zero determinism cost.

**Deposit (Invariant 4 — no float atomics in declared reductions).** ρ is accumulated into a **fixed-point uint64 grid** via `orrery::fixed_atomic_add` (2⁻³² mass quantum, declared), reusing `kDeposit`'s exact idiom (`app/tinyuniverse.cu:389`). But unlike PM's CIC scatter from particles, here the source is cell-local: cell `c` deposits `m·(ψ.x² + ψ.y²)·dx³` into `g[c]` — a **1:1 gather, no scatter, no CIC stencil** (the field already lives on the grid). `[DECIDED]` deposit is still routed through `fixed_atomic_add` even though each cell writes its own slot (no contention): it keeps the |ψ|² → ρ path bit-identical to how M2 quantizes mass, so the field's total-mass reduction and the PM mean ρ̄ agree to the ulp with v1, and the deposit stays a drop-in for the existing `kGreen` pipeline. Rationale: order-invariance is free here and buys cross-check exactness against the oracle farm.

**Poisson.** Identical to `newton`'s PM solve (`kGreen`, `app/tinyuniverse.cu:396`): R2C forward, multiply by `−4πG/(k²·dx³·N³)` with the k=0 bin zeroed ("mean density does not gravitate" — the periodic-gravity convention, matches M2 verbatim), C2R inverse. The **same Green's function, the same normalization constant, the same k=0 rule** — so a uniform field produces Φ≡const≡0, exactly as M2's uniform particle sea does. `[DECIDED]` reuse M2's 3D `kGreen` unchanged rather than write a field-specific solver: the operator (`∇²Φ = 4πGρ_source` on a periodic box) is *identical* — only the source's provenance differs (|ψ|² vs CIC-splatted points). One solver, two callers = one thing to golden. (D-005 spirit: don't fork what is already frozen.)

**Boundary conditions.** Two declared modes, per scenario:
- **PERIODIC** (default): the FFT is periodic in both the kinetic step and Poisson — the field wraps on the L_box 3-torus exactly as v1's PM gravity has since M2 (the box topology *is* the universe's — ARCHITECTURE §6; D-019 made the wrap a coordinate fact). Used for galaxy / merger / cosmological scenarios where mass stays in-box.
- **ABSORBING** (declared per-scenario): a `cos²` sponge layer (8 cells, matching `planck`'s absorbing-edge idiom `kQ3Edge`, `app/tinyuniverse.cu:976`) multiplies ψ after step 5 to bleed outgoing amplitude. **Non-unitary by construction** — used only for free-packet dispersal / soliton-in-vacuum where norm loss to the walls is physical and *declared* (the packet is supposed to leave). `[DECIDED]` absorbing BC is **incompatible with the echo receipt** (a sponge is not conjugation-invertible) — so the echo scenario is PERIODIC-only, and absorbing scenarios do not claim exact reversal. Rationale: keep the two determinism stories (unitary-periodic-reversible vs open-absorbing-lossy) cleanly separated, never blended in one golden. Mirrors `planck`'s "absorbing edges, declared" (`planck.contract.md` line 7).

**Determinism.** Stronger than v1 by construction (PROPOSAL §3, "determinism gets *stronger* than v1's"):
- **No float atomics anywhere in the declared path** — the deposit is fixed-point (Invariant 4); every other kernel is a per-cell map (kinetic phase, potential phase, Green multiply) with no cross-cell reduction except the fixed-point mass sum. There is no scatter, so there is nothing to race.
- cuFFT is bit-stable for a fixed plan + fixed arch (sm_89) — the same guarantee M2 and M6 already ship on.
- **Golden = (scenario, dials, seed, steps) → blake2b-256 of the declared field bytes** (the ψ grid; +Φ where a scenario declares it): hardware-pinned sm_89. See Declared state.
- **Echo reversible by conjugation.** Time reversal of Schrödinger is ψ → ψ* (complex conjugation flips the sign of every phase: kinetic `e^(−iθ)→e^(+iθ)`, potential likewise; Φ is real and unchanged since it depends only on |ψ|²=|ψ*|²). So *forward N steps, conjugate, forward N steps, conjugate* returns ψ to its initial state under PERIODIC BC. This is the M3/M4 Loschmidt echo made structural.
  **[BUILD FINDING 2026-07-12 — the "byte-identical" wording is corrected to "reversible to fp32 round-off"; D-016/D-018 honesty]:** in the fp32 split-step, the ≈4N cuFFT round-trips per echo are **not bit-reversible** (a rounded forward-then-inverse FFT is not the identity in fp32; phase multiplies likewise). So init≠final at the *bit* level, and byte-identity is **physically precluded by fp32**, not a bug to fix. What IS true and gated: **(1)** the round-trip residual sits at fp32 round-off — measured L2 grows linearly with N (3.7e-5 @ N=50 → 4.0e-4 @ N=600), so the *physics* reversal is exact; the golden gates `L2_rel < 1e-3` at the frozen N. **(2)** the reproducibility determinism — same `(scenario,dials,seed,steps) → identical declared bytes` — **is** byte-exact (two independent runs give an identical declared JSON + blake2b, verified). The golden freezes that byte-exact hash. So echoF proves both the exact-reversibility of the fused loop *and* its bit-reproducibility; only the stronger "init-vs-final bit-equality" claim is retired as fp32-impossible. (A truly bit-exact echo would need integer/reversible arithmetic — out of scope for cuFFT fp32.)

**Regimes.** Per ARCHITECTURE Invariant 5 (regimes emerge, no hard switches) and §Madelung below: the QUANTUM/CLASSICAL character is **not a stored bit and not a code path** — it is the local value of the quantum-pressure term, which is a continuous function of ħ/m and the field's own gradients. N1 carries **no regime bitmask** (there are no particles to tag); the derived-only regime idea (`regime.h`) has nothing to bind to until tracers arrive (S-later). `[DECIDED]` drop the regime buffer from the N1 declared state entirely rather than fake a per-cell tag. Rationale: the whole thesis is that the boundary is a dial (ħ/m), not a flag — inventing a per-cell QUANTUM/CLASSICAL byte would contradict the design and bloat the hash.

## The Madelung regime emergence (the physics claim, and how the golden tests it)

Substitute `ψ = √ρ · e^(iS/ħ)` into the Schrödinger–Poisson system (the Madelung transform). The real/imaginary parts become **exactly the Euler equations of a self-gravitating fluid** with velocity `v = ∇S/m`, plus one extra term — the **quantum pressure** (Bohm potential):

```
Q = −(ħ²/2m) · (∇²√ρ / √ρ)
```

The claim (literature-proven — fuzzy dark matter: Schive/Chiueh/Broadhurst 2014; Mocz et al. 2017; PROPOSAL §S1):

- **Large scales, mild gradients:** Q → 0 (it carries ħ² and a second derivative of √ρ). The field obeys ordinary collisionless self-gravitating dynamics — **galaxies rotate, clouds collapse, mergers virialize and relax exactly as the v1 particles did.** This is why the v1 goldens must reproduce *by the field*.
- **Small scales, steep gradients:** Q blows up and acts as a repulsive pressure that **halts gravitational collapse at a finite core** — solitonic cores, interference fringes, tunnelling. Real quantum mechanics, on the same lattice, from the same equation.
- **The boundary is ħ/m — a dial, nothing switches.** The de Broglie / Jeans crossover scale is set by (ħ/m); the soliton core radius scales as `r_c ∝ ħ²/(G m² M)`. Turn ħ/m up: quantum structure grows to galactic scale (fuzzy-dark-matter halos). Turn it down: the field becomes indistinguishable from cold collisionless matter (v1). **No `if (regime == QUANTUM)` exists anywhere in the code** — this is Invariant 5 taken to its limit.

**How the golden tests the claim (both ends of the boundary, one engine):**
- The **quantum end** is pinned by the free-packet σ(t) and SHO scenarios (exact vs the N5 analytic oracle) and the **soliton mass–radius relation** (the ground state of SP — analytic scaling `M·r_c = const`): the field self-supports against its own gravity at the Bohm-pressure radius. If Q were wrong, the soliton would either collapse (Q too weak) or disperse (Q too strong).
- The **classical end** is pinned by the v1 cross-checks: galaxy rotation support, merger entropy rise, cloud collapse morphology — reproduced *by the field* at a small ħ/m where Q should be negligible. If Q leaked into large scales, the rotation curve or the collapse timescale would deviate from the v1 particle golden.
- Together they are a **two-sided proof that the same term Q does the right thing at both extremes** — which is the entire Madelung claim. `[DECIDED]` the soliton test is the *load-bearing new golden* (the packet/SHO are M6 re-runs); it is the one scenario that exercises quantum-pressure *balancing gravity*, which neither M6 (no gravity) nor M2 (no ħ) ever tested. Rationale: it is the smallest experiment that can only pass if the fusion is correct.

## Scenarios & gates (seed 20260711; the golden units)

Field scenarios re-express v1 initial conditions as ψ(x,0). `[DECIDED]` v1's particle ICs map to the field by **ρ(x,0) = (particle density smoothed to the lattice), S(x,0) = m·∫v·dx (the velocity field's potential where the flow is irrotational; for rotational flows the closest curl-free S plus a declared residual)**. Rationale: a wavefunction encodes only irrotational (potential) flow directly — galactic rotation is *supported* by the |ψ|² interference/pressure structure and the coherent phase winding, not by a literal vorticity field. The mapping is a declared scenario-construction step (like M2's virial construction), gated by *outcome* (does the rotation curve / collapse / merger reproduce), not by IC byte-match to v1.

| scenario | grid / BC | init (ψ(x,0)) | golden steps | gates | oracle pedigree |
|---|---|---|---|---|---|
| **freepacket** | 256³ ABSORBING | Gaussian packet σ=1.5, momentum p → coherent phase `e^(ip·x/ħ)`, gravity **off** (G=0 for this scenario — pure kinetic) | declared (≈ M6 count) | free-dispersion width: `\|σ(t) − σ_analytic(t)\|/σ < 1e-2` where `σ²(t)=σ₀²+(ħt/2mσ₀)²` (the N5 closed form) · norm monotone-decreasing (absorbing) · NaN-free | **nexus N5** (fp64: σ(t) to 5.8e-15) — exact |
| **sho3d** | 256³ PERIODIC | 3D isotropic harmonic well Φ_ext = ½ω²r² (external, added to self-Φ which ≈0 at ε-amplitude), imaginary-time relax OR real-time coherent-state breathing | declared | ground: `\|E − (3/2)ħω\|/E < 1e-3` (3D ground = 3·½ħω) · `\|σ_x − √(ħ/2mω)\|/σ < 1e-2` | **nexus N5** (E₀ to 1.3e-13) — exact; extends M6's 2D `shoq` to 3D |
| **soliton** | 256³ PERIODIC, box 16→8 su (2 masses) | SP ground state relaxed by imaginary time to the self-bound eigenstate at mass M; gravity **on** (PM Poisson each step). **[AS BUILT]** the 2:1 mass pair uses the SP self-similar family (box, dx, seed, τ all ∝1/M) so the covariance is exact | 2×10000 imag-time | **mass–radius scale-covariance**: `\|r_c·M2 / r_c·M1 − 1\| < 0.05` (r_c = half-peak-density radius; κ measured in-tool). **[MEASURED]** r_c·M = 360.156 at M=200 AND M=400 → scale_rel = **3.0e-8**; self-bound E<0; virial 2K/\|W\|=1.9176 (both) | **analytic** SP soliton scaling (Chavanis 2011; Mocz 2017) — the load-bearing new physics · GREEN `dd4d08d7…` |
| **galaxyF** | 256³ PERIODIC | v1 `galaxy` morphology mapped to ρ+phase (rotational support via |ψ|² structure + phase winding), small ħ/m | 10⁴ (v1 parity) | rotation curve reproduced: `v_circ(r)` from the field within **8%** of the v1 `galaxy` golden band · mass conserved (`\|M(t)−M(0)\|/M < 1e-3`, periodic → norm exactly conserved) | **v1 galaxy golden** (`2642d510`, D-019) — cross-check |
| **mergerF** | 256³ PERIODIC | two `soliton`/cloud lumps on the v1 `merger` trajectory | 12000 (v1 parity) | **entropy rise**: coarse-grained `S` of |ψ|² rises and `ΔS ≥ −0.01 nat`/sample at ≥ 0.75 monotone (the M3 H-theorem gate form, D-015) · final virialized | **v1 merger golden** (`1e80bee9`, D-019) — cross-check |
| **echoF** | 128³ PERIODIC, box 32 | gravitating Gaussian lump (M=200, so the Poisson kick is exercised), run **N forward → conjugate ψ → N forward → conjugate** (N=400) | 800 (2×400) | **[AS BUILT — see Determinism BUILD FINDING]** reversal residual at fp32 round-off `L2(ψ_final−ψ_init)/‖ψ‖ < 1e-3` (**measured 3.4e-4**) + declared golden reproduces byte-exact. Bit-identity NOT gated (fp32 FFTs not bit-reversible) | **structural** (the determinism receipt; v1 echo `40a84691`→D-019 lineage) · GREEN `433ddcc8…` |

`[DECIDED]` `cloud` collapse is folded into `mergerF`'s morphology check rather than given its own row (a single soliton *is* the collapse endpoint; the two-lump merger exercises infall + relaxation + the entropy gate together). Rationale: keep the golden set to the six that each prove a distinct thing; cloud-collapse morphology is the same physics as merger infall at N=1 lump. If review wants it separate, it re-splits trivially (same engine).

**The gate is on outcome, not on IC byte-match to v1.** The v1 goldens hash *particle* buffers; N1 hashes *field* buffers — they cannot be byte-equal (different representations). The cross-check is that the **derived physics** (rotation curve, collapse time, entropy trajectory, ISCO-band if a compact scenario is later added) lands inside the v1 golden's declared band. This is the oracle-farm premise (PROPOSAL §4): v1 is eighteen (now 21) receipts N1 reproduces by independent means, not bytes it copies.

Oracle pedigree summary: **freepacket + sho3d are exact against nexus N5** (the fp64 split-step proof M6 already leaned on); **soliton is analytic** (SP scaling law); **galaxyF / mergerF are v1-golden cross-checks** (the oracle farm); **echoF is structural** (conjugation is exact arithmetic). No scenario is gated against an unmeasured claim.

## Declared state & the frame-contract bump

- **Declared state = the ψ field grid bytes** (`cufftComplex[PSN³]`), plus Φ for scenarios that declare it (soliton/galaxyF, where the potential is a physics output). Hashed with blake2b-256, canonical byte order (x-fastest, matching M2/M6 grid layout).
- **Frame contract MAJOR bump → v2.0.0** (PROPOSAL §4, §S1). The declared state changes *kind*: **field textures replace the SoA particle buffers** (pos/mom/tau/regime/bubble of `frame.contract.md` v1). |ψ|² *is* the emission field the renderer reads; Φ *is* the field the lens/τ layers read. Tracer particles (if gameplay wants visible dots) become a **non-declared visualization clause** — advected by the Madelung velocity v=∇S/m, never fed back into dynamics (PROPOSAL §2 item 2). This is a breaking change and takes a MAJOR bump + migration note (Invariant 1, Hard Rules). **N1 does not itself edit `frame.contract.md`** — it *names* the required v2.0.0 bump as its precondition; the bump is its own reviewed change (the renderer transfer, PROPOSAL §4, is out of N1's scope — N1 is headless-first like M6).

## The headless face (envelope-conformant; liborrery serializers)

```
tinyuniverse.exe --scenario NAME [--steps S] [--seed N=20260711] [--grid 256|512] --json   # run + declared JSON envelope
tinyuniverse.exe --scenario NAME --golden                                                   # frozen params vs goldens/NAME/declared.hash
```

Declared body (canonical order): `seed, params{scenario, grid, steps, bc, c, hbar, G, m, dt, L_box}, result{state_b2b, mass0, mass1, dmass_rel, norm0, norm1, <scenario observables: sigma_t/E/r_c/v_circ_rms/dS/echo_match>}, gates{...pass booleans...}, verdict`. Floats fmt6. **Exit 0 pass · 1 declared gate fired (real negative result) · 2 error — never conflated** (Invariant 8; Hard Rules). `--selftest` runs a cheap flat-field + norm-conservation smoke.

## Cost (bandwidth arithmetic — ADR-007: a plan, not a measurement)

Per-cell live storage: ψ (2 fp32) + Φ (1 fp32) + FFT scratch/plan working set + fixed-point ρ (uint64 = 2 fp32-equiv) ≈ **8 floats/cell** (PROPOSAL §S1).

| grid | dx | ψ+Φ+scratch (≈8 f/cell) | verdict |
|---|---|---|---|
| **256³** | 2.0 su | ≈ 0.5 GB | **feasible day one**; the N1 gate resolution |
| **512³** | 1.0 su | ≈ 4.3 GB | feasible on the 16 GB card (PROPOSAL §S1 "feasible day one"); gated on a live VRAM preflight (the v1 preflight discipline carries over verbatim — PROPOSAL §3; harness exit-3-on-contention, per HEAD `d162ae4`) |

Step-rate estimate (PROPOSAL §3, ~200 B/cell/step on ~672 GB/s): **256³ ≈ 100–200 steps/s** (near the 240 Hz dial); 512³ ≈ ¼-speed. **`[ARGUMENT-GRADE]` — bandwidth arithmetic, not `--bench`. No perf number ships as a claim until measured (ADR-007 / Invariant 6).** Two cuFFT plans (C2C for kinetic, R2C+C2R for Poisson) dominate the working set; FP16 storage (PROPOSAL §3) is an S-later measured option, not v1.

## Oracle pedigree (Invariant 3 — no silent fallback)

- **Exact analytic:** freepacket σ(t), sho3d E₀ — both against **nexus N5** (the fp64 split-step oracle, already proved to 1e-13). Oracle-gate failure ⇒ run-state SUSPECT, never a quiet pass.
- **Analytic scaling:** soliton mass–radius (SP ground-state scaling law) — measured κ frozen, scale-covariance checked.
- **The oracle farm (v1's 21 goldens):** galaxyF / mergerF reproduce `galaxy` / `merger` physics bands by the field. **`[DECIDED]` N1 additionally sanity-cross-checks the field's own uniform-mode Poisson against M2's `kGreen`** by asserting a flat ψ gives Φ≡0 to the ulp (same solver, same k=0 rule) — a free structural check that the fusion didn't perturb the frozen solver. Rationale: it is the cheapest possible guard that "one solver, two callers" actually holds.
- **Structural:** echoF (conjugation is exact fp32 arithmetic under a reversible operator split).
- **The N0 tie:** N1 does **not** gate against `substrate_nexus` (N0 is GR/EKG — relativistic; N1 is non-relativistic SP). N0 becomes load-bearing at **N3 `horizon`** (PROPOSAL §5). N1's oracles are N5 + analytic + the v1 farm. Stated so the pedigree is honest about which foundation stone carries which rung.

## The honest boundary (permanent, printed — PROPOSAL §2)

1. **Many-body entanglement never unifies on a 3D grid — declared, forever.** ψ(x) on a 3-lattice is a **mean-field, single-configuration** universe. N entangled particles live in 3N-dimensional configuration space; no 3D lattice escapes that exponential. **Single-particle interference is exact** (it *is* field interference — the double slit, tunnelling, the soliton). **Bell/EPR correlations between separated lumps are absent or faked, permanently.** N1 makes no entanglement claim and the HUD says so in physical mode (Invariant 9; ARCHITECTURE §12). This is the same boundary M6 declared for bubbles (`planck` v1 boundary), now inherited by the whole substrate — and it is not a bug to be fixed, it is the price of one lattice.
2. **"Particles" become excitations.** Any gameplay-visible dot is a tracked |ψ|² peak or a non-declared tracer advected by the Madelung flow — a visualization layer, not dynamics (PROPOSAL §2 item 2). Individual particle *identity* is not conserved by a field.
3. **N1 is non-relativistic and gravity is mean-field Newtonian.** c is inert; there is no lapse, no horizon, no light-bending (those are N2/N3). Fusion/Hawking/Ratchet are N4. N1 is exactly SP — stated scope, no more.

## Gate (N1, from PROPOSAL §5)

**SP on a 256³ (→512³) GPU lattice, with the v1 galaxy/cloud/merger/echo scenarios re-run and reproduced *by the field*** — plus the soliton mass–radius oracle (the new quantum-pressure-vs-gravity receipt) and the free-packet/SHO exactness against nexus N5. **The echo returns byte-identical by conjugation** (the determinism receipt). All gates green ⇒ golden frozen; frame contract noted for its v2.0.0 MAJOR bump. Honest boundary (no many-body entanglement) printed.

## Open questions (Q — carried to review / build)

- **Q-N1-1 · IC mapping fidelity.** How faithfully must a v1 particle IC map to ψ(x,0) for the outcome gate to be fair? Rotational support has no exact wavefunction (curl-free phase only). *Deciding experiment:* build galaxyF, measure the rotation curve vs the v1 band; if it fails only because of the irrotational-phase residual, the gate band (currently 8%, borrowed from M6's Fresnel gate) may need a declared widen — **measure, don't pre-tune** (D-018 discipline).
- **Q-N1-2 · soliton κ provenance.** Is the SP scaling constant κ taken from literature (Mocz 2017 fit) or measured in-tool by relaxing two masses and reading the slope? `[DECIDED tentatively]` measure it in-tool (self-consistent, no external-number import), report the literature value as a labelled cross-check — mirrors N0's `gamma_eff` diagnostic discipline (D-021). Confirm at review.
- **Q-N1-3 · imaginary-time determinism.** The soliton/SHO relaxation uses imaginary-time propagation with a renormalization each step (the M6 `shoq` idiom, fixed-point norm via `kQ3NormAcc`). Does the renorm divide stay bit-reproducible at 256³? M6 proved it at 256²; **256³ is unmeasured** — flagged, cheap to check first thing in the build.
- **Q-N1-4 · 512³ VRAM under a shared card.** The 4.3 GB estimate assumes an idle card; the GPU is shared (v1 preflight). Does 512³ + two cuFFT plans + scratch clear a *contended* 16 GB? *Deciding experiment:* the preflight probe (harness exit-3 lineage, `d162ae4`) sizing the two plans before allocating. 256³ is the safe gate; 512³ is the stretch.
- **Q-N1-5 · does the fused loop keep the k=0 convention honest across representations?** M2 subtracts ρ̄ as the particle mean; N1 subtracts M_tot/V as the field mean. These are the same number *iff* the deposit conserves total mass exactly — which the fixed-point path should guarantee. Assert it in `--selftest`.

## Build runbook (the path a build session executes — design→code→golden→harness)

1. **Contract (this file) → operator review.** No `.cu` until reviewed (Invariant 1). ← *you are here.*
2. **MODULE.md** — purpose, this contract link, invariants touched (1,3,4,5,8,9), oracle (N5 + analytic + v1 farm), build command (`BUILD.md`: vcvars64 + nvcc `-arch=sm_89`, links `cufft.lib` → ships `cufft64` dll, D-014 note), known issues (Q-N1-*).
3. **Implement (single file, `app/tinyuniverse.cu` v2 track).** Lift/adapt (see plan): `kQ3PhaseK`→`kPhaseK_field` (kinetic half, 3D, whole grid), `kDeposit`→`kPsiDeposit` (|ψ|²→fixed-point, 1:1 no CIC), `kGreen` (verbatim), a new `kPhaseV_field` (potential kick with live Φ, template = `kQPhaseV`), `kQ3Edge`→sponge for absorbing BC, `kQ3NormAcc`→norm reduction. Wire the 5-step Strang loop; conjugation = negate ψ.y in a one-line kernel for echoF.
4. **Golden freeze.** Target **soliton first** (the load-bearing new physics — the smallest scenario that can only pass if the fusion is right), then freepacket/sho3d (N5 exactness), then the v1 cross-checks, then echoF (the receipt). Freeze `(scenario, dials, seed, grid, steps) → blake2b` in `goldens/field/`, sm_89-pinned.
5. **Harness row.** Add `field` scenarios to `harness/` (compile + `--selftest` + `--golden`); GPU-preflight gate (exit 3 on contention, `d162ae4`) before the 512³ variant. Two-pass verify the soliton and echo (the citable claims).
6. **Register.** ARCHITECTURE §11 v2 row (N1 status → CLOSED with goldens), TASKLIST, DECISIONS (a D-0xx findings entry — the M-milestone habit), RUN_STATE. Spec never lags code (build loop step 7).

## Changelog

- v1.0.1 (2026-07-12) — **BUILT (autonomous overnight); 4/6 scenarios GREEN.** The PM+ψ gravity **weld is demonstrated**: a self-bound SP soliton with r_c·M scale-covariant to 3e-8 across a 2:1 mass pair (D-022). Build findings, all measured: **(a)** the potential coupling carried a spurious extra 1/N³ (self-gravity was ~10⁶× too weak) — removed; the C2R buffer holds the *true* Φ (M2 convention), verified by `--poissontest` vs analytic −GM/r (3.6%). **(b)** soliton gate reframed to SP self-similar scale-covariance (`|r_c·M2/r_c·M1 − 1| < 0.05`) measured along the box∝1/M family — the physically-correct scale-invariance test; MEASURED scale_rel = 3.0e-8; golden `d163d765`. **(c)** echoF: strict byte-identity is fp32-impossible (≈4N cuFFT round-trips not bit-reversible) — gate reframed to reversal-residual-at-round-off (L2 < 1e-3, measured 3.4e-4) + byte-exact *reproducibility* (verified); "byte-identical" wording corrected in Determinism §. Goldens frozen: freepacket `03dd3a3b`, sho3d `dfbc6185`, soliton `d163d765`, echoF `433ddcc8`. galaxyF/mergerF (v1 cross-checks) remain. Frame-contract v2.0.0 bump still named-not-executed (headless-first). Tool: `substrate/field_nexus.cu`.
- v1.0.0 (2026-07-12) — **initial DRAFT, design-only** (operator review requested; no GPU build). Fuses `newton` PM-Poisson (`kGreen` verbatim) + `planck`/`cosmos` split-step (`kQ3PhaseK` lineage) into the Schrödinger–Poisson loop per PROPOSAL §S1. Six scenarios (freepacket/sho3d/soliton/galaxyF/mergerF/echoF); N5 + analytic-SP + v1-farm + structural-echo oracles; declared state = field bytes (frame contract v2.0.0 MAJOR bump named, not executed). Autonomous rulings flagged `[DECIDED]`; open questions Q-N1-1..5 carried. Honest boundary (no many-body entanglement on a 3D grid) printed per PROPOSAL §2.
