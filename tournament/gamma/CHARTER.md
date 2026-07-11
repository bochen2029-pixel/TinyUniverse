# TOURNAMENT — how N0 resolves the Choptuik exponent (or whether the crown should be reframed)

**Owner:** Bo Chen · **Run by:** Claude (Fable 5), human-in-the-loop at each phase boundary · **Isolation:** git worktree `substrate-gamma-tournament` off TINY UNIVERSE `master` @ `8060bd0`. Nothing here touches `master` until the operator reviews it.

## The question (neutral)

TINY UNIVERSE v2 N0 (`substrate/substrate_nexus.cpp`) demonstrated the Choptuik **Type-II transition** but could **not** resolve the universal critical exponent γ on a uniform polar-areal grid (D-021: the self-similar curvature caps at the grid; refining goes chaotic — p\* wandered 0.40→0.356 from N=800→1600). **How should N0 resolve γ ≈ 0.374 — or should its crown be reframed to a deepest-GR receipt a CPU-modest grid can actually nail, deferring γ to the GPU N3 `horizon` stage?**

## Why this fork first — it validates the machine

γ ≈ 0.374 is **ground truth** (Choptuik 1993; Gundlach γ=0.3737, echo period Δ≈3.44). A tournament + ORRERY loop that takes my failed uniform grid and lands 0.374 (or cleanly reframes the crown) **proves the methodology on a problem with a known answer** before we trust it on the genuinely-open forks (the N1–N4 substrate architecture). Oracle before the thing it gates — the project's own law, applied to the method.

## The incumbent (the champion to beat) and the graveyard

- **Incumbent:** `substrate_nexus.cpp` v1.0.0 — spherical EKG, polar-areal, uniform grid, constrained log-metric evolution, RK4. Golden `13aa73e5`. It ships the transition honestly; its crown (precise γ) is deferred.
- **Graveyard (seeded):** *uniform-ultrafine polar-areal.* Corpse reason (D-021, measured): a fixed grid caps the self-similar central curvature (≈(φ/dr)²); refining does not converge, it makes near-critical chaotic. Re-proposal requires a **genuinely new argument**, not "just use more points."

## The fixed constraints (the box every approach must live in)

Any winning approach must obey TINY UNIVERSE / ORRERY doctrine, or it does not ship:
1. **Determinism-or-it-doesn't-ship** — (params, seed) → byte-identical declared output; a frozen `(params)→blake2b` golden.
2. **Golden-gated, honest** — the RAYFORMER / D-016 / D-021 rule: a reported number is backed by an anchor (ground truth), metamorphic relations (invariance), and redundant recovery — or it is declared, never faked.
3. **CPU, ≤ ~5 min** for the golden run (N0 is the GPU-free foundation stone). Single-file-where-possible.
4. **Exit codes** 0 pass · 1 declared gate fired · 2 error — never conflated.

## The three refuter lenses (phase 2 worldviews; each states what counts as a KILL)

- **DETERMINISM** — every mechanism must be buildable as a deterministic, golden-freezable, single-file CPU tool inside the time budget. A step that needs a nondeterministic solver, an unbounded runtime, or a research-grade AMR library that doesn't fit the doctrine is a kill of the claim that leans on it.
- **RESOLVABILITY** — can the approach's grid **actually see** the self-similar exponent, *measured not argued*? A pre-registered pass-bar (γ within tol at declared R²) that no buildable instance can reach is a kill. **This is where you drive ORRERY** (see arming) to run the critical-point search and show the number, rather than assert it.
- **ORACLE-HONESTY** — is the claimed γ anchored (matches 0.374), metamorphic-stable (invariant under grid/param changes it should be), and redundantly recovered (two independent observables agree — e.g. mass-scaling and subcritical-curvature scaling)? A number that is really grid-noise dressed as a result is a kill (the D-018/D-021 discipline).

## Adjudication rules (phase 3)

- **Rule every conflict with a reason; never average mechanisms.**
- **Respect the graveyard.** Every deletion ships a **pre-registered reinstatement trigger** (the reversibility lemma as house discipline).
- **Steal list:** a good idea from a losing approach is harvested regardless of the approach's fate.
- **Internal convergence is weak evidence.** Every verdict in this tournament is argument-grade until its **pre-registered deciding experiment runs** (via ORRERY / a built tool). The register holds the doubt.

## ORRERY arming — every subagent gets this block

You are armed with **ORRERY** (`C:\ORRERY`), a headless, **deterministic** simulation instrument you *call* to make quantitative claims **evidence-grade, not argument-grade**. **Read-only use — run the tools, never edit ORRERY source.**

**Universal envelope (every tool):** `<tool> [--param V] --seed N [--json|--csv PATH] [--selftest] [--golden]` · exit `0` pass / `1` a declared gate fired (a real negative result) / `2` error · `--json` → one JSON envelope on stdout · determinism: (params,seed) → byte-identical declared output.

**How to call (verified working):**
- Python tools directly: `python C:\ORRERY\tools\autotune\autotune.py <flags> --json` (run with CWD `C:\ORRERY` if a tool reads `goldens/`).
- The catalogue via the MCP surface (PowerShell here-string dodges JSON-quote escaping):
  ```
  Set-Location C:\ORRERY
  $req = @'
  {"jsonrpc":"2.0","id":1,"method":"tools/call","params":{"name":"list_tools","arguments":{}}}
  '@
  python tools\mcp\mcp.py --once $req
  ```
  MCP tools: `list_tools`, `describe_contract {tool}`, `run_tool {tool,params}`, `sweep {tool,sweep,metric,lo,hi,target,...}` (drives autotune), `golden_status`. Every run/sweep response embeds the declared **blake2b** — **cite it**.

**The catalogue (9 tools; Python ones run now with no GPU, CUDA ones need the shared card):**
| tool | lang | what it does — when to reach for it |
|---|---|---|
| **autotune** | py ✅ | sweep one param, locate an **argmax (band)** or **level-crossing (critical point)** against a **pre-registered `--target`** (G-OFF-TARGET fires if elsewhere). Can subprocess a *real* tool as the objective and read a JSON metric. **It found ratchet's ρ_c = 0.2581 vs analytic 0.25 blind — this is the tool for a p\*/γ search.** |
| **mcts** | exe | generic root-parallel search over a supplied action/param space |
| **ratchet** | exe | the inscription threshold `(1−p)ρ=p` at scale (the substrate's Ratchet lattice, not Choptuik) |
| **posit** | py ✅ | parsimony auditor — posit-budget count (an Occam check on a design) |
| **algebra** | exe | crossed-product entropy (Type II/III) |
| **someone / mcts / orreryd / mcp** | — | self-model criterion / search / campaign daemon / this surface |

**The rule:** every quantitative claim you make is either backed by an ORRERY run (cite the tool, params, and declared blake2b) or is flagged `[ARGUMENT-GRADE]`. If the deciding experiment for your claim needs a tool that does not exist yet (e.g. a double-null EKG gear), say so precisely and specify the tool's contract — do not fabricate a number.

## Phases & deliverables

- **phase0** — `assumption_ledger.md`: the forks I settled solo in N0, each with its argument, its pre-registered deciding experiment, and an assigned attack stance. *(authored by the lead; the operator reviews before phase 1.)*
- **phase1** — one doc per approach-persona (double-null, CSS/self-similar, AMR, spectral, reframe-the-crown, …): a rigorous approach with a reversibility lemma, a pre-registered falsifier, a costed deciding experiment, accurate cited literature, and whatever cheap ORRERY experiments it can run now.
- **phase2** — each approach × the three lenses: KILLED / WOUNDED / SURVIVED claims (mechanism-precise), a STEAL LIST, a weighted scorecard — grounded in ORRERY runs where a claim is quantitative.
- **phase3** — `hybrid.md`: rule the forks with reasons; the graveyard; the chosen approach + its pre-registered deciding experiment.
- **phase4** — judges score {resolves γ to tol? · obeys doctrine? · generalizes to N1–N4?}; then **build the winner and measure γ for real** — the validation.
