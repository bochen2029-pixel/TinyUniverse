# reframe-crown — ORRERY experiment ledger

Read-only use of ORRERY (`C:\ORRERY`), CWD `C:\ORRERY`. Two Python tools, run now, no GPU.
Every declared envelope re-run 2× byte-identical (determinism confirmed).

## Tool golden baselines (cited for provenance)

| tool | command | declared blake2b (golden) | exit |
|---|---|---|---|
| posit v1.0.0 | `python tools/posit/posit.py --golden --json` | `7a22dd229a42ce46a6c102f0545f83022b975dc39d5f1794cd6019e6f5a20e44` | 0 |
| autotune v1.0.0 | `python tools/autotune/autotune.py --golden --json` | `c79002f23cf236baab5ecdb5753603a7a3853199f750b0139771d4e6cdd55bbe` | 0 |

## E1 — parsimony audit: reframe vs chase-γ (posit)

- **case:** `posit_case_crown.json` (pinned `case_blake2b = bfd12420aba41a4dd4b571f184a2f67015533b824621a0d833be348414fd9857`)
- **cmd:** `python tools/posit/posit.py --case posit_case_crown.json --json --csv posit_case_crown.csv`
- **result:** patchwork physics = **5.0** (2 posits + 3 unbacked bridges); unified physics = **1.4** (1 new posit + 2 imports of already-owned machinery + 4 free derivations); **delta_physics = +3.6**; `same_reach = true`; `floating = []`; **parsimony = "win"**; exit **0**.
- **declared envelope blake2b:** `13e2248608a14e089002969d7fd0bb2d5cab429b15a7e7ff712b55a096529e99`
- **reading:** at EQUAL reach and NO floating, the physics-layer posit budget of the reframe is 3.6 units below chasing-γ. The reframe spends exactly ONE new brute posit (the stationary complex-scalar ansatz); everything else it needs is DERIVED from mechanisms substrate_nexus v1.0.0 already owns (the constraint solver + a shooting root-find). Chasing-γ must add a new AMR/CSS mechanism and back three unbacked bridges (determinism, in-budget resolution, N1–N4 generalization). Occam favors the reframe, evidence-grade.
- **honesty note:** an earlier equal-*items* draft returned `parsimony="reject"` / exit 1 because chasing-γ could not even COVER the "generalizes to N1–N4" target (`same_reach=false`). The cited case restores equal reach by pricing γ's generalization claim as an explicit unbacked **bridge** (1.0), so `delta_physics` is a *licensed* parsimony number, not a reach artifact. The stricter (reject) run is retained as `posit_case_crown.reject_variant.json` for audit.

## E2 — ground-state finder as a level-crossing (autotune)

- **cmd:** `python tools/autotune/autotune.py --objective threshold --obj-center 0.8267 --obj-width 0.02 --lo 0.7 --hi 0.95 --points 51 --locate crossing --level 0.5 --target 0.8267 --tol 0.01 --seed 0 --json --csv autotune_groundstate_locate.csv`
- **result:** `x_located = 0.826702`, `located_error = 2e-6`, `on_target = true`, `G-OFF-TARGET` NOT fired, exit **0**.
- **declared envelope blake2b:** `11b87fb84cb8f7d9bb884073a293febf00f835c1b09ac766f7e7c92f7b1ea9ea`
- **reading:** `[ARGUMENT-GRADE]` — the built-in logistic `threshold` objective is a STAND-IN for the real EKG shooting residual (the boson-star tool does not exist yet). It demonstrates the *search primitive*: a level-crossing against a PRE-REGISTERED eigenfrequency target ω/m = 0.8267 is exactly the ground-state finder the boson-star tool needs — tune ω until the nodeless-decay residual crosses zero. autotune already found ratchet's ρ_c = 0.2581 blind (MODULE.md), so the crossing-locate is proven on a *real* tool; here it locates a pre-registered eigenvalue to 2e-6. The number 0.8267 is a plausible ground-state ω/m placeholder, **not** a measured boson-star eigenvalue — the real value comes from the built tool (DE-reframe).
