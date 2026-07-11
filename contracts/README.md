# contracts/ — the contract discipline

A contract is the sacred, semver'd interface of a module or scenario: CLI + types + ranges, declared-output schema, exit-code semantics, determinism clause, oracle. **No source until the contract exists and is operator-reviewed.** Breaking change = MAJOR bump + migration note + golden supersession under review.

Three contract kinds in this repo:

1. **Module contracts** (`<module>.contract.md` [+ `.schema.json` if it has a headless face]) — nexus, newton, arrow, einstein, gargantua, planck, cosmos.
2. **Scenario contracts** — a named initial condition + dial set, runnable headless; these are the golden units. Scenarios ride the universal envelope below.
3. **THE frame contract** (`frame.contract.md`) — the sim↔renderer seam. Singular, and the most breaking-change-averse file in the repo.

## The universal envelope (every headless face obeys)

```
tinyuniverse.exe --scenario NAME [--dials PATH] --steps N --seed S [--trace PATH]
                 [--json | --csv PATH] [--selftest] [--golden]
tiny_nexus.exe   [--dials PATH] [--seed S] [--json] [--selftest] [--golden]

stdout (--json): one JSON object matching the module's schema:
  { "tool":"...", "version":"1.0.0", "seed":S, "dials":{...}, "params":{...},
    "results":{...declared fields...}, "verdict":"pass|fail", "notes":"..." }

exit: 0 pass · 1 declared-gate-fired (real negative result) · 2 error   (never conflated)
--selftest : internal battery, print PASS/FAIL, exit 0/1
--golden   : run frozen params, print declared-output hash, exit 0/1 vs goldens/
determinism: (dials, seed, trace, steps) → byte-identical declared JSON
             ("notes", timings, and all visuals excluded from the hash domain)
```

Hash domain and serializers: liborrery envelope (blake2b-256 over canonical `{seed,params,results,gates,verdict}`, floats `%.6f`, −0 normalized) — D-005, byte-compatible with ORRERY's.
