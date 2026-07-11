# goldens/ — frozen behavior

One directory per module/scenario: `goldens/<name>/` holds `(dials, seed, trace, steps) → blake2b-256 of declared output`, plus a NOTE.md when a golden is superseded (operator-signed, with why).

Rules: hardware-pinned sm_89 for CUDA modules (re-baseline = signed NOTE.md); MSVC-pinned for CPU modules (nexus); visuals and timings are never in the hash domain; a rewrite that cannot reproduce its golden is not a rewrite, it is a proposal.
