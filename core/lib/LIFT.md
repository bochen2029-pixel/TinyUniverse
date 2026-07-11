# core/lib — the liborrery lift (D-005)

**Lifted:** 2026-07-11, **verbatim**, from `C:\ORRERY\lib` at ORRERY commit **`d56c4c703df56fd6e1ec406b4eb051489105bfc9`**. Documentation of the invariants lives in `C:\ORRERY\lib\MODULE.md` (referenced, not copied — ORRERY is never edited from here).

Per-file SHA-256 verified identical source↔copy at lift time:

| file | sha256 (first 16) |
|---|---|
| envelope.h | 6AA11C7AFF8B0CF3 |
| envelope.cpp | E89700D1BC4A9B8A |
| rng.cuh | 85E91C5B0F27DB1E |
| reduce.cuh | A43824FA0B9894D5 |
| regime.h | 1815863E94CFB00E |
| ckpt.h | 3A7B8965C4F432CC |
| selftest.cu | CAA87AA1A0CEB88E |

**Verbatim rule (Invariant 7):** any change to these files is a fork — named, justified in DECISIONS.md, and golden-superseding for every consumer. Re-lifting from a newer ORRERY commit updates this file's pin + hashes.

**Verification:** `selftest.cu` is ORRERY's KAT battery (42 checks); its golden-integration tail expects ORRERY's `goldens/` layout, so the canonical green run lives in ORRERY. TU-side verification: `harness/hash_compat.cpp` compiles `envelope.cpp` in this repo and cross-checks BLAKE2b-256 against tiny_nexus's self-contained implementation (result recorded in `goldens/nexus/NOTE.md`).

**First consumers:** `app/tinyuniverse.cu` (rng.cuh — counter-keyed init), `harness/hash_compat.cpp` (envelope).
