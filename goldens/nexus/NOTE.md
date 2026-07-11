# goldens/nexus — freeze record

**Frozen:** 2026-07-11 · **golden.hash:** `ad64f810d44ef542b8c7fe259bd89eeed2548c74ea0ab897bf3d16faa58335e6`
**Params:** contract v0 dial defaults (m_p=1, c=20, ħ=0.5, G=2e-3, dt=1/240, L_box=512, k_B=1) · seed 20260711.
**Platform pin:** MSVC 2022 `/O2` x64, Windows 11, fp64 CPU (nexus is the CPU-double oracle; no GPU pin).
**Determinism evidence:** N11 in-process double-run byte-identical; 3× out-of-process `--golden` → `GOLDEN OK ad64f810` (exit 0) each time.
**Artifact:** `runs/nexus_v1.0.0_freeze.json`.

**Anticipated supersession:** the M1 liborrery envelope lift replaces the self-contained BLAKE2b; if the canonical bytes differ, this golden is re-baselined under a signed entry here (MODULE.md known-limitation #1).

**2026-07-11 (M1) — RESOLVED, GOLDEN STANDS:** `harness/hash_compat.cpp` cross-checked tiny_nexus's self-contained BLAKE2b-256 against liborrery's RFC-KAT'd `orrery::blake2b_hex` on 4 vectors (empty, "abc", fox, 276-byte multi-block): **byte-compatible on all**. No supersession; `ad64f810` remains the frozen golden. MODULE.md known-limitation #1 is closed.
