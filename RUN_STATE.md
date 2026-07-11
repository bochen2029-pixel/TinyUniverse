# RUN_STATE.md

**As of:** 2026-07-11 · **Milestone:** M0 `nexus` · **State:** scaffold complete, spec-first — **no module source exists, by design.**

## Current task

M0 gate: **operator review** of `contracts/nexus.contract.md` (v1.0.0-rc1) — specifically:
1. The dial defaults v0 (ARCHITECTURE §6: c=20, ħ=0.5, G=2e-3, m=1, dt=1/240, L=512) — nexus N1 asserts their consequences; changing them later is a MINOR bump until first golden, MAJOR after.
2. The battery N1–N11 — anything missing that a physics claim will later lean on?
3. `contracts/frame.contract.md` v0.1.0-DRAFT — shape check only (freezes at M1).

## Next concrete action

On review approval: implement `nexus/tiny_nexus.cpp` per contract (single file, C++17 fp64, MSVC/clang/g++ parity), battery green, freeze golden, commit as save point, advance RUN_STATE to M1.

## Standing context

Repo founded from the 2026-07-11 brainstorm sessions (concept → rendering strategy → quarry verification). Quarry paths and lift list: README.md table. Sibling doctrine source: `C:\ORRERY`. Agent memory for this project also exists outside the repo (Claude memory dir) — the repo is authoritative where they disagree.
