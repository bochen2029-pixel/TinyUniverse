# harness/ — compile-as-verification

From M0: `verify.py` (pattern: `C:\ORRERY\harness\verify.py`) — compile every module per its MODULE.md command → run every `--selftest` → run every `--golden` → red/green summary, exit 0/1. A claim that cannot compile or fails its golden is not in the project.

Budgets: selftest < 30 s per module; full golden suite < 5 min. Two-pass verification (fresh cold-context agent re-runs contract-conformance + golden with no memory of the build) is mandatory before any result is cited publicly.

**GPU preflight:** `verify.py` checks free VRAM via `nvidia-smi` first and exits **3** (environment-not-ready — never conflated with a red golden) if under ~2 GB free. This machine's card frequently runs other workloads; a contended card gives false REDs (cudaMalloc exit 2) or shared-memory-spilled garbage timings.
