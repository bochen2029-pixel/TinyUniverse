# harness/ — compile-as-verification

From M0: `verify.py` (pattern: `C:\ORRERY\harness\verify.py`) — compile every module per its MODULE.md command → run every `--selftest` → run every `--golden` → red/green summary, exit 0/1. A claim that cannot compile or fails its golden is not in the project.

Budgets: selftest < 30 s per module; full golden suite < 5 min. Two-pass verification (fresh cold-context agent re-runs contract-conformance + golden with no memory of the build) is mandatory before any result is cited publicly.
