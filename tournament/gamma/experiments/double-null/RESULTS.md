# double-null advocate — ORRERY experiment record (phase 1, γ fork)

**Tool:** ORRERY `autotune` v1.0.0 (`C:\ORRERY\tools\autotune\autotune.py`), read-only. Driver: `run_autotune.py` (this dir). CWD `C:\ORRERY`.

**What these runs certify.** *Not* γ, and *not* the double-null EKG solver (which does not exist yet — its contract is in `../../phase1/double-null.md`). They certify the **DE-search half** of the claim: that a critical-point locate against a **pre-registered** target is deterministic, byte-reproducible, and **self-certifying** — it fires a gate (exit 1) when handed a wrong p*, so the search cannot silently rubber-stamp a fabricated critical amplitude. This is the `autotune --locate crossing` mechanism the deciding experiment will point at the real collapse flag.

**Provenance checks (run first).**
- `autotune --selftest` → SELFTEST PASS (6/6: blake2b KAT, argmax recover, threshold-crossing recover, off-target fires, determinism).
- `autotune --golden` → `GOLDEN OK blake2b=c79002f23cf236baab5ecdb5753603a7a3853199f750b0139771d4e6cdd55bbe`.

**Runs (declared blake2b computed with autotune's own serializer = the tool's contract hash):**

| tag | command (abridged) | x_located | on_target | gate fired | exit | declared blake2b |
|---|---|---|---|---|---|---|
| charter-arming | `--objective threshold --obj-center 0.5 --obj-width 0.05 --lo 0 --hi 1 --points 41 --locate crossing --level 0.5 --target 0.5 --tol 0.02 --seed 0` | 0.500000 | true | no | 0 | `db490a3111b770996fcec05b9bc59981f35395980be5f1ad617101dd13552c80` |
| DN-locate-pstar | `--objective threshold --obj-center 0.334 --obj-width 0.01 --lo 0.30 --hi 0.37 --points 41 --locate crossing --level 0.5 --target 0.334 --tol 0.005` | 0.334000 | true | no | 0 | `26ec96dea28874a60730c15628eea966d610667e5a56beae5074617cd8c25ac5` |
| DN-locate-wrongtarget | as above but `--target 0.360` | 0.334000 | **false** | **yes (G-OFF-TARGET)** | **1** | `371b72f65d9a25b200886890d5f9bdc131f53302bddcccf5307394c79fd52838` |

**Determinism.** The charter-arming declared object re-hashes byte-identical across two runs (`db490a31…` both times). `run_autotune.py` returns exit 0 iff all three sub-runs behave exactly as pre-registered.

**Honest scope label.** Any γ number in the approach doc is `[ARGUMENT-GRADE]` / `[ARGUMENT-GRADE]` until the real double-null EKG tool is built and its own golden freezes. `obj-center 0.334` is a *stand-in sigmoid location*, chosen only to have a concrete pre-registered target — it is **not** a measured critical amplitude and **not** a γ figure.

**Reproduce:** `Set-Location C:\ORRERY; python C:\TinyUniverse-tournament\tournament\gamma\experiments\double-null\run_autotune.py`
