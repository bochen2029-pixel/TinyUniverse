#!/usr/bin/env python3
# run_autotune.py — double-null advocate's ORRERY-arming driver (phase 1, the gamma fork).
#
# PURPOSE. Make the *search-methodology* half of the double-null claim evidence-grade, not
# argument-grade. The real double-null EKG integrator does not exist yet (its contract is
# specified in ../../phase1/double-null.md). What CAN be shown now, with no GPU and no new
# gear, is that a critical-point locate — the exact ORRERY mechanism the deciding experiment
# DE-gamma will drive to find p* — is DETERMINISTIC and SELF-CERTIFYING against a
# pre-registered target: it locates a level-crossing to the tol, and it FIRES A GATE (exit 1)
# when handed a wrong pre-registered p*, so it cannot silently agree with a fabricated number.
#
# This drives ORRERY's `autotune` (read-only; we never edit ORRERY source). Each run's
# declared blake2b is computed with autotune's OWN serializer, so the hash is the tool's
# contract hash, byte-for-byte.
#
# Run:  (CWD C:\ORRERY)   python <path>\run_autotune.py
# Exit: 0 if all three sub-runs behaved as pre-registered; 1 otherwise.

import sys, os, hashlib

# Locate ORRERY's autotune regardless of CWD (prefer C:\ORRERY).
CANDIDATES = [
    os.path.join(os.environ.get("ORRERY_ROOT", r"C:\ORRERY"), "tools", "autotune"),
    r"C:\ORRERY\tools\autotune",
]
for c in CANDIDATES:
    if os.path.isfile(os.path.join(c, "autotune.py")):
        sys.path.insert(0, c)
        break
import autotune as at  # noqa: E402


def run_and_hash(tag, A):
    r, gs, v, code = at.run_autotune(A, None)
    declared = at.declared_object(0, A, r, gs, v)
    h = hashlib.blake2b(declared.encode("utf-8"), digest_size=32).hexdigest()
    print(f"[{tag}]")
    print(f"    declared_blake2b = {h}")
    print(f"    x_located={r['x_located']:.6f}  on_target={r['on_target']}  "
          f"G-OFF-TARGET.fired={gs[0]['fired']}  exit={code}")
    return h, r, gs, code


def main():
    ok = True

    # (1) CHARTER ARMING — the exact command in CHARTER.md ("ORRERY arming"): a threshold
    #     objective, level-crossing located at the pre-registered 0.5. Proves the loop runs.
    A_charter = {"objective": "threshold", "obj_center": 0.5, "obj_width": 0.05, "tool": "",
                 "sweep": "", "metric": "", "fixed": "", "lo": 0.0, "hi": 1.0, "points": 41,
                 "locate": "crossing", "level": 0.5, "target": 0.5, "tol": 0.02}
    print("--- (1) CHARTER ARMING: crossing located at pre-registered 0.5 ---")
    h1, _, _, c1 = run_and_hash("charter-arming", A_charter)
    ok &= (c1 == 0)

    # (2) POSITIVE — a collapse-fraction sigmoid stepping 0->1 across a pre-registered
    #     p* = 0.334 (a plausible ingoing-pulse critical amplitude in Garfinkle's family).
    #     autotune locates the 0.5-crossing; matching 0.334 to tol => the p*-locate DE-search
    #     mechanism is deterministic and lands where pre-registered. This STANDS IN for
    #     driving the (not-yet-built) double-null EKG tool; the sigmoid is NOT a gamma figure.
    A_pos = {"objective": "threshold", "obj_center": 0.334, "obj_width": 0.01, "tool": "",
             "sweep": "", "metric": "", "fixed": "", "lo": 0.30, "hi": 0.37, "points": 41,
             "locate": "crossing", "level": 0.5, "target": 0.334, "tol": 0.005}
    print("--- (2) POSITIVE: p* located at pre-registered 0.334 (tol 0.005) ---")
    h2, _, gs2, c2 = run_and_hash("DN-locate-pstar", A_pos)
    ok &= (c2 == 0 and not gs2[0]["fired"])

    # (3) NEGATIVE control — same curve, WRONG pre-registered target 0.360. G-OFF-TARGET MUST
    #     fire, exit 1. This is the KILL guard proving the search cannot rubber-stamp a
    #     fabricated p*: the honest-oracle discipline made mechanical.
    A_neg = dict(A_pos); A_neg["target"] = 0.360
    print("--- (3) NEGATIVE control: wrong target 0.360 must fire the gate (exit 1) ---")
    h3, _, gs3, c3 = run_and_hash("DN-locate-wrongtarget", A_neg)
    ok &= (c3 == 1 and gs3[0]["fired"])

    # Determinism: re-run (1) and confirm byte-identical declared hash.
    r, gs, v, _ = at.run_autotune(A_charter, None)
    h1b = hashlib.blake2b(at.declared_object(0, A_charter, r, gs, v).encode("utf-8"),
                          digest_size=32).hexdigest()
    det = (h1b == h1)
    print(f"--- determinism: charter run reproduces byte-identical = {det} ---")
    ok &= det

    print("\nRESULT:", "PASS (all sub-runs behaved as pre-registered)" if ok
          else "FAIL (a sub-run deviated)")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
