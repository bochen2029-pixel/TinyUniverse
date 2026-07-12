#!/usr/bin/env python3
# harness/verify.py — compile-as-verification (ARCHITECTURE §9).
# Runs every frozen golden cold: tiny_nexus + all app scenarios. Red/green, exit 0/1.
# Run from the repo root. Assumes binaries already built per BUILD.md (build/ dir);
# pass --build to compile first (vcvars64 + cl/nvcc, ~2 min).

import subprocess, sys, os, time

VCVARS = r'"C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat"'
BUILDS = [
    f'{VCVARS} >nul 2>&1 && cl /std:c++17 /EHsc /O2 /W4 /nologo nexus\\tiny_nexus.cpp /Fe:build\\tiny_nexus.exe /Fo:build\\tiny_nexus.obj',
    f'{VCVARS} >nul 2>&1 && cl /std:c++17 /EHsc /O2 /W4 /nologo substrate\\substrate_nexus.cpp /Fe:build\\substrate_nexus.exe /Fo:build\\substrate_nexus.obj',
    f'{VCVARS} >nul 2>&1 && cl /std:c++17 /EHsc /O2 /W4 /nologo substrate\\curve_nexus.cpp /Fe:build\\curve_nexus.exe /Fo:build\\curve_nexus.obj',   # v2 N3 (curve, CPU fp64)
    f'{VCVARS} >nul 2>&1 && nvcc -O3 -arch=sm_89 -Xcompiler "/O2" -o build\\tinyuniverse.exe app\\tinyuniverse.cu core\\lib\\envelope.cpp user32.lib gdi32.lib opengl32.lib cufft.lib',
    f'{VCVARS} >nul 2>&1 && nvcc -O3 -arch=sm_89 -Xcompiler "/O2" -o build\\field_nexus.exe substrate\\field_nexus.cu cufft.lib',   # v2 N1 (SP field)
    f'{VCVARS} >nul 2>&1 && nvcc -O3 -arch=sm_89 -Xcompiler "/O2" -o build\\lapse_nexus.exe substrate\\lapse_nexus.cu cufft.lib',   # v2 N2 (lapse/clock)
]
# CPU fp64 oracles (no GPU) run regardless of card contention; GPU goldens follow the preflight.
CPU_GOLDENS = [
    ("nexus",           [r"build\tiny_nexus.exe", "--golden"]),
    ("substrate_nexus", [r"build\substrate_nexus.exe", "--golden"]),   # v2 N0 (spherical EKG)
]
# v2 N3 curve (geometry: geodesics through the weak-field metric -> exact-GR light bending +
# precession). CPU fp64 geodesic oracle (no GPU) -> runs with the CPU oracles, no preflight.
CURVE_GOLDENS = [
    ("curve_deflect", [r"build\curve_nexus.exe", "--scenario", "deflect", "--golden"]),
    ("curve_precess", [r"build\curve_nexus.exe", "--scenario", "precess", "--golden"]),
    ("curve_shapiro", [r"build\curve_nexus.exe", "--scenario", "shapiro", "--golden"]),
]
GOLDENS = [
    ("kepler",    [r"build\tinyuniverse.exe", "--scenario", "kepler",    "--golden"]),
    ("threebody", [r"build\tinyuniverse.exe", "--scenario", "threebody", "--golden"]),
    ("cloud",     [r"build\tinyuniverse.exe", "--scenario", "cloud",     "--golden"]),
    ("galaxy",    [r"build\tinyuniverse.exe", "--scenario", "galaxy",    "--golden"]),
    ("merger",    [r"build\tinyuniverse.exe", "--scenario", "merger",    "--golden"]),
    ("echo",      [r"build\tinyuniverse.exe", "--scenario", "echo",      "--golden"]),
    ("ratchet",   [r"build\tinyuniverse.exe", "--scenario", "ratchet",   "--golden"]),
    ("detector",  [r"build\tinyuniverse.exe", "--scenario", "detector",  "--golden"]),
    ("keprel",    [r"build\tinyuniverse.exe", "--scenario", "keprel",    "--golden"]),
    ("clocks",    [r"build\tinyuniverse.exe", "--scenario", "clocks",    "--golden"]),
    ("photons",   [r"build\tinyuniverse.exe", "--scenario", "photons",   "--golden"]),
    ("collapse",  [r"build\tinyuniverse.exe", "--scenario", "collapse",  "--golden"]),
    ("isco",      [r"build\tinyuniverse.exe", "--scenario", "isco",      "--golden"]),
    ("hawking",   [r"build\tinyuniverse.exe", "--scenario", "hawking",   "--golden"]),
    ("doubleslit",[r"build\tinyuniverse.exe", "--scenario", "doubleslit","--golden"]),
    ("tunneling", [r"build\tinyuniverse.exe", "--scenario", "tunneling", "--golden"]),
    ("shoq",      [r"build\tinyuniverse.exe", "--scenario", "shoq",      "--golden"]),
    ("circumnav", [r"build\tinyuniverse.exe", "--scenario", "circumnav", "--golden"]),
    ("expand",    [r"build\tinyuniverse.exe", "--scenario", "expand",    "--golden"]),
    ("bubbles",   [r"build\tinyuniverse.exe", "--scenario", "bubbles",   "--golden"]),
]
# v2 N1 field (Schrodinger-Poisson on a 3D lattice) — GPU goldens behind the same
# preflight. freepacket/sho3d are oracle-grounded (vs nexus N5, exact); sho3d
# re-runs 20k imaginary-time iters (~2 min). soliton (the weld) is added once frozen.
FIELD_GOLDENS = [
    ("field_freepacket", [r"build\field_nexus.exe", "--scenario", "freepacket", "--golden"]),
    ("field_sho3d",      [r"build\field_nexus.exe", "--scenario", "sho3d",      "--golden"]),
    # soliton = the PM+psi gravity weld (self-bound SP soliton; r_c*M scale-covariant
    # to 3e-8 across a 2:1 mass pair -> quantum pressure balances self-gravity). 256^3,
    # 2x10000 imaginary-time iters -> ~7 min; the load-bearing new-physics golden.
    ("field_soliton",    [r"build\field_nexus.exe", "--scenario", "soliton",    "--golden"]),
    # echoF = determinism receipt: time-reversal by conjugation (reversible to fp32
    # round-off) + byte-exact reproducibility. 128^3, ~2 s.
    ("field_echoF",      [r"build\field_nexus.exe", "--scenario", "echoF",      "--golden"]),
    # cloudF = classical-limit weld cross-check: an overdense sphere collapses under
    # self-gravity (by the field) at the analytic free-fall time t_ff (v1 cloud/collapse
    # physics; the Madelung Q->0 limit). 128^3, ~6 s.
    ("field_cloudF",     [r"build\field_nexus.exe", "--scenario", "cloudF",     "--golden"]),
    # mergerF = two-body weld cross-check: two self-gravitating psi-lumps attract under
    # mutual gravity (separation shrinks ~2x) and form a denser remnant (v1 merger
    # physics). Proves gravity BETWEEN distinct masses is real. 128^3, ~8 s.
    ("field_mergerF",    [r"build\field_nexus.exe", "--scenario", "mergerF",    "--golden"]),
]
# v2 N2 lapse (the clock: alpha = sqrt(1+2 Phi/c^2), proper time tau = integral alpha dt).
# GPU goldens behind the same preflight (both ~1-2 s). redshift = exact Schwarzschild
# gravitational time dilation (analytic oracle); redshiftPM = the substrate weld (the
# PM-Poisson well of Newtonian depth, through the lapse -> correct gravitational redshift).
LAPSE_GOLDENS = [
    ("lapse_redshift",   [r"build\lapse_nexus.exe", "--scenario", "redshift",   "--golden"]),
    ("lapse_redshiftPM", [r"build\lapse_nexus.exe", "--scenario", "redshiftPM", "--golden"]),
]

def gpu_preflight(min_free_mb=2000):
    # The app needs ~0.8 GB VRAM; a card busy with another workload would fail
    # cudaMalloc (false RED) or spill to shared memory (garbage timings).
    # Exit 3 = environment-not-ready: NOT a red golden, never conflate.
    try:
        r = subprocess.run(["nvidia-smi", "--query-gpu=memory.free",
                            "--format=csv,noheader,nounits"],
                           capture_output=True, text=True, timeout=10)
        free_mb = int(r.stdout.strip().splitlines()[0])
    except Exception:
        print("[preflight] nvidia-smi unavailable - proceeding blind")
        return True
    if free_mb < min_free_mb:
        print(f"[preflight] GPU busy: only {free_mb} MB VRAM free (< {min_free_mb} MB).")
        print("[preflight] Refusing to run goldens on a contended card. Exit 3.")
        return False
    print(f"[preflight] GPU ok: {free_mb} MB VRAM free")
    return True

def run_goldens(goldens):
    red = 0
    for name, cmd in goldens:
        t = time.time()
        r = subprocess.run(cmd, capture_output=True, text=True)
        ok = (r.returncode == 0) and ("GOLDEN OK" in (r.stderr or ""))
        print(f"  {name:<16} {'GREEN' if ok else 'RED':<6} {time.time()-t:6.1f}s"
              + ("" if ok else f"  (exit={r.returncode})"))
        if not ok:
            red += 1
            tail = (r.stderr or "").strip().splitlines()[-3:]
            for ln in tail: print(f"      {ln}")
    return red

def main():
    if "--build" in sys.argv:
        for cmd in BUILDS:
            print(f"[build] {cmd[:80]}...")
            r = subprocess.run(f'cmd /c \'{cmd}\'' if os.name != "nt" else cmd, shell=True)
            if r.returncode != 0:
                print("[build] FAILED"); return 1
    t0 = time.time()
    print("[cpu] fp64 oracles (GPU-independent):")
    red = run_goldens(CPU_GOLDENS)              # nexus + substrate_nexus run under any contention
    red += run_goldens(CURVE_GOLDENS)           # v2 N3 curve (CPU fp64 geodesic oracle) — GPU-independent
    gpu_ok = gpu_preflight()
    if gpu_ok:
        print("[gpu] v1 scenario goldens:")
        red += run_goldens(GOLDENS)
        print("[gpu] v2 N1 field goldens:")
        red += run_goldens(FIELD_GOLDENS)
        print("[gpu] v2 N2 lapse goldens:")
        red += run_goldens(LAPSE_GOLDENS)
    n = len(CPU_GOLDENS) + len(CURVE_GOLDENS) + ((len(GOLDENS) + len(FIELD_GOLDENS) + len(LAPSE_GOLDENS)) if gpu_ok else 0)
    print("-"*40)
    if not gpu_ok:                              # CPU oracles still reported; GPU suite deferred
        print(f"  CPU {'ALL GREEN' if red==0 else f'{red} RED'}; GPU goldens SKIPPED "
              f"(card contended, exit 3)  ({time.time()-t0:.0f}s)")
        return 1 if red else 3
    print(f"  {'ALL GREEN' if red == 0 else f'{red} RED'}  ({n} goldens, {time.time()-t0:.0f}s total)")
    return 0 if red == 0 else 1

if __name__ == "__main__":
    sys.exit(main())
