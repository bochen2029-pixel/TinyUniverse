#!/usr/bin/env python3
# harness/verify.py — compile-as-verification (ARCHITECTURE §9).
# Runs every frozen golden cold: tiny_nexus + all app scenarios. Red/green, exit 0/1.
# Run from the repo root. Assumes binaries already built per BUILD.md (build/ dir);
# pass --build to compile first (vcvars64 + cl/nvcc, ~2 min).

import subprocess, sys, os, time

VCVARS = r'"C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat"'
BUILDS = [
    f'{VCVARS} >nul 2>&1 && cl /std:c++17 /EHsc /O2 /W4 /nologo nexus\\tiny_nexus.cpp /Fe:build\\tiny_nexus.exe /Fo:build\\tiny_nexus.obj',
    f'{VCVARS} >nul 2>&1 && nvcc -O3 -arch=sm_89 -Xcompiler "/O2" -o build\\tinyuniverse.exe app\\tinyuniverse.cu core\\lib\\envelope.cpp user32.lib gdi32.lib opengl32.lib cufft.lib',
]
GOLDENS = [
    ("nexus",     [r"build\tiny_nexus.exe", "--golden"]),
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

def main():
    if not gpu_preflight():
        return 3
    if "--build" in sys.argv:
        for cmd in BUILDS:
            print(f"[build] {cmd[:80]}...")
            r = subprocess.run(f'cmd /c \'{cmd}\'' if os.name != "nt" else cmd, shell=True)
            if r.returncode != 0:
                print("[build] FAILED"); return 1
    red = 0
    t0 = time.time()
    for name, cmd in GOLDENS:
        t = time.time()
        r = subprocess.run(cmd, capture_output=True, text=True)
        ok = (r.returncode == 0) and ("GOLDEN OK" in (r.stderr or ""))
        print(f"  {name:<10} {'GREEN' if ok else 'RED':<6} {time.time()-t:6.1f}s"
              + ("" if ok else f"  (exit={r.returncode})"))
        if not ok:
            red += 1
            tail = (r.stderr or "").strip().splitlines()[-3:]
            for ln in tail: print(f"      {ln}")
    print(f"{'-'*40}\n  {'ALL GREEN' if red == 0 else f'{red} RED'}  ({time.time()-t0:.0f}s total)")
    return 0 if red == 0 else 1

if __name__ == "__main__":
    sys.exit(main())
