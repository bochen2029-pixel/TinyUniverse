# stageB_qr_solve.py — root-find the QR-match residual for the physical eigenvalue k0, then
# G-CONVERGE (k0 stable as dm and h refine) and G-UNIQUE (one physical mode besides the gauge mode).
# Uses stageB_qr.residual (det[c|Q]); root-finds Re(det) along real k (the mode is real).
import numpy as np, math, sys
from scipy.optimize import brentq
import stageB_qr as B

def det_re(bg, k):
    out = B.residual(bg, complex(k))
    if out is None: return float('nan')
    return out[0].real

def find_root(bg, klo, khi):
    f = lambda k: det_re(bg, k)
    k0 = brentq(f, klo, khi, xtol=1e-9, rtol=1e-13, maxiter=200)
    out = B.residual(bg, complex(k0))
    return k0, out[0], out[1]

if __name__ == "__main__":
    z0 = math.exp(-12)
    print("=== Stage B (QR two-sided match): fluid-CSS eigenvalue kappa_0 -> beta ===\n")

    # 1) coarse scan to LOCATE sign-change brackets (visible, honest)
    h0, dm0 = 5e-5, 0.04
    bg = B.get_bg(B.A2STAR, z0, h0, dm0)
    if bg is None: print("no sonic shell"); sys.exit(2)
    print(f"[1] coarse scan  h={h0} dm={dm0} i_m={bg['i_m']}")
    ks = np.arange(0.2, 4.01, 0.05)
    prev = None; brackets = []
    for k in ks:
        out = B.residual(bg, complex(k))
        if out is None: prev = None; continue
        d = out[0].real
        if prev is not None and prev[1]*d < 0:
            brackets.append((prev[0], k))
        prev = (k, d)
    print("   Re(det) sign-change brackets:", [(round(a,3),round(b,3)) for a,b in brackets])

    # 2) refine every bracket
    print("\n[2] refined roots:")
    roots = []
    for (a, b) in brackets:
        try:
            k0, d0, rej0 = find_root(bg, a, b)
            roots.append(k0)
            tag = "GAUGE" if abs(k0 - 0.357) < 0.03 else ("PHYSICAL" if 2.0 < k0 < 4.0 else "?")
            print(f"   k0 = {k0:.8f}   |det|={abs(d0):.2e}  ||rej||={rej0:.2e}  beta=1/k0={1.0/k0:.8f}   {tag}")
        except Exception as e:
            print(f"   bracket ({a:.3f},{b:.3f}) failed: {e}")

    # 3) G-UNIQUE: report all physical (k>1) roots; expect exactly one
    phys = sorted([k for k in roots if k > 1.0])
    gauge = [k for k in roots if k < 1.0]
    print(f"\n[3] G-UNIQUE: gauge modes(k<1)={[round(g,5) for g in gauge]}  physical modes(k>1)={[round(p,6) for p in phys]}")
    if not phys:
        print("   NO physical mode found. Reporting brackets/roots above (honest)."); sys.exit(1)
    kP = phys[0]

    # 4) G-CONVERGE: refine (dm, h) ladder around kP
    print(f"\n[4] G-CONVERGE on physical mode near k={kP:.5f}")
    print(f"{'dm':>8} {'h':>10} {'i_m':>7} {'k0':>14} {'beta':>12} {'||rej||':>11}")
    ladder = [(0.06, 1e-4), (0.04, 5e-5), (0.03, 5e-5), (0.02, 5e-5), (0.02, 2.5e-5), (0.015, 2.5e-5)]
    kvals = []
    for (dm, h) in ladder:
        bgc = B.get_bg(B.A2STAR, z0, h, dm)
        if bgc is None:
            print(f"{dm:8.3f} {h:10.2e}   no shell"); continue
        try:
            lo, hi = kP - 0.2, kP + 0.2
            k0, d0, rej0 = find_root(bgc, lo, hi)
            kvals.append(k0)
            print(f"{dm:8.3f} {h:10.2e} {bgc['i_m']:7d} {k0:14.8f} {1.0/k0:12.8f} {rej0:11.3e}")
        except Exception as e:
            print(f"{dm:8.3f} {h:10.2e}   fail: {e}")
    if kvals:
        spread = max(kvals) - min(kvals)
        kbest = kvals[-1]; beta = 1.0/kbest
        print(f"\n  k0 spread across ladder = {spread:.3e}")
        print(f"  BEST (finest): Re k0 = {kbest:.8f}   beta = 1/Re k0 = {beta:.8f}")
        print(f"  target:        Re k0 = 2.81055255   beta = 0.35580192")
        print(f"  |beta - target| = {abs(beta-0.35580192):.3e}  ({abs(beta-0.35580192)/0.35580192*100:.4f}%)")
        print(f"  G-CONVERGE: {'PASS' if spread < 1e-2 else 'MARGINAL/FAIL'} (spread {spread:.2e})")
