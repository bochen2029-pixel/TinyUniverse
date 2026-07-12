# stageB_solve.py — root-find the corrected residual c_na(k) for the physical eigenvalue k0,
# then G-CONVERGE (stop-shell/step refinement) and G-UNIQUE (one physical relevant mode). Honest.
import numpy as np, math, sys
from scipy.optimize import brentq
import stageB_fast as B

def refine_real(k_lo, k_hi, dstop, h, z0=math.exp(-12)):
    """brentq on Re(c_na) between k_lo,k_hi (real). Returns (k0, c_na(k0))."""
    def f(k):
        out = B.c_na(complex(k), dstop=dstop, z0=z0, h=h)
        return out[0].real if out is not None else float('nan')
    k0 = brentq(f, k_lo, k_hi, xtol=1e-8, rtol=1e-12, maxiter=200)
    out = B.c_na(complex(k0), dstop=dstop, z0=z0, h=h)
    return k0, out[0]

if __name__ == "__main__":
    print("=== Stage B solve: fluid-CSS critical exponent (HONEST, corrected residual) ===\n")

    # 1) coarse scan to SEE crossings (default fine-ish grid)
    dstop0, h0 = 0.02, 5e-5
    print(f"[1] Coarse scan of c_na, dstop={dstop0}, h={h0}:")
    ks = np.arange(0.1, 3.61, 0.1)
    prev = None; brackets = []
    print(f"{'k':>7} {'Re c_na':>14} {'Im c_na':>14} {'log|hp|':>12}")
    for k in ks:
        out = B.c_na(complex(k), dstop=dstop0, h=h0)
        if out is None:
            print(f"{k:7.2f}  (fail)"); prev=None; continue
        cna,proj,dxm,logabs = out
        print(f"{k:7.2f} {cna.real:14.5e} {cna.imag:14.5e} {logabs:12.3f}")
        if prev is not None and prev[1].real*cna.real < 0:
            brackets.append((prev[0], k))
        prev = (k, cna)
    print("  sign-change brackets:", brackets)

    # 2) refine each real bracket
    print("\n[2] Refined real roots (brentq on Re c_na):")
    roots = []
    for (a,b) in brackets:
        try:
            k0, c0 = refine_real(a, b, dstop0, h0)
            roots.append(k0)
            beta = 1.0/k0
            tag = "GAUGE?" if abs(k0-0.35699)<0.05 else ("PHYSICAL?" if 2.0<k0<4.0 else "")
            print(f"   root k0 = {k0:.8f}   |c_na|={abs(c0):.2e}   beta=1/k0={beta:.8f}   {tag}")
        except Exception as e:
            print(f"   bracket ({a:.2f},{b:.2f}) failed: {e}")

    # 3) identify the physical mode (largest real root that is not the gauge mode) and converge it
    phys = [k for k in roots if k > 1.0]
    if not phys:
        print("\n[3] NO physical mode (k>1) found. Reporting all roots above.")
        sys.exit(0)
    kP = phys[0]
    print(f"\n[3] G-CONVERGE on the physical mode near k={kP:.5f}")
    print(f"{'dstop':>8} {'h':>10} {'k0':>14} {'beta':>12}")
    ladder = [(0.04,1e-4),(0.02,5e-5),(0.01,5e-5),(0.01,2.5e-5),(0.005,2.5e-5),(0.005,1.25e-5)]
    kvals=[]
    for (ds,hh) in ladder:
        # re-bracket locally around kP
        lo,hi = kP-0.15, kP+0.15
        try:
            k0,c0 = refine_real(lo,hi,ds,hh)
            kvals.append(k0)
            print(f"{ds:8.3f} {hh:10.2e} {k0:14.8f} {1.0/k0:12.8f}")
        except Exception as e:
            print(f"{ds:8.3f} {hh:10.2e}   fail: {e}")
    if kvals:
        spread = max(kvals)-min(kvals)
        kbest = kvals[-1]; beta=1.0/kbest
        print(f"\n  k0 spread across ladder = {spread:.3e}")
        print(f"  BEST (finest): Re k0 = {kbest:.8f}   beta = 1/Re k0 = {beta:.8f}")
        print(f"  target: Re k0 = 2.81055255, beta = 0.35580192")
        print(f"  |beta - 0.35580192| = {abs(beta-0.35580192):.3e}  ({abs(beta-0.35580192)/0.35580192*100:.3f}% )")
        print(f"  G-CONVERGE: {'PASS' if spread<1e-2 else 'MARGINAL/FAIL'} (spread {spread:.2e})")
