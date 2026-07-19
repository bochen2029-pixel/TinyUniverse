# gamma_scaling.py — THE DIRECT CLASSIC MEASUREMENT: Choptuik mass scaling on the
# validated evolver. M_BH ∝ (p−p*)^γ for supercritical p; γ = the log-log slope; the
# fine-structure wiggle in the residual (period Δ/(2γ) in ln(p−p*), Gundlach/Hod–Piran)
# is an INDEPENDENT Δ handle. Motivation (measured, endgame2.py + RESULTS_dss.md): the
# BVP route is plateau-walled at |r| ~ 0.08 across seven solver formulations — healthy
# structure, no descent. The scaling route needs no BVP at all: the same near-critical
# evolver that gave p* to 1e-14 and Δ_echo 3.216→3.334 does it by direct evolution.
# Honesty bounds: uniform grid (D-021) — the resolvable window is δp ∈ [1e-4, 1e-2] at
# N=1600 (smallest r_H ~ 4·dr); points under-resolved bend the tail and are visible in
# the residuals. Fixed freeze threshold ⇒ multiplicative bias only ⇒ the slope is clean.
# usage: python gamma_scaling.py [N] [lo] [hi] [npts]
import numpy as np, time, sys
from nr_evolve import Ev

PSTAR = 0.03732817692976     # N=1600 (session-1 bisection, rel ~4e-14)

def main(N=1600, lo=-4.0, hi=-2.0, n_pts=26):
    ev = Ev(N=N, rmax=60.0)
    print(f"[scal] N={N}  p*={PSTAR}  log10(dp) in [{lo},{hi}]  {n_pts} pts", flush=True)
    out = []
    for i, ldp in enumerate(np.linspace(lo, hi, n_pts)):
        dp = 10.0**ldp
        t0 = time.time()
        res = ev.run(PSTAR + dp)
        if res['fate'] != 'bh':
            print(f"[{i:2d}] dp={dp:.4e}: fate={res['fate']} (SKIP)", flush=True)
            continue
        out.append((dp, res['MBH'], res['rH'], res['m2r']))
        print(f"[{i:2d}] dp={dp:.4e}: MBH={res['MBH']:.5f}  rH={res['rH']:.4f}  "
              f"m2r={res['m2r']:.3f}  t={res['t']:.1f}  ({time.time()-t0:.0f}s)", flush=True)
        np.save('gamma_scaling.npy', np.array(out))
    arr = np.array(out)
    x = np.log(arr[:, 0]); y = np.log(arr[:, 1])
    A = np.vstack([x, np.ones_like(x)]).T
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    gam = coef[0]
    resid = y - A@coef
    rms = float(np.sqrt(np.mean(resid**2)))
    print(f"\n*** gamma_scaling = {gam:.4f}   (lit 0.374(1); {len(x)} pts, "
          f"lnM rms resid {rms:.4f}) ***", flush=True)
    # fine structure: periodic residual in x = ln(dp), period P = Delta/(2 gamma)
    best = None
    for P in np.linspace(2.0, 8.0, 601):
        B = np.vstack([np.cos(2*np.pi*x/P), np.sin(2*np.pi*x/P), np.ones_like(x)]).T
        cf, *_ = np.linalg.lstsq(B, resid, rcond=None)
        rr = resid - B@cf
        sc = float(rr@rr)
        if best is None or sc < best[1]:
            best = (P, sc, float(np.hypot(cf[0], cf[1])))
    P, _, amp = best
    print(f"*** fine-structure: best period in ln(dp) = {P:.3f}  ->  Delta = 2*gamma*P = "
          f"{2*gam*P:.3f}   (lit 3.4453; wiggle amplitude {amp:.4f}; window holds "
          f"{(x.max()-x.min())/P:.2f} periods) ***", flush=True)

if __name__ == "__main__":
    a = sys.argv[1:]
    main(N=int(a[0]) if len(a) > 0 else 1600,
         lo=float(a[1]) if len(a) > 1 else -4.0,
         hi=float(a[2]) if len(a) > 2 else -2.0,
         n_pts=int(a[3]) if len(a) > 3 else 26)
