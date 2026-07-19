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
        r65, m65 = res['cross'].get(0.65, (np.nan, np.nan))
        r70, m70 = res['cross'].get(0.70, (np.nan, np.nan))
        out.append((dp, res['MBH'], res['rH'], res['m2r'], m70, r70, m65, r65))
        print(f"[{i:2d}] dp={dp:.4e}: M70={m70:.5f} r70={r70:.4f} | M65={m65:.5f} | "
              f"Mfrz={res['MBH']:.5f} m2r={res['m2r']:.3f}  t={res['t']:.1f}  ({time.time()-t0:.0f}s)", flush=True)
        np.save('gamma_scaling.npy', np.array(out))
    arr = np.array(out)

    def fit(ycol, rcol, tag, rmin):
        m = np.isfinite(arr[:, ycol]) & (arr[:, rcol] >= rmin)
        x = np.log(arr[m, 0]); y = np.log(arr[m, ycol])
        A = np.vstack([x, np.ones_like(x)]).T
        coef, *_ = np.linalg.lstsq(A, y, rcond=None)
        resid = y - A@coef
        rms = float(np.sqrt(np.mean(resid**2)))
        print(f"*** gamma[{tag}] = {coef[0]:.4f}   (lit 0.374(1); {int(m.sum())} pts "
              f"rH>={rmin:g}, lnM rms {rms:.4f}) ***", flush=True)
        return x, resid, coef[0]

    rmin = 8*ev.dr        # grid-resolvability floor (measured: rH quantizes below ~8 dr)
    print(f"\n[fit] resolvability cut rH >= 8*dr = {rmin:.3f}", flush=True)
    x, resid, gam = fit(4, 5, 'M70', rmin)
    fit(6, 7, 'M65', rmin)
    fit(1, 2, 'Mfrz', rmin)
    # JOINT slope+wiggle fit (the window holds ~1 period, so the wiggle BENDS a plain
    # slope — fit y = c + gamma*x + A cos + B sin together, scanning the period P)
    m = np.isfinite(arr[:, 4]) & (arr[:, 5] >= rmin)
    xj = np.log(arr[m, 0]); yj = np.log(arr[m, 4])
    best = None
    for P in np.linspace(2.0, 8.0, 601):
        B = np.vstack([xj, np.ones_like(xj), np.cos(2*np.pi*xj/P), np.sin(2*np.pi*xj/P)]).T
        cf, *_ = np.linalg.lstsq(B, yj, rcond=None)
        rr = yj - B@cf
        sc = float(rr@rr)
        if best is None or sc < best[1]:
            best = (P, sc, cf, float(np.sqrt(np.mean(rr**2))))
    P, _, cf, rmsj = best
    gj = cf[0]; amp = float(np.hypot(cf[2], cf[3]))
    print(f"*** JOINT[M70]: gamma = {gj:.4f}  period P = {P:.3f}  ->  Delta = 2*gamma*P = "
          f"{2*gj*P:.3f}   (lit gamma 0.374, Delta 3.4453; wiggle amp {amp:.4f}, "
          f"joint rms {rmsj:.4f}, window {(xj.max()-xj.min())/P:.2f} periods) ***", flush=True)

if __name__ == "__main__":
    a = sys.argv[1:]
    main(N=int(a[0]) if len(a) > 0 else 1600,
         lo=float(a[1]) if len(a) > 1 else -4.0,
         hi=float(a[2]) if len(a) > 2 else -2.0,
         n_pts=int(a[3]) if len(a) > 3 else 26)
