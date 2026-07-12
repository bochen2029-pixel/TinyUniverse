# probe_plateau.py — for fixed k, sweep the stop shell dstop and check whether
#   c_na = [wL_n . hp] / (x_s-x)^(1-2k) / (4k-2)
# is INDEPENDENT of dstop (a plateau). If clean (wL_n kills analytic branches), it must be flat.
# Drift => contamination (analytic admixture leaking through) => extraction not isolating c_na.
import numpy as np, math, pickle, sys
import bg_analyze as B

bg = pickle.load(open(sys.argv[1] if len(sys.argv)>1 else "bg_h1e4.pkl","rb"))
ks = [float(a) for a in sys.argv[2:]] or [0.357, 1.0, 2.0, 2.81055, 3.0]
dstops = [0.05,0.03,0.02,0.012,0.008,0.005,0.003,0.002]
print(f"# nsteps={bg['nsteps']} finalDson={bg['Dnode'][-1]:.2e}")
for k in ks:
    print(f"\nk={k}:  c_na vs stop shell (should be flat if extraction is clean)")
    print(f"{'dstop':>8} {'(xs-x)':>10} {'Re c_na':>15} {'Im c_na':>15} {'|c_na|':>13}")
    for ds in dstops:
        if ds < bg['dmin']*1.2: continue
        r = B.integrate(bg, complex(k), dstop=ds, record_profile=False)
        if r is None: print(f"{ds:8.3f}  fail"); continue
        c=r['cna']
        print(f"{ds:8.3f} {r['dxm']:10.4e} {c.real:15.6e} {c.imag:15.6e} {abs(c):13.4e}")
