# probe_profile.py — model-light eigenvalue probe: measure the local power-law slope of |hp| vs
# (x_s - x) as the perturbation approaches the sonic point.
#   Off-eigenvalue:  |hp| ~ (x_s-x)^(1-2k)  => d log|hp| / d log(x_s-x) -> (1-2k).
#   At an eigenvalue: the non-analytic branch is absent -> |hp| ~ (x_s-x)^0 (bounded), slope -> 0
#     (or a positive analytic slope), i.e. slope DEPARTS from (1-2k) toward 0.
# So define  q(k) = measured_slope(k) - (1-2k). q(k) ~ 0 off-eigenvalue; q(k) -> (2k-1) AT eigenvalue.
# Equivalently the eigenvalue kills the divergence; track (measured_slope) and see where it != (1-2k).
import numpy as np, math, pickle, sys
import bg_analyze as B

bg = pickle.load(open(sys.argv[1] if len(sys.argv)>1 else "bg_h1e4.pkl","rb"))
print(f"# background nsteps={bg['nsteps']} finalDson={bg['Dnode'][-1]:.2e}")
print(f"# probe: local slope of log|hp| vs log(x_s-x) near sonic vs expected (1-2k)")
print(f"{'k':>7} {'meas_slope':>12} {'1-2k':>10} {'slope-(1-2k)':>14}")
for k in [0.2,0.30,0.357,0.40,0.5,0.7,1.0,1.5,2.0,2.5,2.7,2.81055,2.9,3.2]:
    out = B.integrate(bg, complex(k), dstop=bg['dmin']*1.5, record_profile=True)
    if out is None or 'profile' not in out:
        print(f"{k:7.3f}   (fail)"); continue
    prof = out['profile']   # columns: (x_s-x), log|hp|
    if len(prof) < 50: print(f"{k:7.3f}  (short prof {len(prof)})"); continue
    dxm = prof[:,0].astype(float); logabs = prof[:,1].astype(float)
    # fit slope over the last decade approaching sonic (small dxm), where the leading branch dominates
    mask = (dxm>0)&(dxm< np.percentile(dxm,5))   # closest 5% to sonic
    if mask.sum()<10: mask = dxm < np.percentile(dxm,15)
    X = np.log(dxm[mask]); Ylog = logabs[mask]
    A = np.polyfit(X, Ylog, 1)
    slope = A[0]
    exp = 1-2*k.real if hasattr(k,'real') else 1-2*k
    print(f"{k:7.3f} {slope:12.4f} {exp:10.4f} {slope-exp:14.4f}")
