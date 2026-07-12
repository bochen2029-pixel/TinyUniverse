# probe_windows.py — measure the local approach slope of log|hp| in SEVERAL distance windows,
# to separate the true Frobenius branch (intermediate window) from background-integration error
# (very close to sonic, where det~Dson->0 corrupts the background RK4).
import numpy as np, math, pickle, sys
import bg_analyze as B

bg = pickle.load(open(sys.argv[1] if len(sys.argv)>1 else "bg_h1e4.pkl","rb"))
print(f"# nsteps={bg['nsteps']} finalDson={bg['Dnode'][-1]:.2e}")
# windows in (x_s - x): far, mid, near
wins = [("far", 0.05, 0.15), ("mid", 0.01, 0.04), ("near", 0.003, 0.009), ("vnear", 0.0012, 0.0025)]
hdr = "     k    1-2k  " + "  ".join(f"{w[0]:>8}" for w in wins)
print(hdr)
for k in [0.357,0.5,0.7,1.0,1.5,2.0,2.5,2.81055,3.0]:
    out = B.integrate(bg, complex(k), dstop=bg['dmin']*1.3, record_profile=True)
    if out is None or 'profile' not in out: print(f"{k:7.3f}  fail"); continue
    prof = out['profile']; dxm=prof[:,0].astype(float); logabs=prof[:,1].astype(float)
    slopes=[]
    for (_,lo,hi) in wins:
        m=(dxm>=lo)&(dxm<=hi)
        if m.sum()>=8:
            A=np.polyfit(np.log(dxm[m]),logabs[m],1); slopes.append(A[0])
        else: slopes.append(float('nan'))
    print(f"{k:7.3f} {1-2*k:7.3f}  " + "  ".join(f"{s:8.3f}" for s in slopes))
print("\nExpectation: in the window where the leading (1-2k) branch dominates, slope=(1-2k).")
print("If slope stays much milder than (1-2k) in ALL windows -> the divergent branch amplitude is")
print("tiny (c_na~0) OR the seed lacks it; if slope MATCHES (1-2k) in mid/near -> branch present.")
