# probe_wide.py — WIDE real-kappa scan (incl. negative) for ALL zeros of the residual.
# A zero at -2.81 / -0.357 would reveal a sign-convention flip in the kappa term (fixable).
# Uses RAW proj (pre-normalization) sign changes + boundedness dips; robust to the 4k-2 artifact.
import numpy as np, math, pickle, sys
import bg_analyze as B

bg = pickle.load(open(sys.argv[1] if len(sys.argv)>1 else "bg_fast.pkl","rb"))
dstop = float(sys.argv[2]) if len(sys.argv)>2 else 0.006
print(f"# nsteps={bg['nsteps']} finalDson={bg['Dnode'][-1]:.2e} dstop={dstop}")
print(f"{'k':>7} {'sign':>5} {'Re proj':>14} {'log|hp|':>10}")
rows=[]
for k in np.arange(-4.0, 4.01, 0.1):
    if abs(k)<1e-9: k=1e-6
    r = B.integrate(bg, complex(k), dstop=dstop, record_profile=False)
    if r is None: print(f"{k:7.2f}  fail"); rows.append((k,None,None)); continue
    p=r['proj'].real
    rows.append((k,p,r['logabs']))
    print(f"{k:7.2f} {'+' if p>0 else '-':>5} {p:14.5e} {r['logabs']:10.3f}")
print("\nRe(proj) sign changes (candidate eigenvalues):")
for i in range(1,len(rows)):
    a,pa,_=rows[i-1]; b,pb,_=rows[i]
    if pa is not None and pb is not None and pa*pb<0:
        print(f"   between k={a:.2f} and k={b:.2f}")
print("Boundedness (log|hp|) local minima:")
for i in range(1,len(rows)-1):
    if rows[i-1][2] and rows[i][2] and rows[i+1][2]:
        if rows[i][2]<rows[i-1][2] and rows[i][2]<rows[i+1][2]:
            print(f"   min at k={rows[i][0]:.2f}  log|hp|={rows[i][2]:.3f}")
