# probe_gauge.py — does the operator have the KHA gauge mode at k~0.35699? A clean zero there
# validates the perturbation operator. Also fine-scan around 2.81. Use RAW proj (pre-4k-2) and the
# boundedness minimum (log|hp| local min) which is normalization-free.
import numpy as np, math, pickle, sys
import bg_analyze as B

bg = pickle.load(open(sys.argv[1] if len(sys.argv)>1 else "bg_h1e4.pkl","rb"))
dstop = float(sys.argv[2]) if len(sys.argv)>2 else 0.004

def raw(k):
    r = B.integrate(bg, complex(k), dstop=dstop, record_profile=False)
    if r is None: return None
    return r  # dict with cna, proj, dxm, logabs

def scan(ks, label):
    print(f"\n== {label} (dstop={dstop}) ==")
    print(f"{'k':>9} {'Re proj':>14} {'Im proj':>14} {'log|hp|':>11}")
    rows=[]
    for k in ks:
        r=raw(k)
        if r is None: print(f"{k:9.4f}  fail"); continue
        rows.append((k, r['proj'], r['logabs']))
        print(f"{k:9.4f} {r['proj'].real:14.5e} {r['proj'].imag:14.5e} {r['logabs']:11.3f}")
    # sign changes in Re(proj)
    sc=[]
    for i in range(1,len(rows)):
        if rows[i-1][1].real*rows[i][1].real<0: sc.append((rows[i-1][0],rows[i][0]))
    print("  Re(proj) sign changes:", sc)
    # boundedness minima
    mn=[]
    for i in range(1,len(rows)-1):
        if rows[i][2]<rows[i-1][2] and rows[i][2]<rows[i+1][2]: mn.append(rows[i][0])
    print("  log|hp| local minima:", mn)
    return rows

scan(np.arange(0.20, 0.55, 0.02), "gauge-mode region (~0.357)")
scan(np.arange(2.60, 3.05, 0.02), "physical-mode region (~2.81)")
scan(np.arange(0.80, 1.25, 0.02), "the observed ~1.0 crossing")
