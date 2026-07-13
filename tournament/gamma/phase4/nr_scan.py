# nr_scan.py — scan the two-sided Evans |E(kappa)| for its eigenvalue dip, matching further from the
# sonic point (wider dip). Report |E| minima; the physical fluid mode should be kappa=2.81055255.
import numpy as np, sys, time
import nr_evans2 as EV

t_match = float(sys.argv[1]) if len(sys.argv) > 1 else -1.0
ks = np.round(np.arange(0.5, 5.001, 0.1), 3)
t0 = time.time()
print(f"t_match={t_match}: scanning |E(kappa)|")
vals = []
for k in ks:
    v = EV.E(k, t_match=t_match); vals.append((k, abs(v), v.real))
    print(f"  k={k:5.2f}  |E|={abs(v):.4e}", flush=True)
print("\n|E| local minima (candidate eigenvalues):")
for i in range(1, len(vals)-1):
    if vals[i][1] < vals[i-1][1] and vals[i][1] < vals[i+1][1]:
        print(f"   kappa~{vals[i][0]:.2f}  |E|={vals[i][1]:.3e}   (beta=1/k={1/vals[i][0]:.4f})")
print(f"# {time.time()-t0:.1f}s")
