# nr_match.py — fluid-beta eigenvalue by the two-sided match det[ v_center | s1 | s2 | s3 ] on the
# VERIFIED operator with the EXACT sonic Frobenius modes (nr_laurent). Zero <=> the center-regular
# mode lies in the 3-D sonic-analytic span <=> eigenvalue. Physical: kappa=2.81055255. Gauge root
# at kappa~1 (this background's gauge, N_bar'(x_s)~1) is DISCARDED (HKA fn15).
import numpy as np, sys
import hka_pert_core as PC, hka_pert_hka99 as H99
PC.Lnum = H99.Lnum                    # route hka_beta4's center machinery through the verified operator
import hka_beta4 as B
import nr_laurent as NL

B.bg(); B.bg_path(); XS = B.bg()['xs']

def Delta(kappa, xc=-13.0, tm=-0.03, order=18, h=1e-3):
    xm = XS + tm
    vc = B.v_center(kappa, xc)
    c = B.integ(vc, kappa, xc, xm, h=h)                 # center-regular integrated to x_m (unit)
    if c is None: return None
    modes, R, Ls = NL.analytic_modes(kappa, order)
    if len(modes) < 3: return None
    scols = [sum(a[n]*tm**n for n in range(len(a))) for a in modes]
    Q, _ = np.linalg.qr(np.column_stack(scols))         # orthonormal basis of the analytic 3-space
    Q = Q[:, :3]
    rej = np.linalg.norm(c - Q.dot(Q.conj().T.dot(c)))  # rejection: 0 <=> c in span (scale-free)
    det = np.linalg.det(np.column_stack([c, Q]))
    return det, rej

def refine(ka, kb, **kw):
    fa = Delta(ka, **kw)[0].real
    for _ in range(80):
        km = 0.5*(ka+kb); fm = Delta(km, **kw)[0].real
        if fa*fm <= 0: kb = km
        else: ka, fa = km, fm
        if kb-ka < 1e-11: break
    return 0.5*(ka+kb)

if __name__ == "__main__":
    import time
    t0 = time.time()
    if len(sys.argv) > 2 and sys.argv[1] == "refine":
        ks = refine(float(sys.argv[2]), float(sys.argv[3]))
        print(f"kappa = {ks:.8f}   beta = 1/kappa = {1/ks:.8f}   (ref 0.35580192, err {abs(1/ks-0.35580192):.2e})")
    else:
        print(f"match det[v_center|s1|s2|s3], sonic x_s={XS:.5f}, x_m=x_s{-0.03}:")
        print(f"{'kappa':>7} {'Re det':>13} {'|det|':>11} {'rejection':>11}")
        rows = []
        for kap in np.arange(0.3, 6.01, 0.1):
            r = Delta(kap)
            if r is None: print(f"{kap:7.2f}  (fail)"); continue
            det, rej = r; rows.append((kap, det.real, rej))
            print(f"{kap:7.2f} {det.real:13.4e} {abs(det):11.3e} {rej:11.3e}")
        print("\nRe(det) sign changes:")
        for i in range(1, len(rows)):
            if rows[i-1][1]*rows[i][1] < 0:
                print(f"   {rows[i-1][0]:.2f} -> {rows[i][0]:.2f}")
        print("rejection local minima:")
        for i in range(1, len(rows)-1):
            if rows[i][2] < rows[i-1][2] and rows[i][2] < rows[i+1][2]:
                print(f"   kappa~{rows[i][0]:.2f}  rej={rows[i][2]:.3e}")
        print(f"# {time.time()-t0:.1f}s")
