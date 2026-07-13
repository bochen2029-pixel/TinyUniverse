# hka_beta_validate.py — refine the relevant fluid-beta eigenvalue and run the honesty gates
# on the AUTHORITATIVE operator (hka_pert_hka99 via hka_beta_bvp collocation).
#   G-ANCHOR : |beta - 0.35580192| < 4e-3
#   G-CONVERGE : kappa stationary under mesh refinement (|k(N) - k(2N)| < 1e-4)
#   G-UNIQUE : exactly one relevant (Re kappa > 0) mode  (checked by the caller's sweep)
import numpy as np, sys, warnings
warnings.filterwarnings('ignore')
from scipy.integrate import solve_bvp
import hka_beta_bvp as M

BETA_REF = 0.35580192

def solve_at(kguess, npts, tol=1e-9):
    xg = np.linspace(M.XC, M.XM, npts)
    yg = np.zeros((4, xg.size))
    yg[0] = 0.1*np.exp(2*(xg - M.XM)); yg[2] = 1.0; yg[3] = 0.3*np.exp(xg - M.XM)
    s = solve_bvp(M.rhs, M.bc, xg, yg, p=[kguess], tol=tol, max_nodes=300000)
    return (s.p[0] if s.success else None), s

if __name__ == "__main__":
    kg = float(sys.argv[1]) if len(sys.argv) > 1 else 2.81
    # refine at two mesh densities for G-CONVERGE
    k1, s1 = solve_at(kg, 400)
    k2, s2 = solve_at(kg, 800)
    if k1 is None or k2 is None:
        print("solve failed"); sys.exit(2)
    beta = 1.0/k2
    conv = abs(k1 - k2)
    print(f"  kappa (N=400)  = {k1:.8f}")
    print(f"  kappa (N=800)  = {k2:.8f}   <-- reported")
    print(f"  beta = 1/kappa = {beta:.8f}   (ref 0.35580192)")
    print(f"  G-ANCHOR   |beta-ref|   = {abs(beta-BETA_REF):.2e}   {'PASS' if abs(beta-BETA_REF)<4e-3 else 'FAIL'}")
    print(f"  G-CONVERGE |k400-k800|  = {conv:.2e}   {'PASS' if conv<1e-4 else 'FAIL'}")
    print(f"  rel err vs published beta = {abs(beta-BETA_REF)/BETA_REF:.2e}")
