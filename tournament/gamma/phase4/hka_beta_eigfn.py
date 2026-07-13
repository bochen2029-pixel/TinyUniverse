# hka_beta_eigfn.py — extract the eigenfunction for a given kappa and print its profile, to validate
# the recovered mode against HKA99 Fig 2 (curves: Abar_p, Nbar_p/Nbar_ss, obar_p/obar_ss, V_p) and to
# confirm it is a SMOOTH, node-structured physical mode (not a spurious BC artifact).
import numpy as np, sys, warnings
warnings.filterwarnings('ignore')
from scipy.integrate import solve_bvp
import hka_beta_bvp as M

def eigfn(kguess, npts=500, tol=1e-9):
    xg = np.linspace(M.XC, M.XM, npts); yg = np.zeros((4, xg.size))
    yg[0] = 0.1*np.exp(2*(xg - M.XM)); yg[2] = 1.0; yg[3] = 0.3*np.exp(xg - M.XM)
    s = solve_bvp(M.rhs, M.bc, xg, yg, p=[kguess], tol=tol, max_nodes=300000)
    return s

if __name__ == "__main__":
    kg = float(sys.argv[1]) if len(sys.argv) > 1 else 2.81
    s = eigfn(kg)
    if not s.success:
        print(f"solve failed (status {s.status})"); sys.exit(2)
    k = s.p[0]
    print(f"kappa={k:.8f}  beta=1/k={1.0/k:.8f}")
    print(f"{'x':>8} {'Abar_p':>11} {'Nbar_p':>11} {'obar_p':>11} {'V_p':>11}")
    for x in np.linspace(M.XC, M.XM, 11):
        y = s.sol(x)
        print(f"{x:8.3f} {y[0]:11.4f} {y[1]:11.4f} {y[2]:11.4f} {y[3]:11.4f}")
    # count sign changes (nodes) of each component -> fundamental mode is smoothest
    xs = np.linspace(M.XC, M.XM, 400); Y = s.sol(xs)
    nodes = [int(np.sum(np.diff(np.sign(Y[i])) != 0)) for i in range(4)]
    print(f"nodes (Abar,Nbar,obar,V) = {nodes}")
