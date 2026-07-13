# hka_beta_bvp.py — fluid-beta eigenvalue via COLLOCATION (scipy.solve_bvp), on the verified
# authoritative operator L = M_x^{-1}(Gmat - kappa M_s) (hka_pert_hka99). Real system (the relevant
# mode is real: kappa=2.81055255, real eigenfunction, HKA99 Fig 2). No shooting -> no mode collapse.
#
# 4-D ODE  y' = L(x;kappa) y,  y=(Abar_p,Nbar_p,obar_p,V_p),  kappa = unknown eigenvalue parameter.
# BCs (5 = n+1 for one parameter):
#   center x_c:  Abar_p=0, V_p=0        (regular: bounded modes are (0,1,0,0),(0,0,1,0))
#   sonic  x_m:  Nbar_p=0 (gauge),  obar_p=1 (normalization),
#                n.(Gmat - kappa M_s) y = 0   (sonic solvability: RHS in range of the singular M_x)
import numpy as np, sys
from scipy.integrate import solve_bvp
import hka_pert_core as PC, hka_pert_hka99 as H99
PC.Lnum = H99.Lnum
import hka_beta4 as B
B.bg(); B.bg_path(); XS = B.bg()['xs']

def mats(x, kappa):
    c = H99.coeffs(B.bg_state(x))
    Mx = np.array([[1,0,0,0],[0,1,0,0],[0,0,c['Ax'],c['Bx']],[0,0,c['Cx'],c['Dx']]], float)
    Gm = np.array([[c['G1'],0,c['G3'],c['G4']],[c['H1'],0,c['H3'],0],
                   [c['E1'],c['E2'],c['E3'],c['E4']],[c['F1'],c['F2'],c['F3'],c['F4']]], float)
    Ms = np.array([[0,0,0,0],[0,0,0,0],[0,0,c['As'],c['Bs']],[0,0,c['Cs'],c['Ds']]], float)
    return Mx, Gm, Ms

def Lreal(x, kappa):
    Mx, Gm, Ms = mats(x, kappa)
    return np.linalg.solve(Mx, Gm - kappa*Ms)

def rhs(x, y, p):
    k = p[0]
    out = np.empty_like(y)
    for j in range(x.size):
        out[:, j] = Lreal(x[j], k).dot(y[:, j])
    return out

def bc(ya, yb, p):
    k = p[0]
    # sonic analyticity: RHS must lie in range of the SINGULAR fluid block [[Ax,Bx],[Cx,Dx]] at x_s
    # (det = AxDx - BxCx = 0 there). Its exact left-null is m = (Cx, -Ax); embed n = (0,0,Cx,-Ax).
    c = H99.coeffs(B.bg_state(XS - 1e-7))
    _, Gm, Ms = mats(XS - 1e-7, k)
    n = np.array([0.0, 0.0, c['Cx'], -c['Ax']])
    solv = n.dot((Gm - k*Ms).dot(yb))
    return np.array([ya[0], ya[3], yb[1], yb[2] - 1.0, solv])

XC, XM = -13.0, XS - 0.02

if __name__ == "__main__":
    kguess = float(sys.argv[1]) if len(sys.argv) > 1 else 2.81
    xg = np.linspace(XC, XM, 400)
    # initial guess: smooth, regular at center, obar~1, small Abar, some V
    yg = np.zeros((4, xg.size))
    yg[0] = 0.1*np.exp(2*(xg - XM))      # Abar small, ->0 at center
    yg[1] = 0.0                          # Nbar
    yg[2] = 1.0                          # obar ~ 1
    yg[3] = 0.3*np.exp((xg - XM))        # V ->0 at center
    sol = solve_bvp(rhs, bc, xg, yg, p=[kguess], tol=1e-8, max_nodes=200000, verbose=2)
    k = sol.p[0]
    print(f"\nsolve_bvp: success={sol.success}  status={sol.status}")
    print(f"  kappa = {k:.8f}   beta = 1/kappa = {1.0/k:.8f}   (target kappa=2.81055255, beta=0.35580192)")
    print(f"  rel err vs beta: {abs(1.0/k - 0.35580192)/0.35580192:.2e}")
