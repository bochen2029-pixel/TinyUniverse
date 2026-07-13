# hka_beta_match.py — fluid-beta eigenvalue by the reduced-3D forward-shoot 2x2 match determinant,
# on the verified authoritative operator (hka_pert_hka99). This is the robust formulation:
#  * REDUCED 3-D system (Nbar,obar,V), Abar eliminated via the algebraic identity (eq:alg-PP) -> the
#    identity is built in (no drift; the 4-row ODE does NOT preserve it).
#  * FORWARD integration center -> sonic: the growing-at-center e^{-2x} mode DECAYS forward (stable),
#    avoiding the violent +7 eigenvalue near the sonic point that blows up backward integration.
#  * At the sonic side, the 2 center-regular modes must satisfy gauge Nbar(x_m)=0 AND sonic analyticity
#    (solvability) -> a 2x2 determinant Delta(kappa); its zeros (sign changes) are the eigenvalues.
# Physical fluid-beta mode: kappa = 2.81055255, beta = 0.35580192 (HKA99 Table I, unique relevant mode).
import numpy as np, sys
from scipy.integrate import solve_ivp
import hka_pert_core as PC, hka_pert_hka99 as H99
PC.Lnum = H99.Lnum
import hka_beta4 as B, hka_beta_solve as S
B.bg(); B.bg_path(); XS = B.bg()['xs']

def L3(x, k):
    fld = B.bg_state(x); L4 = H99.Lnum(fld, k); idc = S.identity_coeffs(k, x)
    Mrow = -np.array([idc[1], idc[2], idc[3]], complex)/idc[0]
    T = np.array([Mrow, [1, 0, 0], [0, 1, 0], [0, 0, 1]], complex)
    return L4[1:4, :].dot(T)

def lift(y3, k, x):
    idc = S.identity_coeffs(k, x)
    Ab = -(idc[1]*y3[0] + idc[2]*y3[1] + idc[3]*y3[2])/idc[0]
    return np.array([Ab, y3[0], y3[1], y3[2]], complex)

def integ(v3, k, xc, xm):
    s = solve_ivp(lambda x, y: L3(x, k).dot(y), [xc, xm], np.array(v3, complex),
                  method='DOP853', rtol=1e-10, atol=1e-13)
    return s.y[:, -1]

def sonic_mats(k):
    c = H99.coeffs(B.bg_state(XS - 1e-7))
    Gm = np.array([[c['G1'],0,c['G3'],c['G4']],[c['H1'],0,c['H3'],0],
                   [c['E1'],c['E2'],c['E3'],c['E4']],[c['F1'],c['F2'],c['F3'],c['F4']]], complex)
    Ms = np.array([[0,0,0,0],[0,0,0,0],[0,0,c['As'],c['Bs']],[0,0,c['Cs'],c['Ds']]], complex)
    return Gm, Ms, c

def Delta(k, xc=-13.0, xm=None):
    if xm is None: xm = XS - 0.03
    L1 = integ([1, 0, 0], k, xc, xm); L2 = integ([0, 1, 0], k, xc, xm)
    Gm, Ms, c = sonic_mats(k)
    def conds(y3):
        y4 = lift(y3, k, xm); r = (Gm - k*Ms).dot(y4)
        return y3[0], c['Cx']*r[2] - c['Ax']*r[3]         # (gauge Nbar, solvability)
    g1, s1 = conds(L1); g2, s2 = conds(L2)
    return (g1*s2 - g2*s1)

def refine(ka, kb, xc=-13.0, xm=None):
    da = Delta(ka, xc, xm).real
    for _ in range(60):
        km = 0.5*(ka + kb); dm = Delta(km, xc, xm).real
        if da*dm <= 0: kb = km
        else: ka, da = km, dm
        if kb - ka < 1e-10: break
    return 0.5*(ka + kb)

if __name__ == "__main__":
    if len(sys.argv) > 2 and sys.argv[1] == "refine":
        ka, kb = float(sys.argv[2]), float(sys.argv[3])
        ks = refine(ka, kb)
        print(f"kappa = {ks:.8f}   beta = 1/kappa = {1.0/ks:.8f}   (ref 0.35580192, err {abs(1/ks-0.35580192):.2e})", flush=True)
    else:
        lo, hi, step = 0.3, 4.0, 0.1
        print(f"reduced-3D 2x2 match determinant (sign change = eigenvalue):", flush=True)
        prev = None; prevk = None
        k = lo
        while k <= hi + 1e-9:
            d = Delta(k).real
            sc = ""
            if prev is not None and prev*d < 0:
                sc = f"  <== SIGN CHANGE in ({prevk:.2f},{k:.2f})"
            print(f"  kappa={k:.2f}  Delta={d:+.4e}{sc}", flush=True)
            prev, prevk = d, k; k += step
