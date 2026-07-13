# nr_shoot.py — fluid-beta eigenvalue by HKA's ACTUAL method (rflanl.tex sec V, two-point BVP):
#   * At the sonic point the analytic solution is UNIQUE (1-param kappa) after imposing analyticity
#     (a_0 in ker R), the identity (eq:alg-PP), and the GAUGE Nbar_p(0)=0. [3 conds on the 4-vector a_0]
#   * Integrate that solution SONIC -> CENTER. The only growing center mode is (0,0,0,1)e^{-2x} (pure V_p;
#     the bounded modes and the identity-excluded e^{-x} mode all have V_p=0). Regularity (ii) =>
#     eliminate it => V_p(x_center) = 0. The eigenvalue is the SIGN CHANGE of V_p(x_c). kappa=2.81055255.
import numpy as np, sys, time
import hka_pert_core as PC, hka_pert_hka99 as H99
PC.Lnum = H99.Lnum
import hka_beta4 as B, hka_beta_solve as S
import nr_laurent as NL
from scipy.integrate import solve_ivp

B.bg(); B.bg_path(); XS = B.bg()['xs']

def psi_sonic(k, t0=-0.03, order=18):
    """Unique analytic-at-sonic solution, SIGN-CONSISTENT: identity=0, gauge Nbar_p(0)=0, and the HKA
    normalization Abar_p(0)=1 (3 conditions on the 3-D ker R -> unique, no SVD sign ambiguity)."""
    modes, R, Ls = NL.analytic_modes(k, order)
    nm = len(modes)
    a0 = np.array([modes[i][0] for i in range(nm)]).T          # 4 x nm  (ker R basis, t=0)
    idc = np.array(S.identity_coeffs(k, XS), complex)          # identity: idc . Psi = 0
    M = np.array([[complex(np.dot(idc, a0[:, i])) for i in range(nm)],   # identity = 0
                  [a0[1, i] for i in range(nm)],                        # gauge Nbar_p(0) = 0
                  [a0[0, i] for i in range(nm)]], complex)              # Abar_p(0) = 1
    c = np.linalg.solve(M, np.array([0, 0, 1], complex))
    val = np.zeros(4, complex)
    for i in range(nm):
        val += c[i]*sum(modes[i][n]*t0**n for n in range(len(modes[i])))
    return val

def _grow_leftvec(k, xc):
    L = H99.Lnum(B.bg_state(xc), complex(k))
    ev, Vr = np.linalg.eig(L.conj().T)
    j = int(np.argmin(ev.conj().real))     # left eigvec for the e^{-2x} (most negative) mode
    return Vr[:, j]

def E(k, xc=-13.0, t0=-0.03):
    """Sign-consistent growing-mode coefficient at the center (=0 at an eigenvalue). Sign change = kappa*."""
    psi = psi_sonic(k, t0)
    def rhs(x, y): return H99.Lnum(B.bg_state(x), complex(k)).dot(y)
    sol = solve_ivp(rhs, [XS+t0, xc], psi.astype(complex), method='DOP853', rtol=1e-10, atol=1e-13)
    yc = sol.y[:, -1]
    lv = _grow_leftvec(k, xc)
    coeff = np.vdot(lv.conj(), yc)          # projection onto the growing mode (raw, sign-consistent)
    return coeff, np.linalg.norm(yc)

def refine(ka, kb, xc=-13.0):
    fa = E(ka, xc)[0].real
    for _ in range(60):
        km = 0.5*(ka+kb); fm = E(km, xc)[0].real
        if fa*fm <= 0: kb = km
        else: ka, fa = km, fm
        if kb-ka < 1e-9: break
    return 0.5*(ka+kb)

if __name__ == "__main__":
    t0 = time.time()
    if len(sys.argv) > 2 and sys.argv[1] == "refine":
        ks = refine(float(sys.argv[2]), float(sys.argv[3]))
        print(f"kappa = {ks:.8f}   beta = 1/kappa = {1/ks:.8f}   (ref 0.35580192, err {abs(1/ks-0.35580192):.2e})")
    else:
        xc = float(sys.argv[1]) if len(sys.argv) > 1 else -13.0
        print(f"HKA shoot sonic->center (x_c={xc}): V_p(x_c) sign change = eigenvalue")
        print(f"{'kappa':>7} {'Re V_p(xc)/|y|':>16} {'|y(xc)|':>11}")
        prev = None
        for k in np.round(np.arange(0.5, 4.01, 0.1), 3):
            e, gr = E(k, xc); sc = ""
            if prev is not None and prev*e.real < 0: sc = "  <== SIGN CHANGE"
            print(f"{k:7.2f} {e.real:16.4e} {gr:11.2e}{sc}", flush=True)
            prev = e.real
        print(f"# {time.time()-t0:.1f}s")
