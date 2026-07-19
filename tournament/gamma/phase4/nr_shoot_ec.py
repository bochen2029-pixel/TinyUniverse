# nr_shoot_ec.py — HKA's shoot (rflanl.tex sec V, the "first method") on the TRUE EC
# background (nr_ec2, D-032): redundant recovery of kappa_0 = 2.8105525488 -> beta.
#
# The Lyapunov one-step-map spectrum already lands kappa = 2.8105526(3); this is the
# SECOND, independent method (and the high-precision one — HKA: "easier ... to obtain
# accurate values of kappa"). TWO analytic controls in this (sonic) gauge:
#   * the sonic-gauge gauge mode at kappa = -Nbar'_ss(sonic) = +0.355699
#     (the rflanl fn value — a dropped-digit typo prints it as 0.35699), and
#   * the relevant mode at kappa = 2.8105525488.
# Construction identical to nr_shoot.py (unique analytic-at-sonic solution: a0 in ker R
# + identity eq:alg-PP + gauge Nbar_p(0)=0, norm Abar_p(0)=1; integrate sonic->center;
# kill the (0,0,0,1)e^{-2x} expanding mode), with all background inputs = the true EC.
import numpy as np, math, sys, time
import nr_ec2 as EC
import nr_laurent as NL
import hka_pert_hka99 as H99
from scipy.integrate import solve_ivp

G = 4.0/3.0
ec = EC.build_ec(verbose=False)
XS = ec['xs']; BG = ec['bg_state']; V0 = ec['V0']
A0, N0, om0 = EC.sonic_values(V0)

_BRS = None
def branch_series(order=18):
    """EC-branch sonic Taylor series: the bg_series_at branch matching the dense solution."""
    global _BRS
    if _BRS is not None: return _BRS
    brs = EC.bg_series_at(V0, order)
    best, bd = None, 1e99
    for br in brs:
        d = 0.0
        for t in (-0.01, -0.03):
            val = EC.seed_state(br, t)                 # (N, om, V)
            ref = ec['sol'].sol(XS + t)
            d += float(np.abs(val - ref).max())
        if d < bd: bd, best = d, br
    _BRS = (best, bd)
    return _BRS

def identity_coeffs(kappa):
    """eq:alg-PP at the EC sonic point: (kappa-A0) Abar + cN Nbar + com obar + cV V = 0."""
    oV2 = 1 - V0*V0
    cN = 2*G*N0*V0*om0/oV2
    com = 2*om0*(1 + (G-1)*V0*V0 + G*N0*V0)/oV2
    cV = 2*G*om0*(N0*(1 + V0*V0) + 2*V0)/oV2**2
    return np.array([complex(kappa) - A0, cN, com, cV], complex)

def psi_sonic(k, t0=-0.02, order=18):
    br, _ = branch_series(order)
    modes, R, Ls = NL.analytic_modes(k, order, bgser=br)
    nm = len(modes)
    a0 = np.array([modes[i][0] for i in range(nm)]).T          # 4 x nm
    idc = identity_coeffs(k)
    M = np.array([[complex(np.dot(idc, a0[:, i])) for i in range(nm)],
                  [a0[1, i] for i in range(nm)],
                  [a0[0, i] for i in range(nm)]], complex)
    c = np.linalg.solve(M, np.array([0, 0, 1], complex))
    val = np.zeros(4, complex)
    for i in range(nm):
        val += c[i]*sum(modes[i][n]*t0**n for n in range(len(modes[i])))
    return val

def _grow_leftvec(k, xc):
    L = H99.Lnum(BG(xc), complex(k))
    ev, Vr = np.linalg.eig(L.conj().T)
    j = int(np.argmin(ev.conj().real))
    return Vr[:, j]

def E(k, xc=-12.0, t0=-0.02):
    psi = psi_sonic(k, t0)
    def rhs(x, y): return H99.Lnum(BG(x), complex(k)).dot(y)
    sol = solve_ivp(rhs, [XS + t0, xc], psi.astype(complex), method='DOP853',
                    rtol=1e-11, atol=1e-13)
    yc = sol.y[:, -1]
    lv = _grow_leftvec(k, xc)
    return np.vdot(lv.conj(), yc), np.linalg.norm(yc)

def refine(ka, kb, xc=-12.0):
    fa = E(ka, xc)[0].real
    for _ in range(80):
        km = 0.5*(ka + kb); fm = E(km, xc)[0].real
        if fa*fm <= 0: kb = km
        else: ka, fa = km, fm
        if kb - ka < 1e-10: break
    return 0.5*(ka + kb)

if __name__ == "__main__":
    t0w = time.time()
    br, bd = branch_series()
    print(f"[series] EC-branch sonic Taylor vs dense solution: match {bd:.2e} (sum over t=-0.01,-0.03)")
    if len(sys.argv) > 2 and sys.argv[1] == "refine":
        ks = refine(float(sys.argv[2]), float(sys.argv[3]))
        print(f"kappa = {ks:.9f}   beta = 1/kappa = {ks and 1/ks:.9f}")
        print(f"  (relevant ref 2.8105525488 -> beta 0.35580192; sonic-gauge control ref {-(-2+A0-(2-G)*om0):.6f})")
    else:
        xc = float(sys.argv[1]) if len(sys.argv) > 1 else -12.0
        print(f"EC shoot sonic->center (x_c={xc}): growing-mode coefficient; zero = eigenvalue")
        print(f"{'kappa':>7} {'Re coeff':>13} {'|y(xc)|':>10}")
        prev = None
        for k in np.round(np.arange(0.15, 3.61, 0.05), 3):
            e, gr = E(k, xc); sc = ""
            if prev is not None and prev*e.real < 0: sc = "  <== SIGN CHANGE"
            print(f"{k:7.2f} {e.real:13.4e} {gr:10.2e}{sc}", flush=True)
            prev = e.real
        print(f"# {time.time()-t0w:.1f}s")
