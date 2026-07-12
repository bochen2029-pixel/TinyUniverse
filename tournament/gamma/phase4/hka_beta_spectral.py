# hka_beta_spectral.py — Stage-B eigenvalue kappa via GLOBAL CHEBYSHEV COLLOCATION on the CORRECTED
# perturbation operator L(x;kappa) = P(x) + kappa Q(x)  (hka_pert_derive.py -> hka_pert_L.pkl ->
# hka_pert_core.Lnum), which passes the gauge-mode exactness gate (hka_pert.py, residuals ~1e-9).
#
# METHOD. Collocate  Psi'(x) = L(x;kappa) Psi(x)  on x in [x_c, x_r=x_s-delta] on a Chebyshev-Gauss-
# Lobatto grid. The discrete operator is (Dhat - Lhat(kappa)) with Dhat = spectral d/dx (block-Kron I4)
# and Lhat block-diagonal (4x4 per node). Replace 4 rows with the boundary conditions:
#   CENTER  (x_c): kill the two GROWING modes (as x->-inf a mode ~e^{lambda x} blows up for lambda<0;
#     eig(L_center)={-2,-1,0,0}). Rows = the two LEFT eigenvectors of L(x_c;kappa) for lambda in {-1,-2}.
#     (2 conditions -> Psi(x_c) in the lambda=0 plane, i.e. bounded + Abar_p=0.)
#   SONIC  (x_r): the sonic point is a regular singular pt, L ~ R/(x-x_s), R = R_P + kappa R_Q (AFFINE
#     in kappa) with eig(R) = {0,0,0,1-2kappa}. Impose (i) ANALYTICITY = kill the non-analytic
#     t^{1-2kappa} component: w(kappa) . Psi(x_r) = 0, w = LEFT eigenvector of R for eigenvalue 1-2kappa;
#     (ii) HKA gauge fix N̄_p(0) = 0: e_N . Psi(x_r) = 0.  (2 conditions.)
# Total 4 BC + collocation = square M(kappa). A nontrivial solution exists <=> det M(kappa) = 0
# <=> sigma_min(M(kappa)) = 0. Because L is AFFINE in kappa and the center-kill rows are affine, the ONLY
# nonlinearity in kappa is the sonic left-eigenvector w(kappa); we handle it exactly by a robust
# sigma_min(kappa) SWEEP + Brent refine (solver face `find_eigs`). No integration => IMMUNE to the
# exponential mode amplification (the t^{1-2kappa} stiffness) that defeats naive two-sided shooting.
#
# CROSS-CHECK. A linear generalized eigenproblem  (Dhat - Phat) Psi = kappa Qhat Psi  with the sonic
# analyticity row FROZEN at a trial kappa (self-consistent) returns ALL kappa at once (`gep_eigs`); the
# gauge modes MUST appear at kappa ~ 0.35699 (sonic gauge) and kappa = 1 (origin gauge) [HKA fn15] as a
# correctness check on the discretization, and are DISCARDED.
#
# HONESTY: no target anywhere. We report the measured kappa with convergence tables vs (N, x_c, delta),
# and beta = 1/Re(kappa). Eq #s per HKA_beta_equations.md.

import numpy as np, math
import scipy.linalg as sla
import hka_pert_core as PC
import hka_ec as E
import hka_beta4 as B

# ---------------------------------------------------------------------------
# Chebyshev-Gauss-Lobatto differentiation (Trefethen, cheb.m)
# ---------------------------------------------------------------------------
def cheb(N):
    if N == 0:
        return np.array([[0.0]]), np.array([1.0])
    x = np.cos(np.pi * np.arange(N + 1) / N)
    c = np.hstack([2., np.ones(N - 1), 2.]) * (-1) ** np.arange(N + 1)
    X = np.tile(x, (N + 1, 1)).T
    dX = X - X.T
    D = np.outer(c, 1. / c) / (dX + np.eye(N + 1))
    D = D - np.diag(D.sum(axis=1))
    return D, x

# ---------------------------------------------------------------------------
# background + operator sampling (EC background via hka_beta4.bg_state = solve_ivp dense output)
# ---------------------------------------------------------------------------
def ensure_bg():
    B.bg(); B.bg_path()
    return B.bg()['xs']

def _PQ(x):
    st = B.bg_state(x)
    return PC.Pmat(st), PC.Qmat(st)

def _map(N, xc, xr):
    D, t = cheb(N)
    xph = 0.5 * ((xr - xc) * t + (xr + xc))     # t=+1 -> xr (sonic side), t=-1 -> xc (center side)
    Dx = D * (2.0 / (xr - xc))
    return xph, Dx

# ---------------------------------------------------------------------------
# residue R = R_P + kappa R_Q  (affine), by Richardson extrapolation of t*L along the EC background
# ---------------------------------------------------------------------------
_RCACHE = {}
def residue_PQ(xs, dnear=5e-5):
    key = round(dnear, 12)
    if key in _RCACHE:
        return _RCACHE[key]
    def tP(t):
        st = B.bg_state(xs + t); return t * PC.Pmat(st)
    def tQ(t):
        st = B.bg_state(xs + t); return t * PC.Qmat(st)
    t1, t2 = -2 * dnear, -dnear
    RP = tP(t2) + (tP(t2) - tP(t1)) * ((0 - t2) / (t2 - t1))
    RQ = tQ(t2) + (tQ(t2) - tQ(t1)) * ((0 - t2) / (t2 - t1))
    _RCACHE[key] = (RP, RQ)
    return RP, RQ

# ---------------------------------------------------------------------------
# BC functionals
# ---------------------------------------------------------------------------
def center_kill_rows(xc, kappa):
    """Two LEFT eigenvectors of L(xc;kappa) for the growing modes lambda in {-1,-2}."""
    P, Q = _PQ(xc)
    L = P + kappa * Q
    evL, VL = np.linalg.eig(L.T)
    rows = []
    for target in (-1.0, -2.0):
        i = int(np.argmin(np.abs(evL - target)))
        rows.append(np.conj(VL[:, i]))
    return np.array(rows, complex)

def sonic_analytic_row(xs, kappa):
    """LEFT eigenvector w of R = R_P + kappa R_Q for eigenvalue 1-2kappa (the non-analytic exponent);
    w . Psi = 0 kills the non-analytic t^{1-2kappa} component."""
    RP, RQ = residue_PQ(xs)
    R = RP + kappa * RQ
    evL, VL = np.linalg.eig(R.T)
    i = int(np.argmin(np.abs(evL - (1.0 - 2.0 * kappa))))
    return np.conj(VL[:, i])

_ENB = np.array([0., 1., 0., 0.], complex)   # gauge fix N̄_p = 0

# ---------------------------------------------------------------------------
# PRIMARY: sigma_min(M(kappa)) sweep  (M square, all 4 BCs, sonic row uses the true kappa-dep w)
# ---------------------------------------------------------------------------
def build_M(kappa, N, xc, delta):
    kappa = complex(kappa)
    xs = B.bg()['xs']
    xr = xs - delta
    xph, Dx = _map(N, xc, xr)
    n = N + 1
    I4 = np.eye(4)
    M = np.kron(Dx, I4).astype(complex)
    for i in range(n):
        P, Q = _PQ(xph[i])
        M[4 * i:4 * i + 4, 4 * i:4 * i + 4] -= (P + kappa * Q)
    i_son = 0            # node at xr
    i_cen = n - 1        # node at xc
    cr = center_kill_rows(xc, kappa)
    for m in range(2):
        row = np.zeros(4 * n, complex); row[4 * i_cen:4 * i_cen + 4] = cr[m]
        M[4 * i_cen + m, :] = row
    wa = sonic_analytic_row(xs, kappa)
    row = np.zeros(4 * n, complex); row[4 * i_son:4 * i_son + 4] = wa
    M[4 * i_son + 0, :] = row
    row = np.zeros(4 * n, complex); row[4 * i_son:4 * i_son + 4] = _ENB
    M[4 * i_son + 1, :] = row
    return M

def sigma_min(kappa, N, xc, delta):
    M = build_M(kappa, N, xc, delta)
    nrm = np.linalg.norm(M, axis=1); nrm[nrm == 0] = 1.0
    s = sla.svd(M / nrm[:, None], compute_uv=False)
    return s[-1]

def find_eigs(N=48, xc=-14.0, delta=1e-3, kmin=0.2, kmax=5.0, nscan=200,
              sig_gate=0.05, refine=True, verbose=False):
    """Scan sigma_min(kappa); bracket local minima that dip below sig_gate; Brent-refine each."""
    from scipy.optimize import minimize_scalar
    ks = np.linspace(kmin, kmax, nscan)
    sig = np.array([sigma_min(k, N, xc, delta) for k in ks])
    if verbose:
        for k, s in zip(ks, sig):
            print(f"    {k:7.3f}  {s:.4e}")
    roots = []
    for i in range(1, nscan - 1):
        if sig[i] < sig[i - 1] and sig[i] < sig[i + 1] and sig[i] < sig_gate:
            a, b = ks[i - 1], ks[i + 1]
            if refine:
                r = minimize_scalar(lambda k: sigma_min(k, N, xc, delta),
                                    bracket=(a, ks[i], b), method='brent', options=dict(xtol=1e-10))
                roots.append((float(r.x), float(r.fun)))
            else:
                roots.append((float(ks[i]), float(sig[i])))
    return sorted(roots), ks, sig

def label(k):
    if abs(k - 1.0) < 0.04: return "gauge (kappa=1, origin)"
    if abs(k - 0.35699) < 0.03: return "gauge (kappa~0.357, sonic)"
    return "PHYSICAL?"

# ---------------------------------------------------------------------------
# CROSS-CHECK: linear GEP with frozen sonic-analyticity row -> all kappa at once
# ---------------------------------------------------------------------------
def gep_eigs(N=48, xc=-14.0, delta=1e-3, w_kappa=2.8, kmin=0.2, kmax=6.0, imtol=1e-4):
    xs = B.bg()['xs']; xr = xs - delta
    xph, Dx = _map(N, xc, xr)
    n = N + 1; I4 = np.eye(4)
    A0 = np.kron(Dx, I4).astype(complex)
    Bm = np.zeros((4 * n, 4 * n), complex)
    for i in range(n):
        P, Q = _PQ(xph[i])
        A0[4 * i:4 * i + 4, 4 * i:4 * i + 4] -= P
        Bm[4 * i:4 * i + 4, 4 * i:4 * i + 4] = Q
    i_son = 0; i_cen = n - 1
    cr = center_kill_rows(xc, w_kappa)
    for m in range(2):
        A0[4 * i_cen + m, :] = 0.0; Bm[4 * i_cen + m, :] = 0.0
        A0[4 * i_cen + m, 4 * i_cen:4 * i_cen + 4] = cr[m]
    wa = sonic_analytic_row(xs, w_kappa)
    A0[4 * i_son + 0, :] = 0.0; Bm[4 * i_son + 0, :] = 0.0
    A0[4 * i_son + 0, 4 * i_son:4 * i_son + 4] = wa
    A0[4 * i_son + 1, :] = 0.0; Bm[4 * i_son + 1, :] = 0.0
    A0[4 * i_son + 1, 4 * i_son:4 * i_son + 4] = _ENB
    w = sla.eigvals(A0, Bm)
    w = w[np.isfinite(w)]
    real = w[np.abs(w.imag) < imtol * np.maximum(np.abs(w.real), 1.0)].real
    return np.sort(real[(real > kmin) & (real < kmax)])

# ---------------------------------------------------------------------------
if __name__ == "__main__":
    import time, sys
    t0 = time.time()
    xs = ensure_bg()
    print(f"# sonic x_s = {xs:.6f}")
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 40
    xc = float(sys.argv[2]) if len(sys.argv) > 2 else -14.0
    delta = float(sys.argv[3]) if len(sys.argv) > 3 else 1e-3
    print(f"# Chebyshev collocation, corrected L=P+kappa Q.  N={N}  xc={xc}  delta={delta}\n")
    roots, ks, sig = find_eigs(N=N, xc=xc, delta=delta)
    print("sigma_min(kappa) local minima (candidate eigenvalues):")
    for k, s in roots:
        print(f"   kappa = {k:.7f}   sigma_min = {s:.3e}   [{label(k)}]")
    if not roots:
        print("   (none below gate; global min at kappa=%.4f sigma=%.3e)" % (ks[np.argmin(sig)], sig.min()))
    print(f"# {time.time()-t0:.1f}s")
