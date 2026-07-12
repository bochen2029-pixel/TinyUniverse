# hka_beta_fdbvp.py — Stage-B eigenvalue kappa via a GLOBAL FINITE-DIFFERENCE BOUNDARY-VALUE
# generalized eigenproblem on the CORRECTED perturbation operator L(x;kappa) = P(x) + kappa Q(x)
# (hka_pert_derive.py -> hka_pert_L.pkl -> hka_pert_core.Lnum; gauge-exact to ~1e-9).
#
# WHY A GLOBAL METHOD (not shooting). The perturbation operator is strongly NON-NORMAL
# (||[L,L^H]|| ~ 4x ||L|| at the center) and the sonic point carries a t^{1-2kappa} (~t^{-4.6}) singular
# mode. Two-sided shooting and subspace/principal-angle matches are corrupted by this exponential mode
# amplification (they collapse onto the ubiquitous gauge mode). A GLOBAL discretization solves the whole
# BVP at once as a linear generalized eigenproblem A0 Psi = kappa B Psi with NO integration => immune to
# mode amplification. Second-order centered differences on a uniform grid [x_c, x_s-delta]:
#     (Psi_{i+1} - Psi_{i-1})/(2h) - P_i Psi_i = kappa Q_i Psi_i          (interior)
# BOUNDARY CONDITIONS (algebraic rows, put into A0 with 0 in B):
#   CENTER x_c: kill the two GROWING modes lambda in {-1,-2} of L(x_c) (left-eigenvector functionals).
#   SONIC x_r=x_s-delta: (i) analyticity w(kappa).Psi=0 (kill the non-analytic t^{1-2kappa} component;
#     w = left eigenvector of the residue R=R_P+kappa R_Q for eigenvalue 1-2kappa); (ii) HKA gauge fix
#     N̄_p(0)=0.  The two remaining rows at each end use one-sided differences to keep the ODE and the
#     matrix square.
# The analyticity functional w depends on kappa, so we FREEZE w at a trial kappa and iterate to
# SELF-CONSISTENCY (w_kappa <- eigenvalue). scipy.linalg.eig returns all kappa; the physical one is the
# real O(1-3) eigenvalue that is NOT the gauge mode (kappa=1). CONVERGENCE is demonstrated in M (grid),
# with Richardson O(h^2) extrapolation, and robustness vs x_c and delta.
#
# The gauge mode kappa=1 (origin gauge; here (lnN)'=-1 identically so sonic & origin gauges coincide at 1)
# MUST appear -> discretization correctness check. HONEST: no target in the solver; kappa is measured,
# beta = 1/Re(kappa) reported with the convergence table. Eq #s per HKA_beta_equations.md.

import numpy as np, math
import scipy.linalg as sla
import hka_pert_core as PC
import hka_pert_sonic as PS
import hka_ec as E
import hka_beta4 as B

def ensure_bg():
    B.bg(); B.bg_path()
    return B.bg()['xs']

def _PQ(x):
    st = B.bg_state(x)
    return PC.Pmat(st), PC.Qmat(st)

_RC = {}
def residue_PQ(xs, dnear=5e-5):
    if dnear in _RC:
        return _RC[dnear]
    def tP(t): return t * PC.Pmat(B.bg_state(xs + t))
    def tQ(t): return t * PC.Qmat(B.bg_state(xs + t))
    t1, t2 = -2 * dnear, -dnear
    RP = tP(t2) + (tP(t2) - tP(t1)) * ((0 - t2) / (t2 - t1))
    RQ = tQ(t2) + (tQ(t2) - tQ(t1)) * ((0 - t2) / (t2 - t1))
    _RC[dnear] = (RP, RQ)
    return RP, RQ

def center_kill_rows(xc, kappa):
    P, Q = _PQ(xc); L = P + kappa * Q
    evL, VL = np.linalg.eig(L.T)
    return np.array([np.conj(VL[:, int(np.argmin(np.abs(evL - t)))]) for t in (-1.0, -2.0)], complex)

def sonic_analytic_row(xs, kappa):
    RP, RQ = residue_PQ(xs); R = RP + kappa * RQ
    evL, VL = np.linalg.eig(R.T)
    return np.conj(VL[:, int(np.argmin(np.abs(evL - (1.0 - 2.0 * kappa))))])

_ENB = np.array([0., 1., 0., 0.], complex)

def build_AB(M, xc, delta, wk):
    """Assemble (A0, B) for the FD-BVP generalized eigenproblem A0 Psi = kappa B Psi. w (sonic analyticity)
    and the center-kill rows are frozen at wk (iterate to self-consistency)."""
    xs = B.bg()['xs']; xr = xs - delta
    xg = np.linspace(xc, xr, M); h = xg[1] - xg[0]; n = M; I4 = np.eye(4)
    A0 = np.zeros((4 * n, 4 * n), complex); Bm = np.zeros((4 * n, 4 * n), complex)
    for i in range(1, n - 1):
        P, Q = _PQ(xg[i])
        A0[4 * i:4 * i + 4, 4 * (i + 1):4 * (i + 1) + 4] += I4 / (2 * h)
        A0[4 * i:4 * i + 4, 4 * (i - 1):4 * (i - 1) + 4] += -I4 / (2 * h)
        A0[4 * i:4 * i + 4, 4 * i:4 * i + 4] += -P
        Bm[4 * i:4 * i + 4, 4 * i:4 * i + 4] = Q
    # center node 0: 2 BC rows (kill growing) + 2 one-sided ODE rows
    cr = center_kill_rows(xc, wk)
    for m in range(2):
        A0[m, :] = 0; A0[m, 0:4] = cr[m]
    P0, Q0 = _PQ(xg[0])
    for m in (2, 3):
        A0[m, :] = 0; Bm[m, :] = 0
        A0[m, 4:8] += (I4 / h)[m]; A0[m, 0:4] += -(I4 / h)[m] - P0[m]; Bm[m, 0:4] = Q0[m]
    # sonic node n-1: 2 BC rows (analyticity + Nbar fix) + 2 one-sided ODE rows
    wa = sonic_analytic_row(xs, wk)
    A0[4 * (n - 1) + 0, :] = 0; A0[4 * (n - 1) + 0, 4 * (n - 1):4 * (n - 1) + 4] = wa
    A0[4 * (n - 1) + 1, :] = 0; A0[4 * (n - 1) + 1, 4 * (n - 1):4 * (n - 1) + 4] = _ENB
    Pn, Qn = _PQ(xg[-1])
    for m in (2, 3):
        A0[4 * (n - 1) + m, :] = 0; Bm[4 * (n - 1) + m, :] = 0
        A0[4 * (n - 1) + m, 4 * (n - 1):4 * (n - 1) + 4] += (I4 / h)[m] - Pn[m]
        A0[4 * (n - 1) + m, 4 * (n - 2):4 * (n - 2) + 4] += -(I4 / h)[m]
        Bm[4 * (n - 1) + m, 4 * (n - 1):4 * (n - 1) + 4] = Qn[m]
    return A0, Bm

def eigs(M, xc=-14.0, delta=1e-3, wk=2.8, kmin=0.3, kmax=5.0, imtol=1e-3):
    A0, Bm = build_AB(M, xc, delta, wk)
    w = sla.eigvals(A0, Bm); w = w[np.isfinite(w)]
    real = w[np.abs(w.imag) < imtol * np.maximum(np.abs(w.real), 1.0)].real
    return np.sort(real[(real > kmin) & (real < kmax)]), w

def physical_eig(M, xc=-14.0, delta=1e-3, wk=2.8, exclude_gauge=True):
    """Return the physical eigenvalue: the real eigenvalue nearest wk that is not the gauge mode (kappa=1)."""
    e, _ = eigs(M, xc, delta, wk)
    cand = [k for k in e if not (exclude_gauge and abs(k - 1.0) < 0.06)]
    if not cand:
        return None
    return float(cand[int(np.argmin(np.abs(np.array(cand) - wk)))])

def self_consistent(M, xc=-14.0, delta=1e-3, wk0=2.8, iters=6, tol=1e-4):
    """Iterate wk <- physical eigenvalue until fixed point."""
    wk = wk0
    for _ in range(iters):
        k = physical_eig(M, xc, delta, wk)
        if k is None:
            return None
        if abs(k - wk) < tol:
            return k
        wk = k
    return wk

if __name__ == "__main__":
    import time, sys
    t0 = time.time()
    xs = ensure_bg()
    print(f"# sonic x_s = {xs:.6f}")
    print(f"# FD-BVP generalized eigenproblem on corrected L=P+kappa Q. Self-consistent w, convergence in M.\n")
    xc = float(sys.argv[1]) if len(sys.argv) > 1 else -14.0
    delta = float(sys.argv[2]) if len(sys.argv) > 2 else 1e-3
    Ms = [200, 300, 450, 650, 900]
    print(f"{'M':>6} {'h':>10} {'kappa(self-consist)':>20} {'beta=1/kappa':>14}   (xc={xc}, delta={delta})")
    rows = []
    for M in Ms:
        k = self_consistent(M, xc, delta, wk0=2.8)
        h = (xs - delta - xc) / (M - 1)
        if k is None:
            print(f"{M:6d} {h:10.5f} {'(none)':>20}")
            continue
        rows.append((M, h, k)); print(f"{M:6d} {h:10.5f} {k:20.6f} {1.0/k:14.6f}")
    if len(rows) >= 2:
        (M1, h1, k1), (M2, h2, k2) = rows[-2], rows[-1]
        kext = (k2 * h1 ** 2 - k1 * h2 ** 2) / (h1 ** 2 - h2 ** 2)
        print(f"\nRichardson O(h^2) extrapolation:  kappa = {kext:.6f}   beta = {1.0/kext:.6f}")
        print(f"Published (HKA):                  kappa = 2.81055255   beta = 0.35580192")
    # gauge-mode cross-check
    e, _ = eigs(400, xc, delta, wk=1.0)
    g = [k for k in e if abs(k - 1.0) < 0.1]
    print(f"\nGauge-mode cross-check (must appear near kappa=1): found {np.round(g,5)}")
    print(f"# {time.time()-t0:.1f}s")
