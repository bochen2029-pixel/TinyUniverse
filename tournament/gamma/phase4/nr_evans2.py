# nr_evans2.py — fluid-beta eigenvalue by a TWO-SIDED match away from the sonic singularity.
#
# The non-analytic sonic mode ~ |t|^{1-2k} (=|t|^{-4.62} at k=2.81) BLOWS UP as t->0, so forward-shooting
# the center modes to x_m near the sonic point amplifies non-analytic noise ~1e12 and destroys the
# eigenvalue signal. Cure (standard singular-BVP two-sided shoot):
#   * CENTER 2-plane: forward-integrate the two bounded center modes from x_c to an INTERIOR x_match
#     (clean: the excluded e^{-2x},e^{-x} modes DECAY forward).
#   * SONIC 3-plane: take the 3 EXACT analytic Frobenius modes at t0=-0.03 (nr_laurent) and integrate
#     them BACKWARD to x_match (clean: the non-analytic mode DECAYS backward).
# Match at x_match where the amplification |t_match|^{1-2k} is O(10) -> the eigenvalue dip is wide/resolvable.
# Eigenvalue condition (gauge-robust): E(k) = n_sonic . c_phys = 0, n_sonic = normal to the sonic
# 3-plane at x_match, c_phys = the center-regular direction with the KNOWN gauge mode removed.
import numpy as np
import hka_pert_core as PC, hka_pert_hka99 as H99
PC.Lnum = H99.Lnum
import hka_beta4 as B
import nr_laurent as NL
from scipy.integrate import solve_ivp

B.bg(); B.bg_path(); XS = B.bg()['xs']

def L4(x, k): return H99.Lnum(B.bg_state(x), complex(k))

def gauge_at(x, k):
    A, N, om, V, obx, Vx = [float(np.real(z)) for z in B.bg_state(x)]
    ABx, NBx, OBx, VBx = H99._bg_derivs(A, N, om, V, obx, Vx)
    return np.array([ABx, NBx+k, OBx, VBx], complex)

def _mode_init(k, xc):
    L = L4(xc, k); ev, Vr = np.linalg.eig(L); idx = np.argsort(ev.real)
    return Vr[:, idx[2]], Vr[:, idx[3]]     # two bounded (Re~0) center modes

def _integ_subspace(cols, k, xa, xb, nseg=200):
    """Integrate a set of column-vectors from xa to xb, QR-reorthonormalizing each segment (span-stable)."""
    U = np.column_stack(cols).astype(complex); m = U.shape[1]
    xs = np.linspace(xa, xb, nseg+1)
    for i in range(nseg):
        def rhs(x, y):
            Y = y.reshape(m, 4).T; return L4(x, k).dot(Y).T.reshape(-1)
        sol = solve_ivp(rhs, [xs[i], xs[i+1]], U.T.reshape(-1), method='DOP853', rtol=1e-9, atol=1e-12)
        U = sol.y[:, -1].reshape(m, 4).T
        U, _ = np.linalg.qr(U)
    return U

def E(k, xc=-13.0, t0=-0.03, t_match=-0.5, order=18):
    xm = XS + t_match
    # center 2-plane forward to x_match
    u1, u2 = _mode_init(k, xc)
    Uc = _integ_subspace([u1, u2], k, xc, xm)
    # sonic 3-plane: analytic modes at t0, integrate BACKWARD to x_match
    modes, R, Ls = NL.analytic_modes(k, order)
    s_at_t0 = [sum(a[n]*t0**n for n in range(len(a))) for a in modes]
    Us = _integ_subspace(s_at_t0, k, XS+t0, xm)      # backward (xm < XS+t0)
    # normal to the sonic 3-plane
    Uf = np.column_stack([Us[:, 0], Us[:, 1], Us[:, 2]])
    U, S, Vh = np.linalg.svd(Uf); ns = U[:, 3]
    # remove the known gauge mode from the center plane
    g = gauge_at(xm, k)
    gp = Uc.dot(np.linalg.lstsq(Uc, g, rcond=None)[0]); gp = gp/np.linalg.norm(gp)
    c = Uc[:, 0] - np.vdot(gp, Uc[:, 0])*gp
    if np.linalg.norm(c) < 1e-3: c = Uc[:, 1] - np.vdot(gp, Uc[:, 1])*gp
    c = c/np.linalg.norm(c)
    return np.vdot(ns.conj(), c)

if __name__ == "__main__":
    import sys, time
    t0 = time.time()
    tm = float(sys.argv[1]) if len(sys.argv) > 1 else -0.5
    print(f"two-sided Evans, x_s={XS:.5f}, t_match={tm}:  E(kappa)=n_sonic . c_phys")
    for k in [1.0, 2.0, 2.5, 2.7, 2.79, 2.81, 2.8106, 2.83, 2.9, 3.0, 3.5]:
        v = E(k, t_match=tm)
        print(f"  kappa={k:<8}: E={v.real:+.4e}{v.imag:+.4e}j  |E|={abs(v):.3e}", flush=True)
    print(f"# {time.time()-t0:.1f}s")
