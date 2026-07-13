# nr_evans.py — gauge-robust Evans function for the fluid-beta eigenvalue, on the VERIFIED operator.
#
# The center-regular 2-plane at the match point x_m contains the KNOWN gauge mode Psi_g (exact at all
# x) plus the physical center-regular direction c_phys. A solution is sonic-analytic iff n_s4 . (it) = 0,
# where n_s4 = normal to the EXACT 3-D analytic span (nr_laurent modes) at x_m. Since n_s4 . Psi_g = 0
# ALWAYS (gauge is analytic), the eigenvalue condition is purely  E(kappa) = n_s4 . c_phys = 0, with
# c_phys = the center-regular direction with the gauge mode removed. Stable forward integration of the
# center 2-plane (the two bounded center modes DECAY the excluded modes forward), no mode collapse.
import numpy as np
import hka_pert_core as PC, hka_pert_hka99 as H99
PC.Lnum = H99.Lnum
import hka_beta4 as B
import nr_laurent as NL
from scipy.integrate import solve_ivp

B.bg(); B.bg_path(); XS = B.bg()['xs']

def L4(x, k): return H99.Lnum(B.bg_state(x), complex(k))

def gauge_at(x, k):
    A, N, om, V, obx, Vx = B.bg_state(x)
    ABx, NBx, OBx, VBx = H99._bg_derivs(float(np.real(A)), float(np.real(N)),
                                        float(np.real(om)), float(np.real(V)), float(np.real(obx)), float(np.real(Vx)))
    return np.array([ABx, NBx+k, OBx, VBx], complex)

def center_plane_init(k, xc):
    """The 2 bounded (Re lambda ~ 0) center modes of L4 at x_c."""
    L = L4(xc, k); ev, Vr = np.linalg.eig(L)
    idx = np.argsort(ev.real)          # most-negative = e^{-2x},e^{-x} blow-ups; largest ~0 = bounded
    return Vr[:, idx[2]], Vr[:, idx[3]], ev[idx]

def integ_plane(u1, u2, k, xc, xm, nseg=400):
    """Forward-integrate the 2-plane with periodic QR re-orthonormalization (prevents collapse)."""
    xs = np.linspace(xc, xm, nseg+1)
    U = np.column_stack([u1, u2]).astype(complex)
    for i in range(nseg):
        a, b = xs[i], xs[i+1]
        def rhs(x, y):
            Y = y.reshape(2, 4).T; dY = L4(x, k).dot(Y); return dY.T.reshape(-1)
        sol = solve_ivp(rhs, [a, b], U.T.reshape(-1), method='DOP853', rtol=1e-9, atol=1e-12)
        U = sol.y[:, -1].reshape(2, 4).T
        U, _ = np.linalg.qr(U)         # re-orthonormalize the span
    return U[:, 0], U[:, 1]

def n_sonic(k, tm, order=18):
    modes, R, Ls = NL.analytic_modes(k, order)
    scols = np.column_stack([sum(a[n]*tm**n for n in range(len(a))) for a in modes])  # 4x3
    U, S, Vh = np.linalg.svd(scols)
    return U[:, 3]                      # left sing vec for the (zero) 4th singular value = normal

def E(k, xc=-13.0, tm=-0.03, order=18):
    xm = XS + tm
    u1, u2, evs = center_plane_init(k, xc)
    u1, u2 = integ_plane(u1, u2, k, xc, xm)
    g = gauge_at(xm, k); g = g/np.linalg.norm(g)
    ns = n_sonic(k, tm, order)
    # remove gauge from the plane: c_phys = plane direction orthogonal to g
    P = np.column_stack([u1, u2])
    gp = P.dot(np.linalg.lstsq(P, g, rcond=None)[0])   # g projected into the plane
    gp = gp/np.linalg.norm(gp)
    c1 = u1 - np.vdot(gp, u1)*gp
    if np.linalg.norm(c1) < 1e-3: c1 = u2 - np.vdot(gp, u2)*gp
    c_phys = c1/np.linalg.norm(c1)
    return np.vdot(ns.conj(), c_phys)   # = n_s . c_phys (scale-free)

if __name__ == "__main__":
    import sys
    print(f"x_s={XS:.5f}, x_m=x_s-0.03. Evans E(kappa)=n_sonic . c_phys (gauge removed):")
    print(f"center-mode indicial (Re) at x_c=-13:", np.round(center_plane_init(2.81, -13.0)[2].real, 3))
    for k in [1.0, 2.5, 2.7, 2.80, 2.81, 2.8106, 2.82, 2.9, 3.2]:
        val = E(k)
        print(f"  kappa={k:<8}: E={val.real:+.4e}{val.imag:+.4e}j  |E|={abs(val):.3e}")
