# newton_drive.py — EXACT-NEWTON relaxation solver (session-2 recipe #1): finite-difference
# sparse Jacobian on the known block-tridiagonal structure, DIRECT sparse linear solves
# (no lsmr trust-region regularization — measured: trf+lsmr stalls at xtol while the basis
# holds physical structure), damped backtracking line search. This is the solver class
# Gundlach actually used (app D: "relaxation algorithm [Press]" = Newton on the discretized
# algebraic system).
# usage: python newton_drive.py [M KE KO]        (default 48 14 13)
import numpy as np, sys, time
from scipy.sparse import csr_matrix, vstack
from scipy.sparse.linalg import spsolve
from scipy.optimize._numdiff import approx_derivative, group_columns
import nr_relax as R
import nr_evolve as E

def newton(rx, u0, pin=None, freeze_delta=None, max_iter=30, rtol=1e-11, verbose=True):
    """Damped exact Newton. pin=(c,w) appends the amplitude row (overdetermined by 1 ->
    normal-equations step). freeze_delta: hold u[0] fixed (drop its column)."""
    spar = rx.sparsity(pin=(pin is not None)).tocsc()
    groups = group_columns(spar)
    u = u0.copy()
    if freeze_delta is not None:
        u[0] = freeze_delta
    fun = lambda uu: rx.residual(uu, pin=pin)
    r = fun(u)
    hist = [np.linalg.norm(r)]
    for it in range(max_iter):
        J = approx_derivative(fun, u, method='2-point', rel_step=1e-7,
                              sparsity=(spar, groups))
        J = csr_matrix(J)
        if freeze_delta is not None:
            keep = np.ones(rx.nu, bool); keep[0] = False
            Jk = J[:, np.where(keep)[0]]
        else:
            Jk = J
        # solve for the Newton step: square -> direct; overdetermined (pin) -> normal eqs
        if Jk.shape[0] == Jk.shape[1]:
            try:
                dx = spsolve(Jk.tocsc(), -r)
            except Exception:
                dx = spsolve((Jk.T@Jk + 1e-12*csr_matrix(np.eye(Jk.shape[1]))).tocsc(), -Jk.T@r)
        else:
            JT = Jk.T.tocsr()
            dx = spsolve((JT@Jk).tocsc(), -JT@r)
        # backtracking line search
        alpha, ok = 1.0, False
        n0 = np.linalg.norm(r)
        for _ in range(12):
            u_try = u.copy()
            if freeze_delta is not None:
                u_try[1:] = u[1:] + alpha*dx
            else:
                u_try = u + alpha*dx
            r_try = fun(u_try)
            if np.linalg.norm(r_try) < n0*(1.0 - 1e-4*alpha):
                u, r, ok = u_try, r_try, True
                break
            alpha *= 0.5
        n1 = np.linalg.norm(r)
        hist.append(n1)
        if verbose:
            print(f"    it{it:2d}: |r|={n1:.3e}  alpha={alpha:.3f}{'' if ok else '  (NO DECREASE)'}", flush=True)
        if not ok:
            return u, n1, 'stall', hist
        if n1 < rtol:
            return u, n1, 'converged', hist
        if it > 3 and n1 > 0.999*hist[-4]:
            return u, n1, 'slow', hist
    return u, np.linalg.norm(r), 'maxiter', hist

if __name__ == "__main__":
    M_ = int(sys.argv[1]) if len(sys.argv) > 1 else 48
    KE_ = int(sys.argv[2]) if len(sys.argv) > 2 else 14
    KO_ = int(sys.argv[3]) if len(sys.argv) > 3 else 13
    R.configure(M_, KE_, KO_)
    print(f"[newton] truncation (M,KE,KO)=({M_},{KE_},{KO_})", flush=True)
    R.vacuum_control()
    ev = E.Ev(N=800, rmax=60.0)
    u13, info = E.seed_pipeline(ev, 0.03751655962597, verbose=False)
    rx = R.Relax(Nz=40)
    u = R.seed_from_cylinder(rx, info["cyl"], info["fns"]["xi0"], 3.4453)
    print(f"[newton] seed |r| = {np.linalg.norm(rx.residual(u)):.3f}", flush=True)
    DFROZEN = 3.4453
    for c in (0.10, 0.18, 0.28, 0.40):
        t0 = time.time()
        u, rn, status, hist = newton(rx, u, pin=(c, 3.0), freeze_delta=DFROZEN, max_iter=25)
        D2, x2, g2, xp2, xm2 = rx.unpack(u)
        mags = np.abs(xp2).max(axis=0)
        tail = [float(f"{max(mags[2*j], mags[2*j+1]):.4f}") for j in range(R.NO//2)]
        print(f"c={c:.2f}: {status:9s} |r|={rn:.3e}  gdev={np.abs(g2[:,1:]).max():.4f}  "
              f"tail={tail}  ({time.time()-t0:.0f}s)", flush=True)
        np.save(f"newton_{M_}_{KE_}_{KO_}_c{c:.2f}.npy", u)
    print("[newton] ladder done")
