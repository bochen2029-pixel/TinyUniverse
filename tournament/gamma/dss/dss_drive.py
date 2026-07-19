# dss_drive.py — driver: Delta-frozen pinned continuation at configurable truncation.
# usage: python dss_drive.py [M KE KO]      (default 48 14 13)
import numpy as np, sys, time
from scipy.optimize import least_squares
import nr_relax as R
import nr_evolve as E

M_ = int(sys.argv[1]) if len(sys.argv) > 1 else 48
KE_ = int(sys.argv[2]) if len(sys.argv) > 2 else 14
KO_ = int(sys.argv[3]) if len(sys.argv) > 3 else 13
R.configure(M_, KE_, KO_)
print(f"[drive] truncation (M,KE,KO)=({M_},{KE_},{KO_})", flush=True)
R.vacuum_control()

ev = E.Ev(N=800, rmax=60.0)
u13, info = E.seed_pipeline(ev, 0.03751655962597, verbose=False)
rx = R.Relax(Nz=40)
u = R.seed_from_cylinder(rx, info["cyl"], info["fns"]["xi0"], 3.4453)
r0 = rx.residual(u)
print(f"[drive] seed |r| = {np.linalg.norm(r0):.3f}", flush=True)

lo = np.full(rx.nu, -8.0); hi = np.full(rx.nu, 8.0)
lo[0], hi[0] = 3.4453 - 1e-9, 3.4453 + 1e-9        # Delta FROZEN during continuation
spar = rx.sparsity(pin=True)
Be, Bde, Bo, Bdo, Pe, Po = R.bases(3.4453)
Xp_ssh = (rx.unpack(u)[3]@Bo.T)[-1]
c_prev = max(float(np.sqrt(np.mean(Xp_ssh**2))), 0.05)
print(f"[drive] seed SSH RMS(X+) = {c_prev:.4f}", flush=True)

for c in (0.10, 0.18, 0.28, 0.40):
    s = c/c_prev
    F = u[R.NE:].reshape(rx.Nz, R.NPF).copy()
    F[:, R.NE:] *= s                                # X+- harmonics scale ~ amplitude
    F[:, 1:R.NE] *= s*s                             # g deviation ~ amplitude^2
    u2 = u.copy(); u2[R.NE:] = F.ravel()
    t0 = time.time()
    sol = least_squares(lambda uu: rx.residual(uu, pin=(c, 3.0)), u2, method="trf",
                        bounds=(lo, hi), jac_sparsity=spar, tr_solver="lsmr",
                        xtol=1e-15, ftol=3e-16, gtol=1e-14, max_nfev=2200, verbose=0)
    rphys = np.linalg.norm(sol.fun[:-1])
    D2, x2, g2, xp2, xm2 = rx.unpack(sol.x)
    mags = np.abs(xp2).max(axis=0)
    tail = [float(f"{max(mags[2*j], mags[2*j+1]):.4f}") for j in range(R.NO//2)]
    print(f"c={c:.2f}: |r_phys|={rphys:9.4f}  gdev={np.abs(g2[:,1:]).max():.4f}  "
          f"|X+|max={np.abs(xp2).max():.4f}  tail(k=1,3,..)={tail}  "
          f"({time.time()-t0:.0f}s, status={sol.status})", flush=True)
    u = sol.x; c_prev = c
    np.save(f"drive_{M_}_{KE_}_{KO_}_c{c:.2f}.npy", sol.x)
print("[drive] done")
