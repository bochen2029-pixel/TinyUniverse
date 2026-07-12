# Scan V0 along the ingoing sonic cone; watch the desingularized-flow eigenstructure at the sonic
# point. The physical (analytic) sonic passage needs a REAL positive eigenvalue whose eigendirection
# integrates to the regular center. Locate where the complex pair -> real (the analytic branch).
import numpy as np
import hka_desing as D

print(f"{'V0':>8} {'A0':>8} {'N0':>8} {'om0':>8} | eigenvalues of desing Jacobian")
for V0 in np.arange(-0.55, -0.03, 0.02):
    try:
        Y0, J, w, Vc = D.sonic_jacobian(V0)
    except Exception as e:
        print(f"{V0:8.3f}  fail {e}"); continue
    # sort by real part desc
    idx = np.argsort(-w.real)
    w = w[idx]
    tag = "COMPLEX" if np.max(np.abs(w.imag)) > 1e-6 else "real"
    ws = " ".join(f"{wi.real:+.4f}{'' if abs(wi.imag)<1e-9 else f'{wi.imag:+.4f}j'}" for wi in w)
    print(f"{V0:8.3f} {Y0[0]:8.4f} {Y0[1]:8.4f} {Y0[2]:8.4f} | {tag:8s} {ws}")
