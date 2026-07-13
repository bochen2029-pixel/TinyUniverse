# nr_evans_diag.py — localize why the Evans match is flat. Use the KNOWN gauge mode (exact solution,
# center-regular AND sonic-analytic at all k) as ground truth:
#   (1) is the gauge mode in the center 2-plane at x_c?           (tests center-mode selection)
#   (2) is n_sonic . gauge(x_match) ~ 0 ?                          (tests the sonic 3-plane / backward integ)
#   (3) does forward-integrating gauge(x_c) land on gauge(x_match)? (tests the center forward integ)
import numpy as np
import nr_evans2 as EV
import nr_laurent as NL

XS = EV.XS
for k in [1.0, 2.0, 2.81055255, 3.0]:
    xc = -13.0; xm = XS - 0.5
    u1, u2 = EV._mode_init(k, xc)
    Q0, _ = np.linalg.qr(np.column_stack([u1, u2]))
    g_c = EV.gauge_at(xc, k); g_c /= np.linalg.norm(g_c)
    rej_gc = np.linalg.norm(g_c - Q0.dot(Q0.conj().T.dot(g_c)))
    # sonic 3-plane at xm (analytic modes at t0=-0.03 integrated backward)
    modes, R, Ls = NL.analytic_modes(k, 18)
    s_t0 = [sum(a[n]*(-0.03)**n for n in range(len(a))) for a in modes]
    Us = EV._integ_subspace(s_t0, k, XS-0.03, xm)
    U, S, Vh = np.linalg.svd(np.column_stack([Us[:, 0], Us[:, 1], Us[:, 2]]))
    ns = U[:, 3]
    g_m = EV.gauge_at(xm, k); g_m /= np.linalg.norm(g_m)
    dot_gm = abs(np.vdot(ns.conj(), g_m))
    # forward-integrate gauge(x_c) to xm
    G = EV._integ_subspace([EV.gauge_at(xc, k)], k, xc, xm)
    align = abs(np.vdot(G[:, 0], g_m))
    print(f"k={k:<11}: gauge-in-center-plane rej={rej_gc:.2e} | n_sonic.gauge(xm)={dot_gm:.2e} | "
          f"fwd-integ gauge align={align:.6f}")
