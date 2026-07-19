# nr_laurent.py — EXACT Laurent series of L(t;kappa)=M_x^{-1}(Gmat-kappa M_s) about the sonic point,
# by power-series arithmetic on the VERIFIED operator coefficients (hka_pert_hka99) evaluated on the
# EXACT background series (nr_sonic). No DFT. Gives R (t^-1 residue) and L_0,L_1,... to machine
# precision -> pins the indicial exponent mu and feeds an accurate Frobenius recursion.
#
# Structure: M_x = blkdiag(I2, [[Ax,Bx],[Cx,Dx]]); its fluid block is SINGULAR at t=0. With
#   det_fl(t) = Ax*Dx-Bx*Cx = t*(d1 + d2 t + ...)   (d1 != 0),   adj_fl(t)=[[Dx,-Bx],[-Cx,Ax]],
# M_fl^{-1} = adj_fl/det_fl = (1/t) * P(t),  P = adj_fl * (1/(det_fl/t)).  Then the fluid rows of L are
#   (1/t) P(t) . (Gmat-kappa M_s)[fluid rows](t) = (1/t) Q(t),   R_fluid = Q_0, L_j fluid = Q_{j+1}.
# Metric rows are regular: L_j metric = (Gmat-kappa M_s)[metric rows]_j.
import numpy as np, math
import nr_sonic as NS
g = 4.0/3.0
smul = NS.smul; sinv = NS.sinv; sconst = NS.sconst; sderiv = NS.sderiv_shift

def _series_coeffs(kappa, order):
    """All operator coefficients as Taylor series in t (length K), on the exact background."""
    K = order + 4
    A, N, om, V = NS.bg_series(K-1)
    A = np.array(A[:K]); N = np.array(N[:K]); om = np.array(om[:K]); V = np.array(V[:K])
    # pad to length K
    def pad(a):
        b = np.zeros(K, complex); b[:len(a)] = a; return b
    A, N, om, V = pad(A), pad(N), pad(om), pad(V)
    obx = smul(sderiv(om), sinv(om)); Vx = sderiv(V)
    V2 = smul(V, V); NV = smul(N, V); NpV = N + V
    oV2 = sconst(K, 1.0) - V2; S = sinv(oV2); S2 = smul(S, S)
    one = sconst(K, 1.0)
    c = {}
    c['As'] = one.copy(); c['Bs'] = g*smul(V, S); c['Cs'] = (g-1)*V; c['Ds'] = g*S
    c['Ax'] = one + NV; c['Bx'] = g*smul(NpV, S); c['Cx'] = (g-1)*NpV; c['Dx'] = g*smul(one+NV, S)
    AN = smul(A, N); ANV = smul(A, NV); omNV = smul(om, NV); omN = smul(om, N)
    c['E1'] = -((g+2)/2)*ANV
    c['E2'] = ((6-3*g)/2)*NV - ((2+g)/2)*ANV + (2-g)*omNV - smul(NV, obx) - g*smul(smul(N, Vx), S)
    c['E3'] = (2-g)*omNV
    c['E4'] = ((6-3*g)/2)*N - ((2+g)/2)*AN + (2-g)*omN - smul(N, obx) \
              - g*smul(smul(one + 2*NV + V2, Vx), S2)
    c['F1'] = ((2-3*g)/2)*AN
    c['F2'] = (2-g)*(g-1)*omN + ((7*g-6)/2)*N + ((2-3*g)/2)*AN - (g-1)*smul(N, obx) \
              - g*smul(smul(NV, Vx), S)
    c['F3'] = (2-g)*(g-1)*omN
    c['F4'] = -(g-1)*obx - g*smul(smul(N + 2*V + smul(N, V2), Vx), S2)
    c['G1'] = -A; c['G3'] = 2*smul(one + (g-1)*V2, smul(om, S)); c['G4'] = 4*g*smul(smul(om, V), S2)
    c['H1'] = A; c['H3'] = (g-2)*om
    return c, K

def laurent_exact(kappa, order=12):
    """Return [R, L0, L1, ..., L_order] (4x4 complex arrays), exact."""
    c, K = _series_coeffs(kappa, order); k = complex(kappa)
    Z = sconst(K, 0.0)
    # Gmat - kappa*Ms as 4x4 of series
    Gm = [[c['G1'], Z, c['G3'], c['G4']],
          [c['H1'], Z, c['H3'], Z],
          [c['E1'], c['E2'], c['E3'], c['E4']],
          [c['F1'], c['F2'], c['F3'], c['F4']]]
    Ms = [[Z, Z, Z, Z], [Z, Z, Z, Z], [Z, Z, c['As'], c['Bs']], [Z, Z, c['Cs'], c['Ds']]]
    GmM = [[Gm[i][j] - k*Ms[i][j] for j in range(4)] for i in range(4)]   # (Gmat - k Ms)(t)
    # fluid block inverse as Laurent: M_fl^{-1} = (1/t) P(t)
    Ax, Bx, Cx, Dx = c['Ax'], c['Bx'], c['Cx'], c['Dx']
    det = smul(Ax, Dx) - smul(Bx, Cx)                 # ~ t * (d1 + ...)
    assert abs(det[0]) < 1e-9, f"det_fl(0)={det[0]} not ~0"
    Dsh = det[1:].copy()                              # det/t = d1 + d2 t + ...  (length K-1)
    invD = sinv(Dsh)                                  # 1/(det/t)
    invD = np.concatenate([invD, [0.0]])              # back to length K
    adj = [[Dx, -Bx], [-Cx, Ax]]
    P = [[smul(adj[i][j], invD) for j in range(2)] for i in range(2)]   # 2x2 series, P0=adj0/d1
    # Q(t) = P(t) . GmM[fluid rows 2,3] (2x4)
    Gfl = [[GmM[2+i][j] for j in range(4)] for i in range(2)]
    Q = [[smul(P[i][0], Gfl[0][j]) + smul(P[i][1], Gfl[1][j]) for j in range(4)] for i in range(2)]
    # assemble R and L_j
    def coeffmat(jfl, jmet):
        M = np.zeros((4, 4), complex)
        for j in range(4):
            M[0, j] = GmM[0][j][jmet] if jmet < K else 0.0
            M[1, j] = GmM[1][j][jmet] if jmet < K else 0.0
            M[2, j] = Q[0][j][jfl] if jfl < K else 0.0
            M[3, j] = Q[1][j][jfl] if jfl < K else 0.0
        return M
    R = np.zeros((4, 4), complex)
    R[2, :] = [Q[0][j][0] for j in range(4)]
    R[3, :] = [Q[1][j][0] for j in range(4)]
    Ls = []
    for j in range(order+1):
        M = np.zeros((4, 4), complex)
        M[0, :] = [GmM[0][jj][j] for jj in range(4)]
        M[1, :] = [GmM[1][jj][j] for jj in range(4)]
        M[2, :] = [Q[0][jj][j+1] for jj in range(4)]
        M[3, :] = [Q[1][jj][j+1] for jj in range(4)]
        Ls.append(M)
    return R, Ls

def analytic_modes(kappa, order=14):
    """3 analytic Frobenius modes a_0 in ker(R), recursion (nI-R)a_n = sum_j L_j a_{n-1-j}."""
    R, Ls = laurent_exact(kappa, order)
    U, S, Vh = np.linalg.svd(R)
    ker = [np.conj(Vh[i]) for i in range(4) if abs(S[i]) < 1e-7*max(S.max(), 1)]
    modes = []
    for a0 in ker:
        a = [np.array(a0, complex)]
        for n in range(1, order+1):
            rhs = np.zeros(4, complex)
            for j in range(0, n):
                if j < len(Ls): rhs += Ls[j].dot(a[n-1-j])
            a.append(np.linalg.solve(n*np.eye(4) - R, rhs))
        modes.append(a)
    return modes, R, Ls

def bg_fld(t, order=18):
    A, N, om, V = NS.bg_series(order)
    def ev(a): return sum(a[k]*t**k for k in range(len(a)))
    def dev(a): return sum(k*a[k]*t**(k-1) for k in range(1, len(a)))
    omv = ev(om)
    return (ev(A), ev(N), omv, ev(V), dev(om)/omv, dev(V))

def gate(kappa=2.81055255, order=18, tm=-0.03, verbose=True, tol=1e-8):
    """Frobenius GATE: each analytic mode must solve the DIRECT operator ODE Psi'=L Psi at t_m.
    tol is the pass threshold — 1e-8 at the gated t_m=-0.02/-0.03; larger |t_m| are truncation-limited
    (finite series radius ~0.12), gate them with a matching looser tol or not at all."""
    import hka_pert_hka99 as H99
    modes, R, Ls = analytic_modes(kappa, order)
    fld = bg_fld(tm, order); L = H99.Lnum(fld, complex(kappa))
    worst = 0.0
    for i, a in enumerate(modes):
        Psi = sum(a[n]*tm**n for n in range(len(a)))
        dPsi = sum(n*a[n]*tm**(n-1) for n in range(1, len(a)))
        res = np.linalg.norm(dPsi - L.dot(Psi))/max(np.linalg.norm(Psi), 1e-30)
        worst = max(worst, res)
        if verbose: print(f"  mode{i}: |Psi'-L Psi|/|Psi| at t={tm} = {res:.2e}")
    if verbose:
        print(f"  >>> kappa={kappa} order={order} t_m={tm}: {len(modes)} modes, worst={worst:.2e} "
              f"({'PASS' if worst < tol else 'FAIL'} <{tol:g})")
    return worst

if __name__ == "__main__":
    print("EXACT Laurent residue indicial exponent mu (nonzero eig of R):")
    for kap in (2.81055255, 1.0, 0.35699):
        R, Ls = laurent_exact(kap, 8)
        ev = np.linalg.eigvals(R); nz = ev[np.argmax(np.abs(ev))]
        print(f"  kappa={kap:<10}: eig(R)={np.sort_complex(ev).round(6)}  mu={nz.real:+.6f}  1-2k={1-2*kap:+.6f}")
    print("\nFrobenius analytic-mode GATE (direct operator ODE residual; gated at t_m=-0.02/-0.03):")
    for tm in (-0.02, -0.03):
        gate(2.81055255, order=18, tm=tm)
    print("  (t_m=-0.05 radius-context row - order-18 truncation dominates, looser tol:)")
    gate(2.81055255, order=18, tm=-0.05, tol=1e-5)
