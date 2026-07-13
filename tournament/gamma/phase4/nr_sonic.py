# nr_sonic.py — EXACT background sonic Taylor series for the Evans-Coleman radiation-fluid CSS
# solution, by an order-by-order power-series ("jet") recursion about the sonic point x_s.
#
# WHY: the Stage-B Frobenius/Laurent machinery needs the background (A,N,om,V)(t), t=x-x_s, as an
# ACCURATE analytic series. hka_pert_sonic.bg_series_near_sonic fits a one-sided polyfit of the
# SINGULAR resolved ODE right up against the sonic point -> inaccurate -> the residue mu of L was
# never pinned (three disagreeing values). This module computes the coefficients ANALYTICALLY.
#
# STRUCTURE (this is the crux): the reduced 3-D background is
#     N' = N(-2 + A - (2-g)om),   A = A_of(N,om,V)  [constraint 4.2, algebraic]
#     M(N,om,V) (om',V')^T = b(N,om,V)              [fluid pair 4.1d,e]
# with M's determinant  = g*om*(1-V^2)*Dson(N,V)  -> 0 at the sonic point (M_0 is RANK 1).
# Expand Y(t)=Y0+sum_k c_k t^k. The metric row gives N_{m+1} explicitly. The fluid row at order m,
#     M_0 (m+1) w_{m+1} = r_m := b_m - sum_{i>=1} M_i (m-i+1) w_{m-i+1},   w=(om,V),
# is SINGULAR: it needs the solvability  n_L^T r_m = 0  (this is the L'Hopital cond 4.6, and it is
# what FIXES the free null-direction coefficient alpha_m carried in w_m), and then determines
#     w_{m+1} = pinv(M_0) r_m/(m+1) + alpha_{m+1} n_R   (alpha_{m+1} free, fixed at the next order).
# We resolve alpha_m each step by evaluating r_m affinely in alpha_m (two evals). Pure numpy.
#
# VALIDATED in __main__ against the accurately-integrated hka_ec background (must agree ~1e-10) and
# the c_1 slope is cross-checked against the desingularized-flow eigenvector (hka_desing).
import numpy as np, math
import hka_ec as E

g = 4.0/3.0

# ---- truncated power-series arithmetic (coeff arrays, index k = t^k) ----
def sconst(K, val=0.0):
    a = np.zeros(K, complex); a[0] = val; return a
def smul(a, b):
    K = len(a); c = np.zeros(K, complex)
    for i in range(K):
        ai = a[i]
        if ai == 0: continue
        for j in range(K - i):
            c[i+j] += ai * b[j]
    return c
def sinv(a):
    K = len(a); b = np.zeros(K, complex); b[0] = 1.0/a[0]
    for n in range(1, K):
        s = 0.0
        for k in range(1, n+1): s += a[k]*b[n-k]
        b[n] = -s/a[0]
    return b
def sderiv_shift(a):
    """coeffs of d/dt of series a: (a')_k = (k+1) a_{k+1}."""
    K = len(a); d = np.zeros(K, complex)
    for k in range(K-1): d[k] = (k+1)*a[k+1]
    return d

def _fluid_MB(A, N, om, V):
    """Series (matrix M11,M12,M21,M22 and rhs b1,b2) of the om-cleared fluid pair, from Y-series."""
    NV  = smul(N, V)
    V2  = smul(V, V)
    oneMV2 = sconst(len(N), 1.0) - V2
    S   = sinv(oneMV2)                    # 1/(1-V^2)
    NpV = N + V
    omS = smul(om, S)
    oneP_NV = sconst(len(N), 1.0) + NV
    M11 = oneP_NV
    M12 = g*smul(NpV, omS)
    M21 = (g-1.0)*NpV
    M22 = g*smul(oneP_NV, omS)
    RHS_d = (3*(2-g)/2)*NV - ((2+g)/2)*smul(A, NV) + (2-g)*smul(NV, om)
    RHS_e = (2-g)*(g-1)*smul(N, om) + ((7*g-6)/2)*N + ((2-3*g)/2)*smul(A, N)
    b1 = smul(om, RHS_d); b2 = smul(om, RHS_e)
    return (M11, M12, M21, M22), (b1, b2)

def _A_series(N, om, V):
    V2 = smul(V, V); NV = smul(N, V)
    S = sinv(sconst(len(N), 1.0) - V2)
    term1 = 2.0*smul(om, smul(sconst(len(N), 1.0) + (g-1)*V2, S))
    term2 = 2.0*g*smul(NV, smul(om, S))
    return sconst(len(N), 1.0) + term1 + term2

def _c1_ref_fluid():
    """(om_1,V_1) of c_1 from the desingularized-flow analytic eigenvector (for BRANCH SELECTION at
    the order-1 quadratic; the branch's own quadratic root sets the precise value)."""
    import hka_desing as D
    V0f = -1/math.sqrt(3)
    _, _, w, Vc = D.sonic_jacobian(V0f)
    Y0a = D.sonic_point(V0f); eps = 1e-7; dg = np.zeros(4)
    for j in range(4):
        Yp = Y0a.copy(); Yp[j] += eps; Ym = Y0a.copy(); Ym[j] -= eps
        dg[j] = (D.det_fluid(Yp) - D.det_fluid(Ym))/(2*eps)
    idx = int(np.argmax(w.real)); lam = w[idx].real; v = np.real(Vc[:, idx])
    c1 = lam*v/dg.dot(v)
    return (c1[2], c1[3])

def bg_series(order=12):
    """Return (A,N,om,V) Taylor-coeff arrays (length order+1) in t=x-x_s about the EC sonic point.

    The order-1 fluid solvability is QUADRATIC in the null-coefficient alpha_1 (two analytic branches
    through the sonic point) because M_1*w_1 is a product of two alpha_1-linear factors; we solve the
    quadratic and pick the Evans-Coleman branch (nearest the desingularized eigenvector). Orders >=2 are
    linear in alpha_m -> standard solvability."""
    K = order + 2
    A0, N0, om0, V0 = 1.5, 2/math.sqrt(3), 0.75, -1/math.sqrt(3)
    N = sconst(K, N0); om = sconst(K, om0); V = sconst(K, V0)
    A = _A_series(N, om, V)
    (M11, M12, M21, M22), _ = _fluid_MB(A, N, om, V)
    M0 = np.array([[M11[0], M12[0]], [M21[0], M22[0]]], complex)
    U, s, Vh = np.linalg.svd(M0)
    n_R = np.conj(Vh[-1]); n_L = U[:, -1]; M0pinv = np.linalg.pinv(M0)
    c1ref = _c1_ref_fluid()

    def r_of(Nc, omc, Vc):
        Ac = _A_series(Nc, omc, Vc)
        (m11, m12, m21, m22), (b1, b2) = _fluid_MB(Ac, Nc, omc, Vc)
        wpo = sderiv_shift(omc); wpv = sderiv_shift(Vc)
        Phi1 = smul(m11, wpo) + smul(m12, wpv) - b1
        Phi2 = smul(m21, wpo) + smul(m22, wpv) - b2
        return Phi1, Phi2  # r_m = -(Phi1[m],Phi2[m])

    def nLr(m):
        P1, P2 = r_of(N, om, V); return n_L.conj().dot(-np.array([P1[m], P2[m]]))

    p_next = None
    for m in range(0, order+1):
        if m == 1:
            # QUADRATIC solvability n_L^T r_1(alpha_1) = 0: sample at 3 alphas, fit parabola, pick EC branch
            samp = []
            for a in (-2.0, 0.0, 2.0):
                om[1], V[1] = p_next[0]+a*n_R[0], p_next[1]+a*n_R[1]
                samp.append(nLr(1))
            C = samp[1]; Aq = (samp[0]+samp[2]-2*C)/8.0; Bq = (samp[2]-samp[0])/4.0
            roots = np.roots([Aq, Bq, C]) if abs(Aq) > 1e-12 else np.array([-C/Bq])
            best, bd = 0.0, 1e99
            for rt in roots:
                wr = (p_next[0]+rt*n_R[0], p_next[1]+rt*n_R[1])
                d = abs(wr[0]-c1ref[0]) + abs(wr[1]-c1ref[1])
                if d < bd: bd, best = d, rt.real
            om[1], V[1] = p_next[0]+best*n_R[0], p_next[1]+best*n_R[1]
        elif m >= 2:
            # LINEAR solvability n_L^T r_m(alpha_m) = 0
            om[m], V[m] = p_next[0], p_next[1]; r0 = nLr(m)
            om[m], V[m] = p_next[0]+n_R[0], p_next[1]+n_R[1]; r1 = nLr(m)
            denom = r1 - r0
            alpha = -r0/denom if abs(denom) > 1e-11 else 0.0
            om[m], V[m] = p_next[0]+alpha*n_R[0], p_next[1]+alpha*n_R[1]
        if m < order:
            A = _A_series(N, om, V)
            rhsN = smul(N, sconst(K, -2.0) + A - (2.0-g)*om)
            N[m+1] = rhsN[m]/(m+1)
            P1, P2 = r_of(N, om, V); r_m = -np.array([P1[m], P2[m]])
            w_next = (M0pinv.dot(r_m))/(m+1)
            p_next = (w_next[0], w_next[1])
            om[m+1], V[m+1] = 0.0, 0.0
    A = _A_series(N, om, V)
    return A[:order+1], N[:order+1], om[:order+1], V[:order+1]

def eval_bg(coeffs, t):
    A, N, om, V = coeffs
    def ev(a): return sum(a[k]*t**k for k in range(len(a)))
    return np.array([ev(A), ev(N), ev(om), ev(V)])

if __name__ == "__main__":
    import sys
    from scipy.integrate import solve_ivp
    order = int(sys.argv[1]) if len(sys.argv) > 1 else 12
    coeffs = bg_series(order)
    A, N, om, V = coeffs
    print(f"sonic Taylor series (order {order}), t=x-x_s:")
    print(f"  A: {[f'{c.real:+.6g}' for c in A[:6]]}")
    print(f"  N: {[f'{c.real:+.6g}' for c in N[:6]]}")
    print(f"  om:{[f'{c.real:+.6g}' for c in om[:6]]}")
    print(f"  V: {[f'{c.real:+.6g}' for c in V[:6]]}")

    # c_1 cross-check vs desingularized-flow eigenvector
    import hka_desing as D
    Y0, J, w, Vc = D.sonic_jacobian(-1/math.sqrt(3))
    detgrad = np.zeros(4)
    Y0a = D.sonic_point(-1/math.sqrt(3)); eps=1e-6
    for j in range(4):
        Yp=Y0a.copy(); Yp[j]+=eps; detgrad[j]=(D.det_fluid(Yp)-D.det_fluid(Y0a))/eps
    print("\nc_1 (series):", np.array([A[1],N[1],om[1],V[1]]).real.round(6))
    # eigenvector with smallest |lambda| that is the analytic (real) direction; report all for eyeball
    for i in range(4):
        vi = np.real_if_close(Vc[:,i]);
        print(f"   desing eig lam={w[i].real:+.4f}{w[i].imag:+.4f}j dir(A,N,om,V)={np.real(vi).round(4)}")

    # VALIDATION vs accurately-integrated EC background
    r = E.shoot_to_sonic(1.0, 0.375, rtol=1e-12, atol=1e-14); xs = r['x']
    z0=math.exp(-16.0); Y0i=E.center_state3(1.0,0.375,z0)
    sol=solve_ivp(E.rhs3,[-16.0, xs-1e-8], Y0i, method='DOP853', rtol=1e-12, atol=1e-14, dense_output=True)
    def refY(x):
        Nr,omr,Vr=sol.sol(x); return np.array([E.A_of(Nr,omr,Vr),Nr,omr,Vr])
    print(f"\nVALIDATION vs integrated EC background (x_s={xs:.6f}):")
    worst=0.0
    for t in [-0.02,-0.05,-0.08,-0.12,-0.18,-0.25]:
        err=np.abs(eval_bg(coeffs,t).real - refY(xs+t)); worst=max(worst,err.max())
        print(f"  t={t:+.3f}: max|series-ref|={err.max():.2e}   ({', '.join(f'{e:.1e}' for e in err)})")
    print(f"  >>> worst over t in [-0.25,-0.02]: {worst:.2e}   (< 1e-9 => series is EXACT-grade)")
