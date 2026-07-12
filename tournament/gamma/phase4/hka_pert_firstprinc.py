# hka_pert_firstprinc.py — FIRST-PRINCIPLES perturbation operator L(x;kappa), derived from the covariant
# conservation law nabla_a T^{ab}=0 (fluid_sderiv.pkl, built by derive_linear.py directly from the full
# (t,r) stress-energy with the s-derivative chain rule d/dtau -> (-1/tau)(d_x + d_s)) PLUS the two metric
# slope equations (pure d_x). This is TRANSCRIPTION-FREE: no HKA (5.5)-(5.13), no equations.md eq(18) —
# the operator is the linearization of the SAME equations whose background (hka_ec) is verified to machine
# precision. The s-derivatives (the perturbation kappa-coupling) come from the physics, not a table.
#
# Method: linearize en,mo (fluid, with s-derivs) and the metric eqs about A->A(1+eps ap), N->N(1+eps nn),
# om->om(1+eps oo), V->V+eps vv with the mode ~e^{kappa s}: (log)_s -> eps kappa (profile). Collect the
# LINEAR (O(eps)) coefficient rows. The metric eqs give apx,nnx directly (kappa-free). The two fluid eqs
# are linear in (opx,vvx) [x-derivs] and (ap,nn,oo,vv,apx,nnx) [and kappa via s-terms]; we extract the 2x2
# x-derivative matrix M2 = [[d en/d opx, d en/d vvx],[d mo/d opx, d mo/d vvx]] and the RHS coefficient rows,
# invert M2 NUMERICALLY per background point (fast; no slow symbolic solve). L = P(x) + kappa Q(x), affine.
#
# This module lambdifies the (fast-to-extract) symbolic coefficient pieces once, then Lnum(x,kappa) is
# numeric. Eq #s per HKA_beta_equations.md; physics per KHA95 / equations.md 1.2.

import sympy as sp, numpy as np, pickle, os, math
import hka_beta4 as B

_PKL = "hka_pert_firstprinc_coeffs.pkl"

def _build_symbolic():
    d = pickle.load(open("fluid_sderiv.pkl", "rb"))
    en = sp.sympify(d['en']); mo = sp.sympify(d['mo'])
    Axv = sp.sympify(d['Axv']); Nxv = sp.sympify(d['Nxv'])
    A0, N0, o0, V0 = sp.symbols('A0 N0 om0 V0')
    Ax, Nx, ox, Vx = sp.symbols('A_x N_x om_x V_x')      # background x-derivs (ox = om', Vx = V')
    As, Ns, os_, Vs = sp.symbols('A_s N_s om_s V_s')     # s-derivs
    kappa, eps = sp.symbols('kappa eps')
    ap, nn, oo, vv = sp.symbols('a_p n_p o_p v_p')       # profiles (log for A,N,om; additive for V)
    apx, nnx, opx, vpx = sp.symbols('apx npx opx vpx')   # x-derivs of profiles
    # background x-derivatives (substituted at the end so L depends only on A,N,om,V, ox, Vx).
    # NB: fluid_sderiv.pkl still carries explicit r (from tau=r/E, E=1); r drops out at r=1 (scale-free).
    _r = sp.Symbol('r', positive=True)
    subs_bg = {Ax: Axv, Nx: Nxv, _r: 1}
    en = en.subs(_r, 1); mo = mo.subs(_r, 1)

    def perturb(expr):
        rep = {
            A0: A0 * (1 + eps * ap), N0: N0 * (1 + eps * nn), o0: o0 * (1 + eps * oo), V0: V0 + eps * vv,
            Ax: Ax * (1 + eps * ap) + A0 * eps * apx, Nx: Nx * (1 + eps * nn) + N0 * eps * nnx,
            ox: ox * (1 + eps * oo) + o0 * eps * opx, Vx: Vx + eps * vpx,
            As: A0 * eps * kappa * ap, Ns: N0 * eps * kappa * nn, os_: o0 * eps * kappa * oo, Vs: eps * kappa * vv}
        return expr.subs(rep, simultaneous=True)

    lin = lambda R: sp.expand(sp.diff(R, eps).subs(eps, 0))
    Len = sp.expand(lin(perturb(en)).subs(subs_bg))
    Lmo = sp.expand(lin(perturb(mo)).subs(subs_bg))

    # metric perturbation rows -> apx, nnx (kappa-free), solved symbolically (cheap: single linear eq each)
    Ma = (Ax * (1 + eps * ap) + A0 * eps * apx) / (A0 * (1 + eps * ap)) - \
         (1 - A0 * (1 + eps * ap) + 2 * o0 * (1 + eps * oo) * (1 + (V0 + eps * vv) ** 2 / 3) / (1 - (V0 + eps * vv) ** 2))
    Mb = (Nx * (1 + eps * nn) + N0 * eps * nnx) / (N0 * (1 + eps * nn)) - \
         (-2 + A0 * (1 + eps * ap) - sp.Rational(2, 3) * o0 * (1 + eps * oo))
    LMa = sp.expand(lin(Ma).subs(subs_bg)); LMb = sp.expand(lin(Mb).subs(subs_bg))
    apx_s = sp.solve(LMa, apx)[0]; nnx_s = sp.solve(LMb, nnx)[0]

    # substitute apx,nnx into the fluid rows (cheap subs, NOT a coupled solve)
    Len = sp.expand(Len.subs({apx: apx_s, nnx: nnx_s}))
    Lmo = sp.expand(Lmo.subs({apx: apx_s, nnx: nnx_s}))

    # Now Len,Lmo are linear in opx,vpx and the profiles (+kappa). Extract 2x2 x-derivative matrix and the
    # profile/kappa coefficient rows by .coeff (fast). We invert the 2x2 NUMERICALLY later.
    prof = [ap, nn, oo, vv]
    # M2 (2x2): coeff of opx,vpx in [Len,Lmo]
    M2 = sp.Matrix([[Len.coeff(opx), Len.coeff(vpx)], [Lmo.coeff(opx), Lmo.coeff(vpx)]])
    # RHS rows: -(kappa-free profile part) and -(kappa profile part). fluid eq: M2 (opx,vpx)^T + C0.prof + kappa Ck.prof = 0
    # => (opx,vpx)^T = -M2^{-1}(C0 + kappa Ck).prof
    def prof_row(expr, withk):
        e = expr.coeff(kappa, 1) if withk else expr.coeff(kappa, 0)
        return [sp.simplify(e.coeff(p)) for p in prof]
    C0_en = prof_row(Len, False); Ck_en = prof_row(Len, True)
    C0_mo = prof_row(Lmo, False); Ck_mo = prof_row(Lmo, True)
    # apx,nnx profile rows (kappa-free): P rows 0,1
    apx_row = [sp.simplify(sp.expand(apx_s).coeff(p)) for p in prof]
    nnx_row = [sp.simplify(sp.expand(nnx_s).coeff(p)) for p in prof]

    args = (A0, N0, o0, V0, ox, Vx)
    lam = lambda e: sp.lambdify(args, e, 'numpy')
    data = dict(
        M2=[[lam(M2[0, 0]), lam(M2[0, 1])], [lam(M2[1, 0]), lam(M2[1, 1])]],
        C0_en=[lam(x) for x in C0_en], Ck_en=[lam(x) for x in Ck_en],
        C0_mo=[lam(x) for x in C0_mo], Ck_mo=[lam(x) for x in Ck_mo],
        apx_row=[lam(x) for x in apx_row], nnx_row=[lam(x) for x in nnx_row],
    )
    return data

_D = None
def _load():
    global _D
    if _D is None:
        _D = _build_symbolic()
    return _D

def PQ(x):
    """Return (P,Q) 4x4 numeric blocks of L = P + kappa Q at physical x along the EC background."""
    D = _load()
    A, N, om, V, ombar_x, Vx = B.bg_state(x)
    ox = ombar_x * om     # om' (not log)
    a = (A, N, om, V, ox, Vx)
    M2 = np.array([[D['M2'][0][0](*a), D['M2'][0][1](*a)],
                   [D['M2'][1][0](*a), D['M2'][1][1](*a)]], complex)
    M2i = np.linalg.inv(M2)
    C0 = np.array([[D['C0_en'][j](*a) for j in range(4)],
                   [D['C0_mo'][j](*a) for j in range(4)]], complex)
    Ck = np.array([[D['Ck_en'][j](*a) for j in range(4)],
                   [D['Ck_mo'][j](*a) for j in range(4)]], complex)
    # (opx,vpx) rows = -M2i (C0 + kappa Ck) . prof
    Pfluid = -M2i.dot(C0)      # 2x4  (kappa^0)
    Qfluid = -M2i.dot(Ck)      # 2x4  (kappa^1)
    P = np.zeros((4, 4), complex); Q = np.zeros((4, 4), complex)
    P[0] = [D['apx_row'][j](*a) for j in range(4)]
    P[1] = [D['nnx_row'][j](*a) for j in range(4)]
    P[2] = Pfluid[0]; P[3] = Pfluid[1]
    Q[2] = Qfluid[0]; Q[3] = Qfluid[1]
    return P, Q

def Lnum(x, kappa):
    P, Q = PQ(x)
    return P + complex(kappa) * Q

if __name__ == "__main__":
    import time
    B.bg(); B.bg_path(); xs = B.bg()['xs']
    t0 = time.time()
    P, Q = PQ(-1.0)
    print(f"# first-principles operator built in {time.time()-t0:.1f}s")
    A, N, om, V, ox, Vx = B.bg_state(-1.0); g = 4 / 3; oV2 = 1 - V * V
    print("Abar row P[0]:", np.round(P[0].real, 5),
          " expect (G1,0,G3,G4):", np.round([-A, 0, 2 * (1 + (g - 1) * V * V) * om / oV2, 4 * g * om * V / oV2 ** 2], 5))
    print("Nbar row P[1]:", np.round(P[1].real, 5), " expect (A,0,(g-2)om,0):", np.round([A, 0, (g - 2) * om, 0], 5))
    print("Q Abar/Nbar rows (should be 0):", np.round(Q[0].real, 5), np.round(Q[1].real, 5))
    # exponents
    def R(kappa):
        def tL(t): return t * Lnum(xs + t, kappa)
        t1, t2 = -1e-4, -5e-5; v1, v2 = tL(t1), tL(t2)
        return v2 + (v2 - v1) * ((0 - t2) / (t2 - t1))
    print("center eig(L) k=2.81:", np.round(np.sort(np.linalg.eigvals(Lnum(-14.0, 2.81055255)).real), 4), " want {-2,-1,0,0}")
    print("sonic eig(R) k=2.81:", np.round(np.sort(np.linalg.eigvals(R(2.81055255)).real), 4), " want {0,0,0,-4.6211}")
    import hka_pert as PT
    for kbar in [1.0, 2.81055255]:
        res = []
        for xx in [-8, -4, -1, xs - 0.05]:
            Psi = PT.gauge_mode(xx, kbar, B.bg_state); h = 1e-4
            dP = (PT.gauge_mode(xx + h, kbar, B.bg_state) - PT.gauge_mode(xx - h, kbar, B.bg_state)) / (2 * h)
            res.append(np.linalg.norm(dP - Lnum(xx, kbar).dot(Psi)) / np.linalg.norm(dP))
        print(f"gauge residual kbar={kbar}: {[f'{r:.1e}' for r in res]}")
