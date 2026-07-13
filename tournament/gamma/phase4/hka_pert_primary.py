# hka_pert_primary.py — Stage-B perturbation operator L(x;kappa), built by linearizing the
# PRIMARY KHA95 EOM (gr-qc/9503007, eq:EOM) DIRECTLY, keeping d/ds -> kappa.
#
# WHY: the prior operator mixed linear combinations — the reduced flow Jacobian J was built from
# the code's *recombined* Eq4 (hka_ec fluid_slopes: Cx=(g-1)(N+V), Dx=g(1+NV)) while the s-coupling
# As was transcribed from the *primary* Eq3/Eq4. The d/ds coupling is combination-dependent, so
# J and K were inconsistent and the gauge-mode gate failed. Here BOTH come from the SAME primary
# equations, in mixed variables (Abar,Nbar,obar,Vp) = (dlnA, dlnN, dlnom, dV) to match PC/hka_beta4.
#
# Gauge mode (HKA 5.20), mixed vars:  Psi_g = (Abar_x, Nbar_x + kbar, obar_x, V_x).  Must solve
# Psi_g' = L(kbar) Psi_g for EVERY kbar (pure coordinate freedom) -> the exactness gate.
import sympy as sp, numpy as np, math

s, x, eps, k = sp.symbols('s x epsilon kappa')
g = sp.Rational(4, 3)
Ab = sp.Function('Ab'); Nb = sp.Function('Nb'); ob = sp.Function('ob'); Vb = sp.Function('Vb')
Abar = sp.Function('Abar'); Nbar = sp.Function('Nbar'); obar = sp.Function('obar'); Vpf = sp.Function('Vpf')
Es = sp.exp(k*s)
A  = Ab(x)*sp.exp(eps*Abar(x)*Es)
N  = Nb(x)*sp.exp(eps*Nbar(x)*Es)
om = ob(x)*sp.exp(eps*obar(x)*Es)
V  = Vb(x) + eps*Vpf(x)*Es
oV2 = 1 - V**2
Ax = sp.diff(A, x); As = sp.diff(A, s)
Nx = sp.diff(N, x); Ns = sp.diff(N, s)
ox = sp.diff(om, x); os = sp.diff(om, s)
Vx = sp.diff(V, x); Vs = sp.diff(V, s)

# primary KHA95 EOM (each expression = 0):
E1 = Ax/A - (1 - A + 2*om/oV2*(1 + V**2/3))
E2 = Nx/N - (-2 + A - 2*om/3)
E3 = (os + (1 + N*V)*ox)/om + 4*(V*Vs + (N + V)*Vx)/(3*oV2) \
     - N*V*Ax/(3*A) + 4*V*Nx/3 + 2*N*(1 + 4*om/(9*oV2))
E4 = (4*V*os + (4*V + N + 3*N*V**2)*ox)/om + 4*((1 + V**2)*Vs + (1 + V**2 + 2*N*V)*Vx)/oV2 \
     + N*oV2*Ax/A + 4*(1 + V**2)*Nx + 2*N*(1 + 3*V**2)

# background-value + background-derivative symbols; perturbation + its x-deriv symbols
A0, N0, o0, V0 = sp.symbols('A0 N0 o0 V0', positive=True)
ABx, NBx, OBx, VBx = sp.symbols('ABx NBx OBx VBx')      # background dlnA/dx, dlnN/dx, dlnom/dx, dV/dx
Ap, Np, op0, Vp = sp.symbols('Ap Np op Vp')             # Abar, Nbar, obar, V_p
Apx, Npx, opx, Vpx = sp.symbols('Apx Npx opx Vpx')      # their x-derivatives

def lin(Eq):
    d = sp.diff(Eq, eps).subs(eps, 0)
    d = d / Es
    subs = [(sp.Derivative(Abar(x), x), Apx), (sp.Derivative(Nbar(x), x), Npx),
            (sp.Derivative(obar(x), x), opx), (sp.Derivative(Vpf(x), x), Vpx),
            (sp.Derivative(Ab(x), x), A0*ABx), (sp.Derivative(Nb(x), x), N0*NBx),
            (sp.Derivative(ob(x), x), o0*OBx), (sp.Derivative(Vb(x), x), VBx),
            (Abar(x), Ap), (Nbar(x), Np), (obar(x), op0), (Vpf(x), Vp),
            (Ab(x), A0), (Nb(x), N0), (ob(x), o0), (Vb(x), V0)]
    return sp.expand(d.subs(subs))

L1, L2, L3, L4 = [lin(e) for e in (E1, E2, E3, E4)]
sApx = sp.solve(L1, Apx)[0]
sNpx = sp.solve(L2, Npx)[0]
L3b = L3.subs([(Apx, sApx), (Npx, sNpx)])
L4b = L4.subs([(Apx, sApx), (Npx, sNpx)])
sol = sp.solve([L3b, L4b], [opx, Vpx])
sopx, sVpx = sol[opx], sol[Vpx]

# assemble 4x4 L: (Ap,Np,op,Vp) -> (Apx,Npx,opx,Vpx)
rows = [sApx, sNpx, sopx, sVpx]
vars4 = [Ap, Np, op0, Vp]
Lsym = sp.Matrix([[sp.expand(r).coeff(v) for v in vars4] for r in rows])
# sanity: no leftover pert vars in the "constant" part (rows must be linear/homogeneous)
args = (A0, N0, o0, V0, ABx, NBx, OBx, VBx, k)
Lfun = [[sp.lambdify(args, Lsym[i, j], 'numpy') for j in range(4)] for i in range(4)]

def _bg_derivs(A, N, om, V, ombar_x, Vx_):
    ABx_ = 1 - A + 2*om/(1 - V*V)*(1 + V*V/3.0)     # dlnA/dx from Eq1
    NBx_ = -2 + A - 2*om/3.0                         # dlnN/dx from Eq2
    return ABx_, NBx_, ombar_x, Vx_

def Lnum(fld, kappa):
    """fld = (A,N,om,V,ombar_x,V_x) as PC.bg_fields returns. Returns 4x4 complex L."""
    A, N, om, V, ombar_x, Vx_ = fld
    ABx_, NBx_, OBx_, VBx_ = _bg_derivs(A, N, om, V, ombar_x, Vx_)
    a = (A, N, om, V, ABx_, NBx_, OBx_, VBx_, complex(kappa))
    return np.array([[complex(Lfun[i][j](*a)) for j in range(4)] for i in range(4)])

if __name__ == "__main__":
    import hka_beta4 as B
    B.bg(); B.bg_path(); xs = B.bg()['xs']
    print("PRIMARY-EOM operator: gauge-mode exactness gate (residual ~0 for a CORRECT L, all kbar):")
    def gauge(x, kbar):
        A, N, om, V, obx, Vx_ = B.bg_state(x)
        ABx_, NBx_, OBx_, VBx_ = _bg_derivs(A, N, om, V, obx, Vx_)
        return np.array([ABx_, NBx_ + kbar, OBx_, VBx_], complex)
    xlist = [-8.0, -4.0, -1.0, xs - 0.05]
    for kbar in [0.35699, 1.0, 2.81055255]:
        out = []
        for xx in xlist:
            Psi = gauge(xx, kbar)
            dP = (gauge(xx + 1e-6, kbar) - gauge(xx - 1e-6, kbar)) / 2e-6
            res = dP - Lnum(B.bg_state(xx), complex(kbar)).dot(Psi)
            out.append(np.linalg.norm(res) / max(np.linalg.norm(dP), 1e-12))
        print(f"  kbar={kbar:>10}: rel residuals = {[f'{r:.1e}' for r in out]}")
