# sonic_series.py — the sonic-point local structure (the missing 2nd-order piece).
#
# At x=0: Dson(N,V)=0 (singular locus). For local analyticity the Cramer numerators must vanish:
#   numV(A,N,om,V)=0, numOm(A,N,om,V)=0.  With N fixed by Dson=0 given V0, these 2 eqns fix (A0,om0).
#   => sonic-point state is a 1-parameter family labeled by V0=V_ss(0).
#
# 1st-order slopes (L'Hopital / separatrix directions):
#   omp=om', Vp=V' are the finite 0/0 limits. Along the trajectory to 1st order,
#     det   ~ (grad det . slope) x   =: Ddet x
#     numV  ~ (grad numV . slope) x  =: DnV  x     ,  numOm ~ DnOm x
#   with slope=(Np,Ap,omp,Vp), Np,Ap known (regular). Then
#     Vp = numV/det -> Vp = DnV/Ddet ,  omp = DnOm/Ddet.
#   Since Ddet,DnV,DnOm are linear in (omp,Vp), these are 2 coupled quadratics -> a quadratic;
#   two roots = the two separatrix slopes of the saddle. Pick the branch to the regular center.
#
# 2nd order: expand each field Y_i(x)=Y0_i + Y1_i x + Y2_i x^2/2. Y1 from above. For Y2, differentiate
# the ODE Y' = F(Y): Y'' = J(Y).Y'. At x=0, N'',A'' are direct; om'',V'' need the next L'Hopital order
# (ratio of 2nd derivatives), which we get by matching the O(x^2) terms of numV=det*Vp etc. We obtain
# om2,V2 by solving the linear system from d/dx[ numV - det*V' ]=0 and d/dx[ numOm - det*om' ]=0 at 0.

import sympy as sp, pickle, numpy as np
from scipy.optimize import fsolve

d = pickle.load(open("css_syms.pkl", "rb"))
rsym = sp.Symbol('r', positive=True)
# numV,numOm,det share an overall r-power (from alpha ~ r/tau); it cancels in the ratios and does
# not affect the vanishing loci or slope ratios. Set r=1 for evaluation (homogeneous scaling).
S = {k: sp.sympify(v).subs(rsym, 1) for k, v in d.items()}
A0, N0, o0, V0 = sp.symbols('A0 N0 om0 V0')
syms = (A0, N0, o0, V0)

numOm = S['numOm']; numV = S['numV']; det = S['det']
Axv = S['Axv']; Nxv = S['Nxv']; Dson = S['Dson']

# lambdified core
f_numOm = sp.lambdify(syms, numOm, 'numpy')
f_numV  = sp.lambdify(syms, numV, 'numpy')
f_det   = sp.lambdify(syms, det, 'numpy')
f_Dson  = sp.lambdify(syms, Dson, 'numpy')
f_Ax    = sp.lambdify(syms, Axv, 'numpy')
f_Nx    = sp.lambdify(syms, Nxv, 'numpy')
f_omx   = sp.lambdify(syms, S['omx'], 'numpy')
f_Vx    = sp.lambdify(syms, S['Vxx'], 'numpy')

# gradients (as lambdified vectors) of numV,numOm,det wrt (N,A,om,V) — order matches slope vector
def grad(expr):
    gs = [sp.diff(expr, v) for v in (N0, A0, o0, V0)]
    return sp.lambdify(syms, gs, 'numpy')
g_numV  = grad(numV); g_numOm = grad(numOm); g_det = grad(det)

def N_of_V_sonic(V):
    # Dson=0: (3V^2-1)N^2 - 4V N + (3-V^2)=0  -> two roots. Choose the one that is the physical
    # ingoing branch (checked downstream). Return both.
    aa = 3*V*V - 1; bb = -4*V; cc = 3 - V*V
    disc = bb*bb - 4*aa*cc
    if disc < 0: return []
    if abs(aa) < 1e-14:
        return [-cc/bb]
    return [(-bb + np.sqrt(disc))/(2*aa), (-bb - np.sqrt(disc))/(2*aa)]

def sonic_state(V, Nbranch):
    """Given V0 and a chosen N-branch value, solve numV=numOm=0 for (A0,om0)."""
    N = Nbranch
    def res(p):
        A, om = p
        return [float(f_numV(A, N, om, V)), float(f_numOm(A, N, om, V))]
    best = None
    for g in [[1.75, 0.375], [1.5, 0.3], [1.1, 0.05], [2.0, 0.5], [1.3, 0.1],
              [1.46, 0.49], [1.2, 0.2], [1.6, 0.35], [1.0, 0.01], [1.9, 0.6]]:
        s = fsolve(res, g, full_output=True)
        if s[2] == 1:
            A, om = s[0]
            r = max(abs(np.array(res(s[0]))))
            if r < 1e-9 and A > 0 and om > 0:
                if best is None or r < best[2]:
                    best = (A, om, r)
    return best  # (A0,om0,resid) or None

def slopes_1st(V, N, A, om):
    """Two separatrix slope-vectors (Np,Ap,omp,Vp) at the sonic point via L'Hopital quadratic."""
    Np = float(f_Nx(A, N, om, V)); Ap = float(f_Ax(A, N, om, V))
    gV = np.array(g_numV(A, N, om, V), float)   # d numV /d(N,A,om,V)
    gO = np.array(g_numOm(A, N, om, V), float)
    gD = np.array(g_det(A, N, om, V), float)
    # unknowns omp,Vp. slope=(Np,Ap,omp,Vp).
    # DnV = gV.slope, DnOm=gO.slope, Ddet=gD.slope.  Conditions: DnV = Vp*Ddet, DnOm = omp*Ddet.
    op, vp = sp.symbols('op vp', real=True)
    slope = sp.Matrix([Np, Ap, op, vp])
    DnV  = (sp.Matrix(gV).T*slope)[0]
    DnOm = (sp.Matrix(gO).T*slope)[0]
    Ddet = (sp.Matrix(gD).T*slope)[0]
    eqs = [sp.expand(DnV - vp*Ddet), sp.expand(DnOm - op*Ddet)]
    sols = sp.solve(eqs, [op, vp], dict=True)
    out = []
    for s in sols:
        try:
            ov = complex(s[op]); vv = complex(s[vp])
        except Exception:
            continue
        if abs(ov.imag) < 1e-9 and abs(vv.imag) < 1e-9:
            out.append((Np, Ap, float(ov.real), float(vv.real)))
    return out

if __name__ == "__main__":
    print("Sonic one-parameter family (label V0), both N-branches, with separatrix slopes:")
    print(f"{'V0':>7} {'N0':>9} {'A0':>9} {'om0':>9} | separatrix (V\', om\') pairs")
    for V in [-0.30, -0.25, -0.20, -0.18, -0.15, -0.10, 0.0, 0.10, 0.18, 0.25]:
        for Nb in N_of_V_sonic(V):
            if Nb <= 0: continue
            st = sonic_state(V, Nb)
            if st is None: continue
            A, om, r = st
            sl = slopes_1st(V, Nb, A, om)
            slt = "  ".join(f"(V'={s[3]:+.3f},om'={s[2]:+.3f})" for s in sl)
            print(f"{V:7.3f} {Nb:9.4f} {A:9.4f} {om:9.4f} | {slt}")
