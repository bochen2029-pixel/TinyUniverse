# linsys_matrix.py — assemble hp'(x) = L(x;k) . hp for the perturbation eigenproblem, symbolically,
# then lambdify L(N,A,om,V, Nb', Ab', omb', Vb', k). Background slopes passed in (they satisfy the
# background ODE). Produces the 4x4 complex matrix used by the shoot.
import sympy as sp, pickle
r = sp.Symbol('r', positive=True)
A0,N0,o0,V0 = sp.symbols('A0 N0 om0 V0')
Ax,Nx,ox,Vx = sp.symbols('A_x N_x om_x V_x')       # background slopes (H') as symbols
As,Ns,os,Vs = sp.symbols('A_s N_s om_s V_s')
k = sp.symbols('k')                                  # kappa
Np_,Ap_,op_,vp_ = sp.symbols('Np Ap opp vpp')        # perturbation amplitudes hp=(Np,Ap,omp,Vp)
NpP,ApP,opP,vpP = sp.symbols('NpP ApP opP vpP')      # their x-derivatives hp'

L = pickle.load(open("linsys.pkl","rb"))
def S(x): return sp.sympify(x)
cx_en=[S(e) for e in L['cx_en']]; cs_en=[S(e) for e in L['cs_en']]
cx_mo=[S(e) for e in L['cx_mo']]; cs_mo=[S(e) for e in L['cs_mo']]
JE_en=[S(e) for e in L['JE_en']]; JE_mo=[S(e) for e in L['JE_mo']]
Axv=S(L['Axv']); Nxv=S(L['Nxv'])

flds=(N0,A0,o0,V0); hp=(Np_,Ap_,op_,vp_); hpP=(NpP,ApP,opP,vpP)
xds=(Nx,Ax,ox,Vx); sds=(Ns,As,os,Vs)

# background-slope substitution: in JE_* the symbols F_x are the background slopes H'; F_s=0.
bgsub = {Ns:0,As:0,os:0,Vs:0}   # background has no s-derivative
# JE_* was d/dF of the full eq; it contains F_x symbols (=H') and F_s(=0). Keep F_x as-is (they ARE
# the background slopes, provided as inputs). Zero the F_s.
JE_en=[e.subs(bgsub) for e in JE_en]
JE_mo=[e.subs(bgsub) for e in JE_mo]

# metric perturbation slopes: Ap' = d/dF[Axv].hp ; Np' = d/dF[Nxv].hp
def linRHS(RHS):
    return sum(sp.diff(RHS,f)*h for f,h in zip(flds,hp))
Ap_expr = sp.cancel(linRHS(Axv))      # = Ap'
Np_expr = sp.cancel(linRHS(Nxv))      # = Np'

# fluid perturbation eqns:  cx . hp' + k cs . hp + JE . hp = 0   (energy, momentum)
def fluidEq(cx,cs,JE):
    lhs = sum(cx[i]*hpP[i] for i in range(4)) + k*sum(cs[i]*hp[i] for i in range(4)) \
          + sum(JE[i]*hp[i] for i in range(4))
    return lhs
E_en = fluidEq(cx_en,cs_en,JE_en)
E_mo = fluidEq(cx_mo,cs_mo,JE_mo)
# substitute NpP=Np_expr, ApP=Ap_expr (metric-determined), solve for (opP,vpP)
E_en = E_en.subs({NpP:Np_expr, ApP:Ap_expr})
E_mo = E_mo.subs({NpP:Np_expr, ApP:Ap_expr})
# 2x2 linear solve for (opP,vpP): coefficients
a11 = sp.cancel(E_en.coeff(opP,1)); a12 = sp.cancel(E_en.coeff(vpP,1))
a21 = sp.cancel(E_mo.coeff(opP,1)); a22 = sp.cancel(E_mo.coeff(vpP,1))
b1  = sp.cancel(-E_en.coeff(opP,0).coeff(vpP,0))
b2  = sp.cancel(-E_mo.coeff(opP,0).coeff(vpP,0))
detf = sp.cancel(a11*a22 - a12*a21)
opP_e = sp.cancel((b1*a22 - a12*b2)/detf)
vpP_e = sp.cancel((a11*b2 - b1*a21)/detf)

# L rows: hp' = (Np', Ap', omp', Vp') each linear in hp=(Np,Ap,omp,Vp). Extract 4x4 matrix.
rows = [Np_expr, Ap_expr, opP_e, vpP_e]
Lmat = sp.Matrix(4,4, lambda i,j: sp.cancel(sp.diff(rows[i], hp[j])))
args = (N0,A0,o0,V0, Nx,Ax,ox,Vx, k)
Lmat_r = Lmat.subs(r,1)
import pickle as pk
pk.dump({'Lmat':sp.srepr(Lmat_r),
         'opP':sp.srepr(sp.cancel(opP_e.subs(r,1))),
         'vpP':sp.srepr(sp.cancel(vpP_e.subs(r,1))),
         'detf':sp.srepr(sp.cancel(detf.subs(r,1)))}, open("Lmat.pkl","wb"))
print("assembled L(x;k).")
# the fluid-sector denominator = the sonic singular factor; check it vanishes at sonic point
den = sp.factor(sp.denom(sp.together(vpP_e.subs(r,1))))
print("den(vpP) factored =", den)
pt={V0:1/sp.sqrt(3),N0:2/sp.sqrt(3),A0:sp.Rational(3,2),o0:sp.Rational(3,4)}
print("den at sonic point =", sp.simplify(den.subs(pt)))
print("wrote Lmat.pkl")
