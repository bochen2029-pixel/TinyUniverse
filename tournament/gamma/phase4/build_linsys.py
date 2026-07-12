# build_linsys.py — construct the linearized perturbation ODE system for h_p(x) e^{k s}.
# Fields F = (N,A,om,V). Background H(x). Perturbation F = H + eps hp(x) e^{k s}, hp=(Np,Ap,omp,Vp).
#
# The system:
#   (metric) two constraint eqns  M1 := A_x - A*(...) = 0 ,  M2 := N_x - N*(...) = 0   (no s-deriv)
#   (fluid)  en(...) = 0 ,  mo(...) = 0   (carry s-derivs)
#
# Linearize. For metric constraints: A_x = A'(x)+eps(Ap' e^{ks}), etc. O(eps) gives
#   Ap' = d/dF[A*(...)] . hp   -> Ap' = (linear in Np,Ap,omp,Vp).  Similarly Np'.
# For fluid: collect O(eps). ox->omp'e^{ks}, os->k omp e^{ks}, etc. Divide out e^{ks}.
# Result: c_en . (Np',Ap',omp',Vp') + (k * s_en . hp) + (Jalg_en . hp) = 0  (and same for mo),
# where c_en are the x-deriv coeffs (evaluated on background), s_en the s-deriv coeffs, Jalg the
# jacobian of the algebraic part. Combined with metric Ap',Np' we solve the 4x4 for hp'.
#
# We produce hp' = L(x;k) . hp  with L a 4x4 matrix (rational in background fields & k). Dump lambdified.
import sympy as sp, pickle
r = sp.Symbol('r', positive=True)
A0,N0,o0,V0 = sp.symbols('A0 N0 om0 V0')
Ax,Nx,ox,Vx = sp.symbols('A_x N_x om_x V_x')
As,Ns,os,Vs = sp.symbols('A_s N_s om_s V_s')
d = pickle.load(open("fluid_sderiv.pkl","rb"))
en = sp.sympify(d['en']); mo = sp.sympify(d['mo'])
Axv = sp.sympify(d['Axv']); Nxv = sp.sympify(d['Nxv'])

flds = (N0,A0,o0,V0)
xds  = (Nx,Ax,ox,Vx)
sds  = (Ns,As,os,Vs)

# --- fluid: split each eq into x-deriv coeffs, s-deriv coeffs, algebraic remainder ---
def split(eq):
    eq = sp.expand(eq)
    cx = [sp.simplify(eq.coeff(xd,1)) for xd in xds]
    cs = [sp.simplify(eq.coeff(sd,1)) for sd in sds]
    alg = eq
    for xd in xds: alg = alg.coeff(xd,0)   # drop terms with that xd
    for sd in sds: alg = alg.coeff(sd,0)
    alg = sp.simplify(alg)
    return cx, cs, alg
cx_en, cs_en, alg_en = split(en)
cx_mo, cs_mo, alg_mo = split(mo)

# background satisfies alg + cx.H' = 0 with H' the background slopes; but for LINEARIZATION we need
# the Jacobian of the WHOLE fluid eq wrt fields at fixed (treating x-deriv & s-deriv as independent
# perturbation contributions). The perturbation of eq E(F,F_x,F_s)=0:
#   dE = sum_i [dE/dF_i] hp_i + [dE/dF_{x,i}](hp_i' ) + [dE/dF_{s,i}](k hp_i)  (times e^{ks}) = 0
# dE/dF_{x,i} = cx_i ; dE/dF_{s,i} = cs_i ; dE/dF_i = full partial (alg + cx.H'+cs.0 differentiated).
# On the background F_s=0. So dE/dF_i must be computed from E with F_x replaced by background slopes
# AFTER differentiation? No: E depends on F_i explicitly AND through nothing else (F_x,F_s are
# independent symbols). So dE/dF_i = partial of E wrt F_i treating F_x,F_s as constants.
def dEdF(eq):
    return [sp.simplify(sp.diff(eq, f)) for f in flds]
JE_en = dEdF(en)   # includes derivative of cx*F_x and cs*F_s pieces wrt F_i -> but those multiply
JE_mo = dEdF(mo)   # F_x,F_s which on background are (H'_i, 0). We'll substitute background slopes.

pickle.dump({
  'cx_en':[sp.srepr(e) for e in cx_en], 'cs_en':[sp.srepr(e) for e in cs_en],
  'cx_mo':[sp.srepr(e) for e in cx_mo], 'cs_mo':[sp.srepr(e) for e in cs_mo],
  'JE_en':[sp.srepr(e) for e in JE_en], 'JE_mo':[sp.srepr(e) for e in JE_mo],
  'Axv':sp.srepr(Axv),'Nxv':sp.srepr(Nxv),
}, open("linsys.pkl","wb"))
print("cx_en (x-deriv coeffs of energy):")
for f,e in zip(xds,cx_en): print(f"  d/d{f}: {sp.simplify(e.subs(r,1))}")
print("cs_en (s-deriv coeffs of energy):")
for f,e in zip(sds,cs_en): print(f"  d/d{f}: {sp.simplify(e.subs(r,1))}")
print("wrote linsys.pkl")
