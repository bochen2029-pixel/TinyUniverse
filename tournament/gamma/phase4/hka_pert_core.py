# hka_pert_core.py — numeric perturbation operator L(x;kappa) and background-field provider along
# the EC solution, for the Stage-B eigenvalue shoot.
#
# L(x;kappa) = P(fields) + kappa Q(fields), from hka_pert_symbols.build() (5.5-5.13). The background
# fields are (A,N,om,V) and the LOG-derivative ombar_x = om'/om and V_x, evaluated along the EC
# background. We lambdify P_ij, Q_ij entrywise (functions of A,N,om,V,ombar_x,V_x) for fast numeric L.
#
# Modes:
#  - CENTER (x->-inf): perturbation regular. md line 104: x->-inf modes {1,1,e^{-2x},e^{-x}(...)},
#    e^{-x} excluded by (5.15); the GROWING e^{-2x} mode must be killed by choosing kappa. Regular
#    (bounded) center perturbation: Abar_p=0 (md line 102-103).
#  - SONIC (x=0-ish): regular singular point (AxDx-BxCx=0). Analytic modes + one non-analytic
#    (indicial exponent 1-2kappa). Analyticity fixes the passage.
# This module provides the machinery; the shoot logic lives in hka_beta.py.
#
# Eq #s per HKA_beta_equations.md.

import numpy as np, sympy as sp, pickle, os, math
import hka_ec as E

G=4.0/3.0; GM1=G-1.0; CS=math.sqrt(GM1)

# ---- load/assemble symbolic L, lambdify P,Q ----
if not os.path.exists("hka_pert_L.pkl"):
    import hka_pert_symbols as _S; _S.build()
_data=pickle.load(open("hka_pert_L.pkl","rb"))
_A,_N,_om,_V=sp.symbols('A N om V', real=True)
_omx,_Vx=sp.symbols('om_x V_x', real=True)     # ombar_x = om'/om (log deriv), V_x
_k=sp.symbols('kappa')
_L=sp.sympify(_data['L'])
_P=_L.subs(_k,0)
_Q=sp.expand(_L-_P).applyfunc(lambda e: sp.expand(e).coeff(_k,1))
_args=(_A,_N,_om,_V,_omx,_Vx)
_fP=[[sp.lambdify(_args,_P[i,j],'numpy') for j in range(4)] for i in range(4)]
_fQ=[[sp.lambdify(_args,_Q[i,j],'numpy') for j in range(4)] for i in range(4)]

def bg_fields(N,om,V):
    """Return (A,N,om,V, ombar_x, V_x) from the reduced background at a state (N,om,V). ombar_x=om'/om."""
    A=E.A_of(N,om,V)
    omx,Vx=E._fluid_slopes(A,N,om,V)   # these are om', V' (x-derivatives)
    ombar_x=omx/om                     # log derivative
    return A,N,om,V,ombar_x,Vx

def Pmat(fld):
    A,N,om,V,ox,Vx=fld
    return np.array([[complex(_fP[i][j](A,N,om,V,ox,Vx)) for j in range(4)] for i in range(4)])
def Qmat(fld):
    A,N,om,V,ox,Vx=fld
    return np.array([[complex(_fQ[i][j](A,N,om,V,ox,Vx)) for j in range(4)] for i in range(4)])

def Lnum(fld,k):
    return Pmat(fld)+k*Qmat(fld)

# ---- aux identity (5.14) check ----
def aux_resid(fld, Psi, k):
    """kappa Abar_p + Abar_p' + (2 g N V om/(1-V^2))(Nbar_p+ombar_p) + (2 g N om(1+V^2)/(1-V^2)^2)V_p.
    Abar_p' from row 0 of L.Psi. Should be ~0 for a true solution (checks the L assembly)."""
    A,N,om,V,ox,Vx=fld; oV2=1-V*V
    L=Lnum(fld,k)
    Abar_p_prime=L[0].dot(Psi)
    Abar,Nbar,ombar,Vp=Psi
    return k*Abar+Abar_p_prime + (2*G*N*V*om/oV2)*(Nbar+ombar) + (2*G*N*om*(1+V*V)/oV2**2)*Vp

if __name__=="__main__":
    # examine L at the sonic point and near the center
    print("perturbation operator L(x;kappa) checks:")
    # sonic point fields (A0,N0,om0,V0)=(3/2,2/sqrt3,3/4,-1/sqrt3)
    N0=2/math.sqrt(3); om0=0.75; V0=-1/math.sqrt(3)
    fld=bg_fields(N0,om0,V0)
    print(f" sonic fields A,N,om,V,ombar_x,V_x = {tuple(round(float(np.real(x)),5) for x in fld)}")
    for k in [1.0, 2.81055255, 0.35699]:
        L=Lnum(fld,complex(k))
        ev=np.linalg.eigvals(L)
        print(f"  kappa={k}: eig(L at sonic)={np.sort_complex(ev).round(4)}")
    # near center: large N, small om,V
    z=math.exp(-6); N=1/z; om=0.375*z*z; V=(E.M_CENTER)*z
    fldc=bg_fields(N,om,V)
    print(f"\n center-ish (z=e^-6) fields A,N,om,V,ombar_x,V_x = {tuple(round(float(np.real(x)),4) for x in fldc)}")
    for k in [2.81055255]:
        L=Lnum(fldc,complex(k))
        ev=np.linalg.eigvals(L)
        print(f"  kappa={k}: eig(L near center)={np.sort_complex(ev).round(4)}  (indicial-like)")
