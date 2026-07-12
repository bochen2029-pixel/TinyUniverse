# hka_pert_symbols2.py — CORRECTED perturbation operator L(x;kappa).
#
# The prior assembly (hka_pert_symbols) used G,H (5.9,5.10) as the Abar,Nbar evolution rows and FAILED
# the gauge-mode exactness test. The auxiliary identity (5.14) — VERIFIED exact on the gauge mode —
# gives the true Abar row (WITH a -kappa Abar_p term). By analogy (and matching H1,H3), the Nbar row
# is the linearization of 4.1b with a -kappa Nbar_p term. The (ombar,V) rows follow (5.5)-(5.8) with
# the kappa Ms coupling; we ALSO include the -kappa shift on the diagonal there (from the s-derivative
# on e^{kappa s}) and VALIDATE the whole thing on the gauge mode (must be exact for all kappa).
#
# Rows (Psi=(Abar,Nbar,ombar,V)_p):
#  Abar_p' = -kappa Abar_p - (2 g N V om/oV2)(Nbar_p+ombar_p) - (2 g N om(1+V^2)/oV2^2) V_p     [5.14]
#  Nbar_p' = -kappa Nbar_p + A Abar_p + (g-2) om ombar_p                                        [lin 4.1b]
#  [[Ax,Bx],[Cx,Dx]] (ombar_p', V_p') = (E.Psi ; F.Psi) + kappa (As,Bs;Cs,Ds)(ombar_p,V_p)      [5.5-5.8, 5.13]
#     -> (ombar_p',V_p') = Minv[ E.Psi + kappa Ms (ombar_p,V_p) ]
# where E,F use ombar_x = om'/om (log deriv) and V_x (background slopes).
#
# We build L and lambdify; validate on the gauge mode in __main__.
import sympy as sp, pickle, os

def build():
    g=sp.Rational(4,3)
    A,N,om,V=sp.symbols('A N om V', real=True)
    omx,Vx=sp.symbols('om_x V_x', real=True)     # ombar_x := om'/om ; V_x (background)
    k=sp.symbols('kappa'); oV2=1-V**2
    # (5.6)
    Ax=1+N*V; Bx=g*(N+V)/oV2; Cx=(g-1)*(N+V); Dx=g*(1+N*V)/oV2
    # (5.5)
    As=sp.Integer(1); Bs=g*V/oV2; Cs=(g-1)*V; Ds=g/oV2
    # (5.7) E1..E4
    E1=-((g+2)/2)*A*N*V
    E2=((6-3*g)/2)*N*V-((2+g)/2)*A*N*V+(2-g)*om*N*V-N*V*omx-g*N*Vx/oV2
    E3=(2-g)*om*N*V
    E4=((6-3*g)/2)*N-((2+g)/2)*A*N+(2-g)*om*N-N*omx-g*(1+2*N*V+V**2)*Vx/oV2**2
    # (5.8) F1..F4
    F1=((2-3*g)/2)*A*N
    F2=(2-g)*(g-1)*om*N+((7*g-6)/2)*N+((2-3*g)/2)*A*N-(g-1)*N*omx-g*N*V*Vx/oV2
    F3=(2-g)*(g-1)*om*N
    F4=-(g-1)*omx-g*(N+2*V+N*V**2)*Vx/oV2**2
    Mx=sp.Matrix([[Ax,Bx],[Cx,Dx]]); Minv=Mx.inv()
    Ms=sp.Matrix([[As,Bs],[Cs,Ds]])
    Evec=sp.Matrix([[E1,E2,E3,E4],[F1,F2,F3,F4]])
    sel=sp.Matrix([[0,0,1,0],[0,0,0,1]])
    L_ov=Minv*Evec + k*(Minv*Ms*sel)     # 2x4, the (ombar',V') rows -- NOTE: no extra -k here; the
    # kappa enters ONLY via +kappa Ms(ombar,V) per (5.13). We test whether an extra -k Psi_ov is needed.
    # Abar row (5.14):
    rowA=[-k, -(2*g*N*V*om/oV2), -(2*g*N*V*om/oV2), -(2*g*N*om*(1+V**2)/oV2**2)]
    # Nbar row (lin 4.1b):
    rowN=[A, -k, (g-2)*om, 0]
    L=sp.Matrix([rowA, rowN, list(L_ov.row(0)), list(L_ov.row(1))])
    L=sp.simplify(L)
    pickle.dump(dict(L=sp.srepr(L)), open("hka_pert_L2.pkl","wb"))
    return L,(A,N,om,V,omx,Vx,k)

if __name__=="__main__":
    import numpy as np, math, time
    t0=time.time(); L,syms=build(); A,N,om,V,omx,Vx,k=syms
    print(f"assembled corrected L in {time.time()-t0:.1f}s; affine in kappa:",
          not any(sp.diff(L,k)[i,j].has(k) for i in range(4) for j in range(4)))
    # numeric lambdify and validate on the gauge mode
    fL=[[sp.lambdify((A,N,om,V,omx,Vx,k),L[i,j],'numpy') for j in range(4)] for i in range(4)]
    import hka_beta4 as B; B.bg(); B.bg_path(); G=4.0/3.0
    def bgfld(x):
        Ab,Nb,ob,Vb,ox,Vxb=B.bg_state(x); return (Ab,Nb,ob,Vb,ox,Vxb)
    def Lnum(x,kap):
        f=bgfld(x); return np.array([[complex(fL[i][j](*f,complex(kap))) for j in range(4)] for i in range(4)])
    def gauge_vec(x,kbar):
        Ab,Nb,ob,Vb,ox,Vxb=B.bg_state(x); oV2=1-Vb*Vb
        Ax=Ab*(1-Ab+2*ob*(1+(1/3.0)*Vb*Vb)/oV2); Nx=Nb*(-2+Ab-(2-G)*ob)
        return np.array([Ax/Ab, Nx/Nb+kbar, ox, Vxb],complex)
    print("\nGAUGE-MODE EXACTNESS TEST (d/dx Psi_gauge - L Psi_gauge; must be ~0):")
    for kbar in [0.35699,1.0,2.81]:
        errs=[]
        for x in [-8.0,-4.0,-1.0,-0.3]:
            Psi=gauge_vec(x,kbar); dP=(gauge_vec(x+1e-5,kbar)-gauge_vec(x-1e-5,kbar))/2e-5
            res=dP-Lnum(x,kbar).dot(Psi); errs.append(np.linalg.norm(res)/np.linalg.norm(dP))
        print(f"  kbar={kbar}: rel residuals at x=-8,-4,-1,-0.3 = {[f'{e:.1e}' for e in errs]}")
