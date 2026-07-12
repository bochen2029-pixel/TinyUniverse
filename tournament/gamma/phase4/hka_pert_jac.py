# hka_pert_jac.py — perturbation operator from the background FLOW JACOBIAN (first principles) + the
# kappa s-coupling, validated on the exact gauge mode.
#
# Background log-flow: d/dx (lnA,lnN,lnom,V) = R(A,N,om,V):
#   R1=(lnA)'=1-A+2om(1+(g-1)V^2)/(1-V^2)                 [4.1a]
#   R2=(lnN)'=-2+A-(2-g)om                                [4.1b]
#   R3=(lnom)'=om'/om   (om' from fluid pair)             [4.1d,e]
#   R4=V'=V'            (from fluid pair)
# L0 = d R / d(lnA,lnN,lnom,V) is the kappa^0 perturbation operator (Jacobian; chain rule:
#   dR/d(lnA)=A dR/dA, etc.).
# kappa-coupling K (from 5.14 + 5.13):
#   Abar,Nbar rows: -kappa on the diagonal.
#   (ombar,V) rows: +kappa Minv Ms (ombar,V), Minv=[[Ax,Bx],[Cx,Dx]]^{-1}, Ms=[[As,Bs],[Cs,Ds]].
# L(x;kappa) = L0(x) + kappa K(x). VALIDATE: gauge mode Psi_g=((lnA)',(lnN)'+kbar,(lnom)',V') must
# solve Psi_g' = L(x;kbar) Psi_g exactly for all kbar.
import numpy as np, sympy as sp, pickle, math

def build():
    g=sp.Rational(4,3)
    lA,lN,lo,V=sp.symbols('lA lN lo V', real=True)     # log variables + V
    A=sp.exp(lA); N=sp.exp(lN); om=sp.exp(lo); oV2=1-V**2
    # background slopes (resolved) as functions of the log vars
    Ax=A*(1-A+2*om*(1+(g-1)*V**2)/oV2)
    Nx=N*(-2+A-(2-g)*om)
    RHS_d=(3*(2-g)/2)*N*V-((2+g)/2)*A*N*V+(2-g)*N*V*om
    RHS_e=(2-g)*(g-1)*N*om+((7*g-6)/2)*N+((2-3*g)/2)*A*N
    a=(1+N*V)/om; b=g*(N+V)/oV2; c=(g-1)*(N+V)/om; d=g*(1+N*V)/oV2
    det=a*d-b*c
    omx=(d*RHS_d-b*RHS_e)/det         # om'
    Vx=(a*RHS_e-c*RHS_d)/det          # V'
    R1=1-A+2*om*(1+(g-1)*V**2)/oV2    # (lnA)'
    R2=-2+A-(2-g)*om                  # (lnN)'
    R3=omx/om                         # (lnom)'
    R4=Vx                             # V'
    R=sp.Matrix([R1,R2,R3,R4])
    L0=R.jacobian([lA,lN,lo,V])       # kappa^0 Jacobian
    # kappa coupling
    Ax_=1+N*V; Bx_=g*(N+V)/oV2; Cx_=(g-1)*(N+V); Dx_=g*(1+N*V)/oV2
    As=sp.Integer(1); Bs=g*V/oV2; Cs=(g-1)*V; Ds=g/oV2
    Mx=sp.Matrix([[Ax_,Bx_],[Cx_,Dx_]]); Minv=Mx.inv(); Ms=sp.Matrix([[As,Bs],[Cs,Ds]])
    MinvMs=Minv*Ms
    K=sp.zeros(4,4)
    K[0,0]=-1; K[1,1]=-1
    # (ombar,V) rows get +Minv Ms acting on (ombar,V)=components 2,3
    K[2,2]=MinvMs[0,0]; K[2,3]=MinvMs[0,1]
    K[3,2]=MinvMs[1,0]; K[3,3]=MinvMs[1,1]
    L0=sp.simplify(L0); K=sp.simplify(K)
    pickle.dump(dict(L0=sp.srepr(L0),K=sp.srepr(K)), open("hka_pert_jac.pkl","wb"))
    return L0,K,(lA,lN,lo,V)

if __name__=="__main__":
    import time
    t0=time.time(); L0,K,syms=build(); lA,lN,lo,V=syms
    print(f"built L0 (Jacobian) + K in {time.time()-t0:.1f}s")
    fL0=[[sp.lambdify((lA,lN,lo,V),L0[i,j],'numpy') for j in range(4)] for i in range(4)]
    fK=[[sp.lambdify((lA,lN,lo,V),K[i,j],'numpy') for j in range(4)] for i in range(4)]
    import hka_beta4 as B; B.bg(); B.bg_path(); G=4.0/3.0
    def logvars(x):
        A,N,om,Vv,ox,Vx=B.bg_state(x); return (math.log(A),math.log(N),math.log(om),Vv)
    def Lnum(x,kap):
        lv=logvars(x)
        L0n=np.array([[complex(fL0[i][j](*lv)) for j in range(4)] for i in range(4)])
        Kn=np.array([[complex(fK[i][j](*lv)) for j in range(4)] for i in range(4)])
        return L0n+kap*Kn
    def gauge_vec(x,kbar):
        A,N,om,Vv,ox,Vx=B.bg_state(x); oV2=1-Vv*Vv
        Ax=A*(1-A+2*om*(1+(1/3.0)*Vv*Vv)/oV2); Nx=N*(-2+A-(2-G)*om)
        return np.array([Ax/A, Nx/N+kbar, ox, Vx],complex)
    print("\nGAUGE-MODE EXACTNESS (d/dx Psi_g - L Psi_g, relative; must be ~0):")
    for kbar in [0.35699,1.0,2.81]:
        errs=[]
        for x in [-8.0,-4.0,-1.0,-0.3]:
            Psi=gauge_vec(x,kbar); dP=(gauge_vec(x+1e-5,kbar)-gauge_vec(x-1e-5,kbar))/2e-5
            res=dP-Lnum(x,complex(kbar)).dot(Psi); errs.append(np.linalg.norm(res)/max(np.linalg.norm(dP),1e-12))
        print(f"  kbar={kbar}: {[f'{e:.1e}' for e in errs]}")
