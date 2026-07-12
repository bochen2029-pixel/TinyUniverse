# hka_pert_jac2.py — perturbation L = L0(flow Jacobian, log vars) + kappa K, with K determined from
# first principles: fluid rows get -kappa Minv Ms (the s-coupling with MINUS), Abar/Nbar rows get
# -kappa on diagonal. Validate on the exact gauge mode.
import numpy as np, sympy as sp, pickle, math

def build():
    g=sp.Rational(4,3)
    lA,lN,lo,V=sp.symbols('lA lN lo V', real=True)
    A=sp.exp(lA); N=sp.exp(lN); om=sp.exp(lo); oV2=1-V**2
    RHS_d=(3*(2-g)/2)*N*V-((2+g)/2)*A*N*V+(2-g)*N*V*om
    RHS_e=(2-g)*(g-1)*N*om+((7*g-6)/2)*N+((2-3*g)/2)*A*N
    a=(1+N*V)/om; b=g*(N+V)/oV2; c=(g-1)*(N+V)/om; d=g*(1+N*V)/oV2
    det=a*d-b*c
    omx=(d*RHS_d-b*RHS_e)/det; Vx=(a*RHS_e-c*RHS_d)/det
    R1=1-A+2*om*(1+(g-1)*V**2)/oV2; R2=-2+A-(2-g)*om; R3=omx/om; R4=Vx
    R=sp.Matrix([R1,R2,R3,R4]); L0=R.jacobian([lA,lN,lo,V])
    Ax_=1+N*V; Bx_=g*(N+V)/oV2; Cx_=(g-1)*(N+V); Dx_=g*(1+N*V)/oV2
    As=sp.Integer(1); Bs=g*V/oV2; Cs=(g-1)*V; Ds=g/oV2
    Minv=sp.Matrix([[Ax_,Bx_],[Cx_,Dx_]]).inv(); Ms=sp.Matrix([[As,Bs],[Cs,Ds]]); MM=Minv*Ms
    L0=sp.simplify(L0)
    pickle.dump(dict(L0=sp.srepr(L0),MM=sp.srepr(sp.simplify(MM))), open("hka_pert_jac2.pkl","wb"))
    return L0,MM,(lA,lN,lo,V)

if __name__=="__main__":
    import time
    t0=time.time(); L0,MM,syms=build(); lA,lN,lo,V=syms
    print(f"built in {time.time()-t0:.1f}s")
    fL0=[[sp.lambdify((lA,lN,lo,V),L0[i,j],'numpy') for j in range(4)] for i in range(4)]
    fMM=[[sp.lambdify((lA,lN,lo,V),MM[i,j],'numpy') for j in range(2)] for i in range(2)]
    import hka_beta4 as B; B.bg(); B.bg_path(); G=4.0/3.0
    def lv(x):
        A,N,om,V,ox,Vx=B.bg_state(x); return (math.log(A),math.log(N),math.log(om),V)
    def gv(x,kbar):
        A,N,om,V,ox,Vx=B.bg_state(x); oV2=1-V*V
        Ax=A*(1-A+2*om*(1+(1/3.0)*V*V)/oV2); Nx=N*(-2+A-(2-G)*om)
        return np.array([Ax/A,Nx/N+kbar,ox,Vx],complex)
    def gvp(x,kbar,dx=1e-6): return (gv(x+dx,kbar)-gv(x-dx,kbar))/(2*dx)
    # test several K choices for Abar/Nbar diagonal + fluid -MM
    def Lnum(x,kap, kA=-1.0, kN=-1.0, fluid_sign=-1.0):
        v=lv(x)
        L=np.array([[complex(fL0[i][j](*v)) for j in range(4)] for i in range(4)])
        MMn=np.array([[complex(fMM[i][j](*v)) for j in range(2)] for i in range(2)])
        K=np.zeros((4,4),complex); K[0,0]=kA; K[1,1]=kN
        K[2,2]=fluid_sign*MMn[0,0]; K[2,3]=fluid_sign*MMn[0,1]
        K[3,2]=fluid_sign*MMn[1,0]; K[3,3]=fluid_sign*MMn[1,1]
        return L+kap*K
    print("\nGauge-mode residuals for K variants (kA,kN,fluid_sign):")
    for (kA,kN,fs) in [(-1,-1,-1),(-1,-1,+1),(-1,0,-1),(0,-1,-1),(-1,-1,0)]:
        errs=[]
        for kbar in [1.0,2.81]:
            e=[]
            for x in [-3.0,-1.0,-0.3]:
                Psi=gv(x,kbar); dP=gvp(x,kbar); res=dP-Lnum(x,complex(kbar),kA,kN,fs).dot(Psi)
                e.append(np.linalg.norm(res)/max(np.linalg.norm(dP),1e-9))
            errs.append(max(e))
        print(f"  kA={kA} kN={kN} fluid_sign={fs}: max rel resid (kbar=1,2.81) = {errs[0]:.2e}, {errs[1]:.2e}")
