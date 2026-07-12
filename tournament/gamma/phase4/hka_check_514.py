# Check the auxiliary identity (5.14) on the gauge mode, and figure out the correct Abar evolution.
# (5.14): kappa Abar_p + Abar_p' = -(2 g N V om/(1-V^2))(Nbar_p+ombar_p) - (2 g N om(1+V^2)/(1-V^2)^2) V_p
# Gauge mode: Abar_p=(lnA)'_x, Nbar_p=(lnN)'_x + kbar, ombar_p=(lnom)'_x, V_p=V'_x, kappa=kbar.
import numpy as np, math
import hka_beta4 as B
G=4.0/3.0
B.bg(); B.bg_path()

def gauge_vec(x, kbar):
    A,N,om,V,ombar_x,Vx=B.bg_state(x); oV2=1-V*V
    Ax=A*(1-A+2*om*(1+(1/3.0)*V*V)/oV2); Nx=N*(-2+A-(2-G)*om)
    return np.array([Ax/A, Nx/N+kbar, ombar_x, Vx],complex), (A,N,om,V)

def rhs_514(x, kbar):
    Psi,(A,N,om,V)=gauge_vec(x,kbar); oV2=1-V*V
    Abar,Nbar,ombar,Vp=Psi
    return -(2*G*N*V*om/oV2)*(Nbar+ombar) - (2*G*N*om*(1+V*V)/oV2**2)*Vp

def lhs_514(x, kbar, dx=1e-5):
    Psi,_=gauge_vec(x,kbar); Psip=(gauge_vec(x+dx,kbar)[0]-gauge_vec(x-dx,kbar)[0])/(2*dx)
    return kbar*Psi[0]+Psip[0]

if __name__=="__main__":
    print("(5.14) on the gauge mode: LHS = kappa Abar_p + Abar_p' vs RHS. Should match if gauge mode &")
    print("(5.14) are consistent (validates the constraint 5.14 and the gauge-mode form).")
    for kbar in [0.35699,1.0,2.81]:
        print(f"\n kbar={kbar}:")
        for x in [-8.0,-4.0,-1.0,-0.3]:
            l=lhs_514(x,kbar); r=rhs_514(x,kbar)
            print(f"   x={x:+.2f}: LHS={l:+.5e}  RHS={r:+.5e}  diff={abs(l-r):.2e}")
