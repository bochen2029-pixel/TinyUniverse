# Integrate inward from the sonic point with scipy solve_ivp (stiff, adaptive) and examine
# whether the trajectory approaches (A=1,V=0). Scan V0 and eigendirection. Use the eigendir
# whose inward orientation heads toward decreasing A / V->0.
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve

def branch1_N(V): return (-2*V+np.sqrt(3)*(V*V-1))/(3*V*V-1)
def metric_slopes(N,A,om,V):
    Ap=A*(1-A+(2*om/(1-V*V))*(1+V*V/3.0)); Np=N*(-2+A-2.0*om/3.0); return Ap,Np
def rhs(x,Y):
    N,A,om,V=Y; Ap,Np=metric_slopes(N,A,om,V); omV2=1-V*V
    m11=(1+N*V);m12=4*(N+V)/(3*omV2);m21=(4*V+N+3*N*V*V);m22=4*(1+V*V+2*N*V)/omV2
    r1=-(-N*V*Ap/(3*A)+4*V*Np/3.0+2*N*(1+4*om/(9*omV2)))
    r2=-( N*omV2*Ap/A+4*(1+V*V)*Np+2*N*(1+3*V*V))
    det=m11*m22-m12*m21
    return [Np,Ap,(r1*m22-m12*r2)/det*om,(m11*r2-r1*m21)/det]
def sonic_state(V0):
    N0=branch1_N(V0)
    def res(p):
        A0,om0=p; Ap,Np=metric_slopes(N0,A0,om0,V0); omV2=1-V0*V0
        m11=(1+N0*V0);m12=4*(N0+V0)/(3*omV2);m21=(4*V0+N0+3*N0*V0*V0);m22=4*(1+V0*V0+2*N0*V0)/omV2
        r1=-(-N0*V0*Ap/(3*A0)+4*V0*Np/3.0+2*N0*(1+4*om0/(9*omV2)));r2=-(N0*omV2*Ap/A0+4*(1+V0*V0)*Np+2*N0*(1+3*V0*V0))
        return [r1*m22-m12*r2, m11*r2-r1*m21]
    for g in [[1.1,0.01],[1.15,0.011],[1.0,0.02],[1.2,0.05],[1.05,0.005],[1.3,0.02]]:
        s=fsolve(res,g,full_output=True)
        if s[2]==1 and abs(res(s[0])[0])+abs(res(s[0])[1])<1e-11 and s[0][0]>0 and s[0][1]>0:
            return N0,s[0][0],s[0][1]
    return N0,None,None

# Jacobian eigen-directions at sonic to step off correctly
def eigdirs(V0,N0,A0,om0):
    Y0=np.array([N0,A0,om0,V0]); eps=1e-6; J=np.zeros((4,4))
    # rhs is singular at sonic; instead use the local expansion: near sonic the analytic
    # slope dY/dx is finite. Estimate it by evaluating rhs slightly off using the metric slopes
    # exact and L'Hopital for om,V. Simpler: numerically get dY/dx just off x=0 by stepping the
    # KNOWN metric slopes and reading om',V' as the finite limit from both sides.
    hs=1e-5
    outs=[]
    for sgn in [+1,-1]:
        # perturb along det gradient a hair then eval rhs
        Ntest=N0; # move V a hair to get off det=0
        Y=np.array([N0,A0,om0,V0+sgn*hs])
        try:
            d=rhs(0,Y); outs.append((sgn,np.array(d)))
        except Exception: pass
    return outs

print("inward integration via solve_ivp (LSODA), look for A->1,V->0:")
for V0 in np.linspace(-0.45,0.45,19):
    N0,A0,om0=sonic_state(V0)
    if A0 is None: continue
    Np0,Ap0=metric_slopes(N0,A0,om0,V0)
    # step off sonic toward center using metric slopes + numeric om',V' just inside
    Yin=np.array([N0,A0,om0,V0]); h0=5e-3
    dj=rhs(0,np.array([N0,A0,om0,V0-1e-4]))  # slopes just inside
    Y=Yin - h0*np.array([dj[0],dj[1],dj[2],dj[3]])*0 + np.array([-Np0*h0,-Ap0*h0,-dj[2]*h0,-dj[3]*h0])
    sol=solve_ivp(rhs,[-h0,-9],Y,method='LSODA',rtol=1e-9,atol=1e-12,dense_output=False,max_step=0.05)
    if sol.success or sol.t[-1]<-1:
        Yf=sol.y[:,-1]
        print(f"V0={V0:+.3f}: reached x={sol.t[-1]:+.2f}  N={Yf[0]:.2f} A={Yf[1]:.4f} om={Yf[2]:.4f} V={Yf[3]:.4f}")
    else:
        print(f"V0={V0:+.3f}: fail at x={sol.t[-1]:.2f}")
