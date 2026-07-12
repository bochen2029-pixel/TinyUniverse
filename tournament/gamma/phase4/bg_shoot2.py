# Lightweight background shoot: residual = signed behavior of dV/dx as the trajectory nears the
# sonic point. For a2 below critical the solution veers one way (V'->+inf), above the other
# (V'->-inf); the analytic solution threads between (finite). Bracket the sign flip in a2.
import sympy as sp, pickle, numpy as np
from scipy.optimize import brentq
d=pickle.load(open("rhs_full.pkl","rb"))
A0,N0,o0,V0=sp.symbols('A0 N0 om0 V0')
f_om=sp.lambdify((A0,N0,o0,V0),sp.sympify(d['omx']),'numpy',cse=True)
f_V =sp.lambdify((A0,N0,o0,V0),sp.sympify(d['Vxx']),'numpy',cse=True)
def Axf(A,N,om,V): return A*(1-A+(2*om/(1-V*V))*(1+V*V/3.0))
def Nxf(A,N,om,V): return N*(-2+A-2.0*om/3.0)
def slopes(Y):
    N,A,om,V=Y; return np.array([Nxf(A,N,om,V),Axf(A,N,om,V),float(f_om(A,N,om,V)),float(f_V(A,N,om,V))])
def Dloc(N,V): return 3*N*N*V*V-N*N+4*N*V-V*V+3
NC=1.5
def run(a2, z0=np.exp(-13), xmax=6.0, h=5e-5):
    v1=1.5/NC; w=1.5*a2
    Y=np.array([NC/z0,1+a2*z0*z0,w*z0*z0,v1*z0]); x=np.log(z0)
    Dprev=Dloc(Y[0],Y[3]); lastVp=0.0
    for i in range(int((xmax-x)/h)):
        k1=slopes(Y);k2=slopes(Y+0.5*h*k1);k3=slopes(Y+0.5*h*k2);k4=slopes(Y+h*k3)
        Yn=Y+(h/6)*(k1+2*k2+2*k3+k4)
        if not np.all(np.isfinite(Yn)) or abs(Yn[3])>=1 or Yn[1]<=0:
            # blew up: sign of last V' indicates side
            return np.sign(lastVp), x, Y
        D=Dloc(Yn[0],Yn[3])
        if D*Dprev<0 and Yn[3]>0.05:   # crossed sonic locus with V>0 -> analytic passage!
            return 0.0, x, Yn   # threaded through
        lastVp=k1[3]; Dprev=D; Y=Yn; x+=h
    return 0.0, x, Y
print(f"NC={NC}: bracket a2 for analytic sonic passage (residual = side of blow-up):")
prev=None
brack=None
for a2 in np.linspace(0.05,0.7,27):
    s,x,Y=run(a2)
    print(f"  a2={a2:.4f}: resid={s:+.0f} stop@x={x:+.3f} N={Y[0]:.3f} A={Y[1]:.4f} om={Y[2]:.5f} V={Y[3]:.4f}")
    if prev is not None and prev[1]*s<0 and prev[1]!=0 and s!=0:
        brack=(prev[0],a2)
    prev=(a2,s)
if brack:
    print("bracket:",brack)
    ac=brentq(lambda a: run(a)[0], brack[0],brack[1],xtol=1e-9,maxiter=80)
    s,x,Y=run(ac)
    print(f"CRITICAL a2*={ac:.8f}: sonic x={x:.4f} N={Y[0]:.5f} A={Y[1]:.5f} om={Y[2]:.5f} V={Y[3]:.5f}")
