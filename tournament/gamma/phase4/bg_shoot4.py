# Analyticity residual = sign of dV/dx evaluated just BEFORE the sonic locus (at fixed small
# det>0). For analytic passage nV->0 so V'=nV/det stays finite/changes sign; off-critical a2
# gives V'->+inf or -inf with definite sign. Bracket the sign flip.
import sympy as sp, pickle, numpy as np
from scipy.optimize import brentq
d=pickle.load(open("rhs_full.pkl","rb"))
A0,N0,o0,V0=sp.symbols('A0 N0 om0 V0')
f_om=sp.lambdify((A0,N0,o0,V0),sp.sympify(d['omx']),'numpy',cse=True)
f_V =sp.lambdify((A0,N0,o0,V0),sp.sympify(d['Vxx']),'numpy',cse=True)
def Axf(A,N,om,V): return A*(1-A+(2*om/(1-V*V))*(1+V*V/3.0))
def Nxf(A,N,om,V): return N*(-2+A-2.0*om/3.0)
def slp(Y):
    N,A,om,V=Y; return np.array([Nxf(A,N,om,V),Axf(A,N,om,V),float(f_om(A,N,om,V)),float(f_V(A,N,om,V))])
def Dloc(N,V): return 3*N*N*V*V-N*N+4*N*V-V*V+3
NC=1.5
def run(a2, z0=np.exp(-10), xmax=3.0, h=2e-4, Dthr=0.05):
    v1=-1.5/NC; w=1.5*a2
    Y=np.array([NC/z0,1+a2*z0*z0,w*z0*z0,v1*z0]); x=np.log(z0)
    for i in range(int((xmax-x)/h)):
        N,A,om,V=Y; D=Dloc(N,V)
        # once we're near the sonic locus (D small +) on the V>0 approach, record V'
        if 0<D<Dthr and V>0.05:
            Vp=float(f_V(A,N,om,V))
            return np.sign(Vp), x, Y, Vp
        k1=slp(Y);k2=slp(Y+0.5*h*k1);k3=slp(Y+0.5*h*k2);k4=slp(Y+h*k3)
        Yn=Y+(h/6)*(k1+2*k2+2*k3+k4)
        if not np.all(np.isfinite(Yn)) or abs(Yn[3])>=1 or Yn[1]<=0:
            return np.sign(k1[3]), x, Y, k1[3]
        Y=Yn;x+=h
    return 0.0,x,Y,0.0
print(f"NC={NC} ingoing. residual = sign(V') at D=0.05 approach:")
prev=None;brack=None
xs=np.linspace(0.02,1.0,25)
for a2 in xs:
    s,x,Y,vp=run(a2)
    print(f"  a2={a2:.4f}: signV'={s:+.0f} V'={vp:+.3e} @x={x:+.3f} N={Y[0]:.3f} A={Y[1]:.4f} om={Y[2]:.5f} V={Y[3]:.4f}")
    if prev and prev[1]*s<0 and prev[1] and s: brack=(prev[0],a2)
    prev=(a2,s)
if brack:
    print("bracket",brack)
    ac=brentq(lambda a: run(a)[0],brack[0],brack[1],xtol=1e-11,maxiter=100)
    s,x,Y,vp=run(ac,Dthr=0.02)
    print(f"CRITICAL a2*={ac:.9f}: near-sonic x={x:.4f} N={Y[0]:.6f} A={Y[1]:.6f} om={Y[2]:.6f} V={Y[3]:.6f} V'={vp:.4e}")
    # save the critical solution profile
    np.save('a2crit.npy', np.array([ac,NC]))
else:
    print("no sign flip found in scan")
