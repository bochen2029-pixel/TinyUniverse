# Shoot the background: find the center parameters so the solution passes analytically through
# the sonic point. Residual = how far the state at the sonic crossing is from the sonic-regular
# manifold. Reduce params: fix nc (absorbs an x-translation / scale), shoot on a2 (density).
import sympy as sp, pickle, numpy as np
from scipy.optimize import brentq, fsolve
d=pickle.load(open("rhs_full.pkl","rb"))
A0,N0,o0,V0=sp.symbols('A0 N0 om0 V0')
omx=sp.sympify(d['omx']); Vxx=sp.sympify(d['Vxx'])
f_om=sp.lambdify((A0,N0,o0,V0),omx,'numpy',cse=True)
f_V =sp.lambdify((A0,N0,o0,V0),Vxx,'numpy',cse=True)
# sonic regularity numerators from the full system:
en_s,mo_s=pickle.load(open("fluid_full.pkl","rb"))
en=sp.sympify(en_s); mo=sp.sympify(mo_s)
Ax,Nx,ox,Vx,Es=sp.symbols('A_x N_x om_x V_x E')
Axv=A0*(1-A0+(2*o0/(1-V0**2))*(1+V0**2/3)); Nxv=N0*(-2+A0-sp.Rational(2,3)*o0)
e1=sp.expand(en.subs({Ax:Axv,Nx:Nxv,Es:1})); e2=sp.expand(mo.subs({Ax:Axv,Nx:Nxv,Es:1}))
a11=e1.coeff(ox);a12=e1.coeff(Vx);a21=e2.coeff(ox);a22=e2.coeff(Vx)
b1=-(e1.subs({ox:0,Vx:0}));b2=-(e2.subs({ox:0,Vx:0}))
nV=sp.simplify(sp.numer(sp.together(a11*b2-b1*a21)))
f_nV=sp.lambdify((A0,N0,o0,V0),nV,'numpy',cse=True)
def Axf(A,N,om,V): return A*(1-A+(2*om/(1-V*V))*(1+V*V/3.0))
def Nxf(A,N,om,V): return N*(-2+A-2.0*om/3.0)
def slopes(Y):
    N,A,om,V=Y; return np.array([Nxf(A,N,om,V),Axf(A,N,om,V),float(f_om(A,N,om,V)),float(f_V(A,N,om,V))])
def Dloc(N,V): return 3*N*N*V*V-N*N+4*N*V-V*V+3

NC=1.5
def integrate(a2, z0=np.exp(-13), xmax=6.0, h=1e-4):
    v1=1.5/NC; w=1.5*a2
    N=NC/z0;A=1+a2*z0*z0;om=w*z0*z0;V=v1*z0;x=np.log(z0)
    Yp=np.array([N,A,om,V])
    for i in range(int((xmax-x)/h)):
        Y=Yp
        k1=slopes(Y);k2=slopes(Y+0.5*h*k1);k3=slopes(Y+0.5*h*k2);k4=slopes(Y+h*k3)
        Yn=Y+(h/6)*(k1+2*k2+2*k3+k4);x+=h
        if not np.all(np.isfinite(Yn)) or abs(Yn[3])>=1 or Yn[1]<=0:
            return ('blow',x,Yp)
        D0=Dloc(Y[0],Y[3]);D1=Dloc(Yn[0],Yn[3])
        if D0*D1<0 and Yn[3]>0.05:
            fr=D0/(D0-D1); Ys=Y+fr*(Yn-Y)
            return ('sonic',x-h+fr*h,Ys)
        Yp=Yn
    return ('end',x,Yp)

# residual: nV at the sonic crossing (analytic passage <=> nV=0 there since det=0 there)
def resid(a2):
    tag,x,Y=integrate(a2)
    if tag!='sonic': return np.nan
    return float(f_nV(Y[1],Y[0],Y[2],Y[3]))

print(f"NC={NC}. Shoot a2 so nV(sonic)=0 (analytic passage):")
a2s=np.linspace(0.02,0.8,30)
prev=None
for a2 in a2s:
    r=resid(a2)
    tag,x,Y=integrate(a2)
    mark=""
    if prev is not None and np.isfinite(r) and np.isfinite(prev[1]) and r*prev[1]<0: mark=" <== SIGN CHANGE"
    print(f"  a2={a2:.3f}: sonic@x={x:+.3f} V={Y[3]:.4f} A={Y[1]:.4f} om={Y[2]:.5f} nV={r:+.4e}{mark}")
    prev=(a2,r)
