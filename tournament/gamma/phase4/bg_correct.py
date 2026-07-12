# Background CSS solver with the CORRECT full-4D derived ODEs. Seed the regular center, shoot
# the free parameter so the solution passes analytically through the sonic point.
import sympy as sp, pickle, numpy as np
from scipy.optimize import brentq
d=pickle.load(open("rhs_full.pkl","rb"))
A0,N0,o0,V0=sp.symbols('A0 N0 om0 V0')
f_om=sp.lambdify((A0,N0,o0,V0),sp.sympify(d['omx']),'numpy',cse=True)
f_V =sp.lambdify((A0,N0,o0,V0),sp.sympify(d['Vxx']),'numpy',cse=True)
def Axf(A,N,om,V): return A*(1-A+(2*om/(1-V*V))*(1+V*V/3.0))
def Nxf(A,N,om,V): return N*(-2+A-2.0*om/3.0)
def slopes(Y):
    N,A,om,V=Y
    return np.array([Nxf(A,N,om,V),Axf(A,N,om,V),float(f_om(A,N,om,V)),float(f_V(A,N,om,V))])
def Dloc(N,V): return 3*N*N*V*V-N*N+4*N*V-V*V+3

# Center series (z=e^x->0): A=1+a2 z^2, V=v1 z? Let's discover leading powers numerically from
# the correct system. dV/dx~1.5/N=1.5 z/nc near center -> V' ~ (1.5/nc) z => V=(1.5/nc) z^... wait
# V'=dV/dx=1.5 z/nc; integrate: V = (1.5/nc) z (since dV/dx=z dV/dz => dV/dz=1.5/nc => V=1.5 z/nc).
# So v1=1.5/nc. And dom/dx=2e-10-ish? that was with om=1e-10 forced. Let's get om scaling:
# dom/dx at center with om=om: proportional to om (homogeneous). So om~z^m free m. Need the actual.
# Just seed at small z0 with a family and integrate; the free parameter is nc (or the density).
def seed(nc, a2, z0):
    v1=1.5/nc
    w=1.5*a2   # guess omega=w z^2, w=1.5 a2 (from A-om relation A'=2a2 z^2, om couples)
    N=nc/z0; A=1+a2*z0*z0; om=w*z0*z0; V=v1*z0
    return np.array([N,A,om,V])

def integrate_out(nc,a2, z0=np.exp(-13), xmax=10.0, h=2e-4):
    Y=seed(nc,a2,z0); x=np.log(z0); traj=[(x,*Y)]; sonic=None
    for i in range(int((xmax-x)/h)):
        Yp=Y.copy()
        k1=slopes(Y);k2=slopes(Y+0.5*h*k1);k3=slopes(Y+0.5*h*k2);k4=slopes(Y+h*k3)
        Y=Y+(h/6)*(k1+2*k2+2*k3+k4);x+=h
        if not np.all(np.isfinite(Y)) or abs(Y[3])>=0.99999 or Y[1]<=0 or Y[1]>1e3: break
        D0=Dloc(Yp[0],Yp[3]);D1=Dloc(Y[0],Y[3])
        if sonic is None and D0*D1<0 and Y[3]>0.05:  # real sonic (not the trivial V~0 center)
            fr=D0/(D0-D1); sonic=(x-h+fr*h,*(Yp+fr*(Y-Yp)))
        traj.append((x,*Y))
    return traj, sonic

print("scan (nc,a2) center params -> trajectory to sonic point:")
print("(v1=1.5/nc is the center velocity slope; a2 sets density)")
for nc in [1.0,1.3,1.5,1.7,2.0]:
    for a2 in [0.05,0.1,0.2,0.4]:
        tr,son=integrate_out(nc,a2)
        e=tr[-1]
        if son:
            xs,N,A,om,V=son
            print(f"  nc={nc} a2={a2}: SONIC x={xs:.3f} N={N:.4f} A={A:.4f} om={om:.5f} V={V:.4f}  end x={e[0]:.2f}")
        else:
            print(f"  nc={nc} a2={a2}: no sonic; end x={e[0]:.2f} A={e[2]:.3f} om={e[3]:.4f} V={e[4]:.4f} ({len(tr)} steps)")
