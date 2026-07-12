# Prototype background shoot using MY covariantly-derived fluid ODEs (which give correct
# center regularity dV/dx->0). Export omx_mine, Vx_mine as fast numpy, integrate from center.
import sympy as sp, pickle, numpy as np
en_s,mo_s=pickle.load(open("fluid_manual.pkl","rb"))
en=sp.sympify(en_s); mo=sp.sympify(mo_s)
A0,N0,o0,V0,Ax,Nx,ox,Vx,Es=sp.symbols('A0 N0 om0 V0 A_x N_x om_x V_x E')
Axv=A0*(1-A0+(2*o0/(1-V0**2))*(1+V0**2/3))
Nxv=N0*(-2+A0-sp.Rational(2,3)*o0)
sol=sp.solve([en.subs({Ax:Axv,Nx:Nxv}),mo.subs({Ax:Axv,Nx:Nxv})],[ox,Vx],dict=True)[0]
omx=sp.simplify(sol[ox]); Vxx=sp.simplify(sol[Vx])
f_om=sp.lambdify((A0,N0,o0,V0),omx,'numpy')
f_V =sp.lambdify((A0,N0,o0,V0),Vxx,'numpy')
def Axf(A,N,om,V): return A*(1-A+(2*om/(1-V*V))*(1+V*V/3.0))
def Nxf(A,N,om,V): return N*(-2+A-2.0*om/3.0)
def slopes(Y):
    N,A,om,V=Y
    return np.array([Nxf(A,N,om,V),Axf(A,N,om,V),f_om(A,N,om,V),f_V(A,N,om,V)])

# Center seed: A=1+a2 z^2, V=v1 z, om=w z^2, N=nc/z with z=e^{x0}. Determine leading coeffs
# by consistency with MY slopes at small z. First find v1,w,a2 relations numerically.
# Probe: near center dV/dx->0, dom/dx->0; sub-leading gives the manifold. Just seed with a
# guess and let a short pre-integration relax? Better: fit seed by requiring slopes ~ ansatz deriv.
# Numerically: pick z0=e^{-6}. Unknown center params: nc, v1, w, a2. 4 unknowns.
# Consistency eqs at z0: dA/dx = Axf must equal d/dx(1+a2 z^2)=2 a2 z^2; similarly others.
from scipy.optimize import fsolve
def seed_res(p, z):
    nc,v1,w,a2=p
    N=nc/z; A=1+a2*z*z; om=w*z*z; V=v1*z
    dN,dA,dom,dV=slopes([N,A,om,V])
    # ansatz x-derivatives: dA/dx=z d/dz(1+a2 z^2)=2 a2 z^2; dV/dx=v1 z; dom/dx=2w z^2; dN/dx=-nc/z
    return [dA-2*a2*z*z, dV-v1*z, dom-2*w*z*z, dN-(-nc/z)]
z0=np.exp(-6)
guesses=[[1.5,-0.5,0.1,0.1],[1.7,-0.3,0.05,0.05],[2,-0.6,0.2,0.1],[1.2,-0.2,0.02,0.02],[1.0,-0.1,0.01,0.01]]
seeds=[]
for g in guesses:
    s=fsolve(seed_res,g,args=(z0,),full_output=True)
    if s[2]==1 and max(abs(np.array(seed_res(s[0],z0))))<1e-8:
        seeds.append(tuple(s[0]))
uniq=[]
for s in seeds:
    if not any(np.allclose(s,q,atol=1e-4) for q in uniq): uniq.append(s)
print("center seeds (nc,v1,w,a2) at z0=e^-6:")
for s in uniq: print("  ",[round(x,5) for x in s])

# integrate outward from center to find sonic point (det of KHA matrix ~ D=0 crossing) & beyond.
def Dloc(N,V): return 3*N*N*V*V-N*N+4*N*V-V*V+3
def integrate_out(seed, xmax=6.0, h=1e-3):
    nc,v1,w,a2=seed; z=z0
    N=nc/z;A=1+a2*z*z;om=w*z*z;V=v1*z;x=np.log(z0)
    traj=[(x,N,A,om,V)]
    n=int((xmax-x)/h)
    for i in range(n):
        Y=np.array([N,A,om,V])
        k1=slopes(Y);k2=slopes(Y+0.5*h*k1);k3=slopes(Y+0.5*h*k2);k4=slopes(Y+h*k3)
        Y=Y+(h/6)*(k1+2*k2+2*k3+k4);x+=h;N,A,om,V=Y
        traj.append((x,N,A,om,V))
        if not np.all(np.isfinite(Y)) or abs(V)>=0.9999 or A<=0 or A>100: break
    return traj
for s in uniq[:3]:
    tr=integrate_out(s)
    print(f"\nseed {[round(x,4) for x in s]}: {len(tr)} steps, reached x={tr[-1][0]:.2f}")
    # print sonic crossing (D changes sign) and endpoint
    for k in range(1,len(tr)):
        D0=Dloc(tr[k-1][1],tr[k-1][4]); D1=Dloc(tr[k][1],tr[k][4])
        if D0*D1<0:
            print(f"   SONIC at x={tr[k][0]:.4f}: N={tr[k][1]:.4f} A={tr[k][2]:.4f} om={tr[k][3]:.4f} V={tr[k][4]:.4f}")
            break
    print(f"   endpoint: x={tr[-1][0]:.2f} N={tr[-1][1]:.3f} A={tr[-1][2]:.4f} om={tr[-1][3]:.5f} V={tr[-1][4]:.4f}")
