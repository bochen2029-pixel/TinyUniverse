# Correct background CSS solver using the covariantly-derived Euler ODEs (proper velocity V).
# Verified: same sonic locus as KHA (D=3N^2V^2-N^2+4NV-V^2+3=0), correct center regularity
# dV/dx->0 as N->inf. Build fast RHS, seed the regular center, shoot the free parameter so the
# solution passes analytically through the sonic point, then continue to the dispersive state.
import sympy as sp, pickle, numpy as np
from scipy.optimize import brentq, fsolve

# --- build RHS ---
en_s,mo_s=pickle.load(open("fluid_manual.pkl","rb"))
en=sp.sympify(en_s); mo=sp.sympify(mo_s)
A0,N0,o0,V0,Ax,Nx,ox,Vx=sp.symbols('A0 N0 om0 V0 A_x N_x om_x V_x')
Axv=A0*(1-A0+(2*o0/(1-V0**2))*(1+V0**2/3)); Nxv=N0*(-2+A0-sp.Rational(2,3)*o0)
sol=sp.solve([en.subs({Ax:Axv,Nx:Nxv}),mo.subs({Ax:Axv,Nx:Nxv})],[ox,Vx],dict=True)[0]
omx=sp.simplify(sol[ox]); Vxx=sp.simplify(sol[Vx])
f_om=sp.lambdify((A0,N0,o0,V0),omx,'numpy',cse=True)
f_V =sp.lambdify((A0,N0,o0,V0),Vxx,'numpy',cse=True)
def Axf(A,N,om,V): return A*(1-A+(2*om/(1-V*V))*(1+V*V/3.0))
def Nxf(A,N,om,V): return N*(-2+A-2.0*om/3.0)
def slopes(Y):
    N,A,om,V=Y
    return np.array([Nxf(A,N,om,V),Axf(A,N,om,V),float(f_om(A,N,om,V)),float(f_V(A,N,om,V))])
def Dloc(N,V): return 3*N*N*V*V-N*N+4*N*V-V*V+3

# --- center seed ---
# Near center z=e^x->0: leading auto-satisfied (N->inf). Determine (a2,w,v1) consistent seed by
# requiring the ODE slopes match the ansatz derivatives at a small z0, with omega>0 (w>0).
# free params: v1 (velocity coeff) and overall density -> we FIX nc by n1? Use 3 unknowns
# (nc,v1,a2) with w=1.5*a2 (from om-resid), seed at z0, solve consistency of dV,dA,dN (3 eqs).
def seed_state(param, z0):
    # param = v1 (center velocity slope). Solve for (nc,a2) s.t. ODE consistent; w=1.5 a2.
    v1=param
    def res(p):
        nc,a2=p
        if a2<=0: a2=1e-9
        w=1.5*a2
        N=nc/z0; A=1+a2*z0*z0; om=w*z0*z0; V=v1*z0
        dN,dA,dom,dV=slopes([N,A,om,V])
        return [dA-2*a2*z0*z0, dN-(-nc/z0)]   # match A' and N' leading; V' auto, om' auto
    for g in [[1.3,0.1],[1.7,0.05],[2.0,0.2],[1.2,0.02],[1.5,0.3]]:
        s=fsolve(res,g,full_output=True)
        if s[2]==1 and max(abs(np.array(res(s[0]))))<1e-9 and s[0][0]>0 and s[0][1]>0:
            nc,a2=s[0]; w=1.5*a2
            N=nc/z0;A=1+a2*z0*z0;om=w*z0*z0;V=v1*z0
            return np.array([N,A,om,V])
    return None

def integrate_from_center(v1, z0=np.exp(-12), xstop=8.0, h=5e-4):
    Y=seed_state(v1,z0)
    if Y is None: return None
    x=np.log(z0); traj=[(x,*Y)]; sonic=None
    n=int((xstop-x)/h)
    for i in range(n):
        Yp=Y.copy()
        k1=slopes(Y);k2=slopes(Y+0.5*h*k1);k3=slopes(Y+0.5*h*k2);k4=slopes(Y+h*k3)
        Y=Y+(h/6)*(k1+2*k2+2*k3+k4);x+=h
        if not np.all(np.isfinite(Y)) or abs(Y[3])>=0.99999 or Y[1]<=0: break
        D0=Dloc(Yp[0],Yp[3]); D1=Dloc(Y[0],Y[3])
        if sonic is None and D0*D1<0:
            frac=D0/(D0-D1); xs=x-h+frac*h
            Ys=Yp+frac*(Y-Yp); sonic=(xs,*Ys)
        traj.append((x,*Y))
    return traj, sonic

# The shoot: as we cross the sonic point, for a GENERIC v1 the solution is singular there
# (dV/dx blows up because det->0 but numerator!=0). The critical v1 makes numerator->0 too
# (analytic passage). Measure the regularity residual = numerator of V' at the sonic crossing.
def sonic_residual(v1):
    r=integrate_from_center(v1)
    if r is None or r[1] is None: return np.nan
    xs,N,A,om,V=r[1]
    # residual: the V' numerator at sonic (should ->0 for analytic). Evaluate f_V*det-> use the
    # fact that near sonic V' ~ numV/det; a regular soln has finite V'. Proxy: |V'| stays bounded.
    # Better: evaluate the actual numerator numV at the sonic state.
    return N,A,om,V,xs

print("scan center velocity slope v1 -> sonic-point state (looking for analytic passage):")
for v1 in np.linspace(0.05,1.2,24):
    r=integrate_from_center(v1)
    if r is None:
        print(f"  v1={v1:.3f}: seed fail"); continue
    traj,sonic=r
    xend=traj[-1][0]
    if sonic:
        xs,N,A,om,V=sonic
        print(f"  v1={v1:.3f}: SONIC x={xs:.3f} N={N:.4f} A={A:.4f} om={om:.5f} V={V:.4f}  (int to x={xend:.2f}, {len(traj)} steps)")
    else:
        print(f"  v1={v1:.3f}: no sonic, ended x={xend:.2f} A={traj[-1][2]:.3f} V={traj[-1][4]:.3f}")
