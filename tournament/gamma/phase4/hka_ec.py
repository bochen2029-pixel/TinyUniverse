# hka_ec.py — Evans-Coleman background: robust center->sonic shoot with scipy solve_ivp + event.
#
# STRUCTURE (established by hka_grazing.py / hka_refine_oi.py):
#  - Launch from the regular center (4.12-4.13) with gauge N_inf and central-density param oi.
#  - The 1-parameter family (in oi) all reach the sonic locus Dson=0 TRANSVERSALLY, with A->3/2
#    (2m/r=1/3) always. The EC solution is the member whose sonic crossing is ANALYTIC:
#    the L'Hopital / row-proportionality residual (4.6) vanishes AT the crossing.
#  - Use solve_ivp (RK45, tight tol) with a terminal event Dson=0 to get the crossing state to
#    ~machine precision (no linear-interpolation noise), then brentq oi on lhop(crossing)=0.
#  - The converged sonic point comes out at (A0,N0,om0,V0) = (3/2, 2/sqrt3, 3/4, -1/sqrt3) EXACTLY
#    (V0=-c_s), matching HKA (4.7-4.9). This IS the Evans-Coleman critical CSS solution.
#
# Eq #s per HKA_beta_equations.md.  Deterministic; no target tuned.

import numpy as np, math
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

G=4.0/3.0; GM1=G-1.0; CS=math.sqrt(GM1)
M_CENTER=-2.0/(3.0*G)

def A_of(N,om,V):
    oV2=1-V*V
    return 1 + 2*om*(1+GM1*V*V)/oV2 + 2*G*N*V*om/oV2

def Dson(N,V):
    return (1+N*V)**2 - GM1*(N+V)**2

def _fluid_slopes(A,N,om,V):
    g=G; oV2=1-V*V
    RHS_d=(3*(2-g)/2)*N*V-((2+g)/2)*A*N*V+(2-g)*N*V*om
    RHS_e=(2-g)*(g-1)*N*om+((7*g-6)/2)*N+((2-3*g)/2)*A*N
    a=(1+N*V)/om; b=g*(N+V)/oV2; c=(g-1)*(N+V)/om; d=g*(1+N*V)/oV2
    det=a*d-b*c
    return (d*RHS_d-b*RHS_e)/det, (a*RHS_e-c*RHS_d)/det

def rhs3(x, Y3):
    N,om,V=Y3
    A=A_of(N,om,V)
    Nx=N*(-2+A-(2.0-G)*om)
    omx,Vx=_fluid_slopes(A,N,om,V)
    return [Nx,omx,Vx]

def lhop(A,N,om,V):
    """Row-proportionality residual (4.6): (1+NV) RHS_e - (g-1)(N+V) RHS_d. =0 at analytic sonic passage."""
    g=G
    RHS_d=(3*(2-g)/2)*N*V-((2+g)/2)*A*N*V+(2-g)*N*V*om
    RHS_e=(2-g)*(g-1)*N*om+((7*g-6)/2)*N+((2-3*g)/2)*A*N
    return (1+N*V)*RHS_e-(g-1)*(N+V)*RHS_d

def center_state3(N_inf, oi, z):
    Vi=M_CENTER/N_inf; N=N_inf/z; om=oi*z*z; V=Vi*z
    return np.array([N,om,V])

def _event_Dson(x, Y3):
    return Dson(Y3[0], Y3[2])
_event_Dson.terminal=True
_event_Dson.direction=+1     # Dson rising through 0 (from center-negative to positive)

def shoot_to_sonic(N_inf, oi, x0=-12.0, rtol=1e-11, atol=1e-12, xmax=3.0):
    """Integrate center->outward; stop at first Dson=0 crossing. Return crossing state + lhop."""
    z0=math.exp(x0); Y0=center_state3(N_inf, oi, z0)
    sol=solve_ivp(rhs3, [x0, xmax], Y0, method='RK45', rtol=rtol, atol=atol,
                  events=_event_Dson, dense_output=False, max_step=0.05)
    if sol.status==1 and len(sol.t_events[0])>0:
        xe=sol.t_events[0][0]; Ye=sol.y_events[0][0]
        N,om,V=Ye; A=A_of(N,om,V)
        return dict(status='sonic', x=xe, N=N, om=om, V=V, A=A, lhop=lhop(A,N,om,V))
    return dict(status=sol.status, x=sol.t[-1], Y=sol.y[:,-1])

def _F(oi, N_inf=1.0, **kw):
    r=shoot_to_sonic(N_inf, oi, **kw)
    return r['lhop'] if r['status']=='sonic' else None

def find_ec(N_inf=1.0, bracket=(0.30,0.45), **kw):
    """brentq oi on lhop(sonic)=0. Returns (oi*, crossing dict)."""
    lo,hi=bracket
    flo=_F(lo,N_inf,**kw); fhi=_F(hi,N_inf,**kw)
    # widen/scan if not bracketed
    if flo is None or fhi is None or flo*fhi>0:
        ois=np.linspace(lo,hi,31); vals=[]
        for o in ois:
            v=_F(o,N_inf,**kw); vals.append(v)
        # find a sign-change pair
        br=None
        for i in range(1,len(ois)):
            if vals[i-1] is not None and vals[i] is not None and vals[i-1]*vals[i]<0:
                br=(ois[i-1],ois[i]); break
        if br is None:
            raise RuntimeError(f"no lhop sign change in {bracket}; vals={[None if v is None else round(v,3e-1) for v in vals]}")
        lo,hi=br
    oistar=brentq(lambda o:_F(o,N_inf,**kw), lo, hi, xtol=1e-12, rtol=1e-13)
    return oistar, shoot_to_sonic(N_inf, oistar, **kw)

if __name__=="__main__":
    import sys, time
    N_inf=1.0
    if len(sys.argv)>1: N_inf=float(sys.argv[1])
    t0=time.time()
    print("h-independent shoot via solve_ivp(RK45, dense event). Converge oi vs tolerance:")
    for rtol in (1e-9,1e-10,1e-11,1e-12):
        try:
            oistar, r = find_ec(N_inf, rtol=rtol, atol=rtol*0.1)
        except Exception as e:
            print(f" rtol={rtol:.0e}: {e}"); continue
        print(f" rtol={rtol:.0e}: oi*={oistar:.9f}  sonic x={r['x']:.5f} "
              f"V={r['V']:.8f} N={r['N']:.8f} om={r['om']:.8f} A={r['A']:.8f} lhop={r['lhop']:.2e}")
    print(f"\n reference exact sonic (V0=-1/sqrt3): V0={-1/math.sqrt(3):.8f} N0={2/math.sqrt(3):.8f} om0=0.75 A0=1.5")
    print(f"# {time.time()-t0:.1f}s")
