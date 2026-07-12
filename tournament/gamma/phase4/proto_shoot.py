# Full background shoot: from the sonic point x=0, integrate inward (x->-inf, the center)
# and require A->1, V->0. The shoot parameter is V0 (with A0,omega0 fixed by regularity).
# The slope at the sonic point is 0/0; get it by L'Hopital (analytic node).
import numpy as np
from scipy.optimize import fsolve, brentq

def branch1_N(V):  return (-2*V + np.sqrt(3)*(V*V-1))/(3*V*V-1)

def metric_slopes(N,A,om,V):
    Ap = A*(1 - A + (2*om/(1-V*V))*(1 + V*V/3.0))
    Np = N*(-2 + A - 2.0*om/3.0)
    return Ap, Np

def slopes(N,A,om,V):
    # returns (N',A',omega',V'); away from sonic uses 2x2 solve
    Ap, Np = metric_slopes(N,A,om,V)
    omV2 = 1 - V*V
    m11=(1+N*V); m12=4*(N+V)/(3*omV2)
    m21=(4*V+N+3*N*V*V); m22=4*(1+V*V+2*N*V)/omV2
    rhs1=-(-N*V*Ap/(3*A)+4*V*Np/3.0+2*N*(1+4*om/(9*omV2)))
    rhs2=-( N*omV2*Ap/A+4*(1+V*V)*Np+2*N*(1+3*V*V))
    det=m11*m22-m12*m21
    u=(rhs1*m22-m12*rhs2)/det
    Vp=(m11*rhs2-rhs1*m21)/det
    return Np, Ap, u*om, Vp, det

def sonic_state(V0):
    # solve (A0,omega0) so both sonic numerators vanish (regular node). branch1 N.
    N0 = branch1_N(V0)
    def res(p):
        A0,om0=p
        Ap,Np=metric_slopes(N0,A0,om0,V0)
        omV2=1-V0*V0
        m11=(1+N0*V0); m12=4*(N0+V0)/(3*omV2)
        m21=(4*V0+N0+3*N0*V0*V0); m22=4*(1+V0*V0+2*N0*V0)/omV2
        rhs1=-(-N0*V0*Ap/(3*A0)+4*V0*Np/3.0+2*N0*(1+4*om0/(9*omV2)))
        rhs2=-( N0*omV2*Ap/A0+4*(1+V0*V0)*Np+2*N0*(1+3*V0*V0))
        return [rhs1*m22-m12*rhs2, m11*rhs2-rhs1*m21]
    # prefer physical (A0>0,om0>0). try several guesses.
    best=None
    for g in [[1,1],[1.5,0.5],[2,1],[1,0.2],[3,2],[1.8,0.2]]:
        try:
            s=fsolve(res,g,full_output=True)
            sol,ier=s[0],s[2]
            if ier==1 and abs(res(sol)[0])+abs(res(sol)[1])<1e-9 and sol[0]>0 and sol[1]>0:
                if best is None or sol[1]<best[1]:  # smallest omega branch (EC-like, low density node)
                    best=sol
        except Exception: pass
    return N0, best

def sonic_slope(V0, N0, A0, om0):
    # L'Hopital: at sonic det=0, slope s=V' solves a quadratic (two roots -> nodal picks one).
    # Use: differentiate det along solution: d det/dx = grad(det).(N',A',om',s). Numerator num_V
    # also ->0; V' = lim num_V/det = (d num_V/dx)/(d det/dx). Set s=V', w=om'. Two eqs (from
    # num_U and num_V L'Hopital). Solve the 2x2-quadratic numerically for (s,w).
    eps=1e-6
    def detf(N,A,om,V):
        omV2=1-V*V
        m11=(1+N*V); m12=4*(N+V)/(3*omV2); m21=(4*V+N+3*N*V*V); m22=4*(1+V*V+2*N*V)/omV2
        return m11*m22-m12*m21
    def numU(N,A,om,V):
        Ap,Np=metric_slopes(N,A,om,V); omV2=1-V*V
        m12=4*(N+V)/(3*omV2); m22=4*(1+V*V+2*N*V)/omV2
        rhs1=-(-N*V*Ap/(3*A)+4*V*Np/3.0+2*N*(1+4*om/(9*omV2)))
        rhs2=-( N*omV2*Ap/A+4*(1+V*V)*Np+2*N*(1+3*V*V))
        return rhs1*m22-m12*rhs2
    def numV(N,A,om,V):
        Ap,Np=metric_slopes(N,A,om,V); omV2=1-V*V
        m11=(1+N*V); m21=(4*V+N+3*N*V*V)
        rhs1=-(-N*V*Ap/(3*A)+4*V*Np/3.0+2*N*(1+4*om/(9*omV2)))
        rhs2=-( N*omV2*Ap/A+4*(1+V*V)*Np+2*N*(1+3*V*V))
        return m11*rhs2-rhs1*m21
    Ap,Np=metric_slopes(N0,A0,om0,V0)
    # partial derivatives via finite difference
    def grad(f):
        dN=(f(N0+eps,A0,om0,V0)-f(N0-eps,A0,om0,V0))/(2*eps)
        dA=(f(N0,A0+eps,om0,V0)-f(N0,A0-eps,om0,V0))/(2*eps)
        do=(f(N0,A0,om0+eps,V0)-f(N0,A0,om0-eps,V0))/(2*eps)
        dV=(f(N0,A0,om0,V0+eps)-f(N0,A0,om0,V0-eps))/(2*eps)
        return dN,dA,do,dV
    gd=grad(detf); gu=grad(numU); gv=grad(numV)
    # unknowns s=V', w=om'. ddx f = gN*Np+gA*Ap+go*w+gV*s
    # V': s*ddx(det)=ddx(numV) -> s*(gd2)=gv2 where gd2=gd[0]*Np+gd[1]*Ap+gd[2]*w+gd[3]*s
    # w:  (w/om0? no) use u=om'/om so w=om'. from numU: u*ddx(det)=ddx(numU); u=w/om0
    from scipy.optimize import fsolve as fs
    def eqs(p):
        s,w=p
        ddx_det=gd[0]*Np+gd[1]*Ap+gd[2]*w+gd[3]*s
        ddx_u=gu[0]*Np+gu[1]*Ap+gu[2]*w+gu[3]*s
        ddx_v=gv[0]*Np+gv[1]*Ap+gv[2]*w+gv[3]*s
        return [s*ddx_det-ddx_v, (w/om0)*ddx_det-ddx_u]
    roots=[]
    for g in [[-1,-1],[1,1],[-0.5,0.5],[0.1,-2],[-2,2],[0.5,-0.5]]:
        try:
            r=fs(eqs,g,full_output=True)
            if r[2]==1 and abs(eqs(r[0])[0])+abs(eqs(r[0])[1])<1e-6:
                roots.append(tuple(r[0]))
        except Exception: pass
    # dedup
    uniq=[]
    for r in roots:
        if not any(abs(r[0]-q[0])<1e-4 and abs(r[1]-q[1])<1e-4 for q in uniq): uniq.append(r)
    return uniq

# integrate from sonic point inward toward center with RK4, small step, avoid x=0 exactly
def integrate_in(V0, slope_choice, x_end=-8.0, h=1e-4):
    N0,best=sonic_state(V0)
    if best is None: return None
    A0,om0=best
    slopes_sonic=sonic_slope(V0,N0,A0,om0)
    if not slopes_sonic: return None
    # pick slope: inward branch. choose the requested index if available
    sc=slope_choice % len(slopes_sonic)
    s_V,s_om=slopes_sonic[sc]
    # step off the sonic point analytically: x=-h0
    h0=1e-3
    N=N0+ (N0*(-2+A0-2*om0/3.0))*(-h0)
    A=A0+ (A0*(1-A0+(2*om0/(1-V0*V0))*(1+V0*V0/3)))*(-h0)
    om=om0+ s_om*(-h0)
    V=V0+ s_V*(-h0)
    x=-h0
    steps=int((x-x_end)/h)
    for i in range(steps):
        def f(N,A,om,V):
            Np,Ap,omp,Vp,det=slopes(N,A,om,V); return np.array([Np,Ap,omp,Vp])
        y=np.array([N,A,om,V])
        k1=f(*y); k2=f(*(y-0.5*h*k1)); k3=f(*(y-0.5*h*k2)); k4=f(*(y-h*k3))
        y=y-(h/6)*(k1+2*k2+2*k3+k4)
        N,A,om,V=y; x-=h
        if not np.all(np.isfinite(y)) or abs(V)>0.999: break
    return x,N,A,om,V

print("shoot inward from sonic point; want A->1, V->0 at center (x->-inf):")
print("V0      slopes_at_sonic        (x_end,N,A,om,V)")
for V0 in [-0.4,-0.35,-0.3,-0.25,-0.2,0.2,0.25,0.3,0.35,0.4]:
    N0,best=sonic_state(V0)
    if best is None:
        print(f"{V0:+.2f}  no regular sonic state"); continue
    A0,om0=best
    sl=sonic_slope(V0,N0,A0,om0)
    print(f"{V0:+.2f}  N0={N0:.3f} A0={A0:.3f} om0={om0:.3f} slopes(V',om')={[(round(a,3),round(b,3)) for a,b in sl]}")
    for sc in range(len(sl)):
        res=integrate_in(V0,sc)
        if res:
            x,N,A,om,V=res
            print(f"        slope#{sc}: x_end={x:.2f} N={N:.3f} A={A:.4f} om={om:.4f} V={V:.4f}")
