# Robust two-sided shoot from the sonic point x=0 with a proper local Taylor step.
# Sonic state (N0,A0,om0,V0) fixed by V0 (both Cramer numerators vanish -> regular node).
# Analytic slope at sonic: L'Hopital gives a QUADRATIC for (V',om'); the "nodal"/analytic
# direction is the one that lets the solution reach BOTH the center and infinity. We pick a
# slope, integrate inward (x<0) and outward (x>0), and read the asymptotics.
#
# Anchors (Evans-Coleman gr-qc/9402041): as x->+inf, A=a^2->~1.145, omega->~0.0109, V->? (outgoing,
# V>0). As x->-inf (center): A->1, V->0. Similarity exponent n~1.1485 (cross-check via N?).
import numpy as np
from scipy.optimize import fsolve

def branch1_N(V): return (-2*V+np.sqrt(3)*(V*V-1))/(3*V*V-1)
def branch0_N(V): return (-2*V-np.sqrt(3)*(V*V-1))/(3*V*V-1)

def metric_slopes(N,A,om,V):
    Ap=A*(1-A+(2*om/(1-V*V))*(1+V*V/3.0)); Np=N*(-2+A-2.0*om/3.0); return Ap,Np

def full_slopes(N,A,om,V):
    Ap,Np=metric_slopes(N,A,om,V); omV2=1-V*V
    m11=(1+N*V);m12=4*(N+V)/(3*omV2);m21=(4*V+N+3*N*V*V);m22=4*(1+V*V+2*N*V)/omV2
    rhs1=-(-N*V*Ap/(3*A)+4*V*Np/3.0+2*N*(1+4*om/(9*omV2)))
    rhs2=-( N*omV2*Ap/A+4*(1+V*V)*Np+2*N*(1+3*V*V))
    det=m11*m22-m12*m21
    u=(rhs1*m22-m12*rhs2)/det; Vp=(m11*rhs2-rhs1*m21)/det
    return Np,Ap,u*om,Vp,det

def sonic_state(V0, Nbranch=branch1_N):
    N0=Nbranch(V0)
    def res(p):
        A0,om0=p; Ap,Np=metric_slopes(N0,A0,om0,V0); omV2=1-V0*V0
        m11=(1+N0*V0);m12=4*(N0+V0)/(3*omV2);m21=(4*V0+N0+3*N0*V0*V0);m22=4*(1+V0*V0+2*N0*V0)/omV2
        rhs1=-(-N0*V0*Ap/(3*A0)+4*V0*Np/3.0+2*N0*(1+4*om0/(9*omV2)))
        rhs2=-( N0*omV2*Ap/A0+4*(1+V0*V0)*Np+2*N0*(1+3*V0*V0))
        return [rhs1*m22-m12*rhs2, m11*rhs2-rhs1*m21]
    for g in [[1.1,0.01],[1.15,0.011],[1.0,0.02],[1.2,0.05],[1.5,0.5],[2,1]]:
        s=fsolve(res,g,full_output=True)
        if s[2]==1 and abs(res(s[0])[0])+abs(res(s[0])[1])<1e-10 and s[0][0]>0 and s[0][1]>0:
            return N0,s[0][0],s[0][1]
    return N0,None,None

def sonic_slopes(V0,N0,A0,om0):
    # L'Hopital via FD partials. Solve (s=V',w=om') from:
    #   s*ddx(det)=ddx(numV);  (w/om0)*ddx(det)=ddx(numU)  where ddx uses N',A' exact + (w,s).
    eps=1e-7
    def detf(N,A,om,V):
        omV2=1-V*V
        return (1+N*V)*(4*(1+V*V+2*N*V)/omV2)-(4*(N+V)/(3*omV2))*(4*V+N+3*N*V*V)
    def numUf(N,A,om,V):
        Ap,Np=metric_slopes(N,A,om,V);omV2=1-V*V
        m12=4*(N+V)/(3*omV2);m22=4*(1+V*V+2*N*V)/omV2
        rhs1=-(-N*V*Ap/(3*A)+4*V*Np/3.0+2*N*(1+4*om/(9*omV2)));rhs2=-(N*omV2*Ap/A+4*(1+V*V)*Np+2*N*(1+3*V*V))
        return rhs1*m22-m12*rhs2
    def numVf(N,A,om,V):
        Ap,Np=metric_slopes(N,A,om,V);omV2=1-V*V
        m11=(1+N*V);m21=(4*V+N+3*N*V*V)
        rhs1=-(-N*V*Ap/(3*A)+4*V*Np/3.0+2*N*(1+4*om/(9*omV2)));rhs2=-(N*omV2*Ap/A+4*(1+V*V)*Np+2*N*(1+3*V*V))
        return m11*rhs2-rhs1*m21
    Ap,Np=metric_slopes(N0,A0,om0,V0)
    def g(f):
        return ((f(N0+eps,A0,om0,V0)-f(N0-eps,A0,om0,V0))/(2*eps),
                (f(N0,A0+eps,om0,V0)-f(N0,A0-eps,om0,V0))/(2*eps),
                (f(N0,A0,om0+eps,V0)-f(N0,A0,om0-eps,V0))/(2*eps),
                (f(N0,A0,om0,V0+eps)-f(N0,A0,om0,V0-eps))/(2*eps))
    gd=g(detf);gu=g(numUf);gv=g(numVf)
    def eqs(p):
        s,w=p
        dd=gd[0]*Np+gd[1]*Ap+gd[2]*w+gd[3]*s
        du=gu[0]*Np+gu[1]*Ap+gu[2]*w+gu[3]*s
        dv=gv[0]*Np+gv[1]*Ap+gv[2]*w+gv[3]*s
        return [s*dd-dv,(w/om0)*dd-du]
    out=[]
    for gg in [[-0.5,-0.5],[0.5,0.5],[-1,1],[1,-1],[0,-2],[-0.3,0.3],[2,2],[-2,-2]]:
        r=fsolve(eqs,gg,full_output=True)
        if r[2]==1 and abs(eqs(r[0])[0])+abs(eqs(r[0])[1])<1e-6:
            if not any(abs(r[0][0]-o[0])<1e-4 and abs(r[0][1]-o[1])<1e-4 for o in out): out.append(tuple(r[0]))
    return out

def integrate(V0,slope,direction, x_far, h=2e-4):
    N0,A0,om0=sonic_state(V0)
    if A0 is None: return None
    Np0,Ap0=metric_slopes(N0,A0,om0,V0)  # exact metric slopes at sonic
    sV,som=slope
    d=np.sign(direction)
    h0=1e-3*d
    # 1st-order Taylor step off sonic
    N=N0+Np0*h0; A=A0+Ap0*h0; om=om0+som*h0; V=V0+sV*h0; x=h0
    n=int(abs((x_far-x)/h))
    traj=[]
    for i in range(n):
        y=np.array([N,A,om,V])
        def f(Y):
            Np,Ap,omp,Vp,det=full_slopes(*Y); return np.array([Np,Ap,omp,Vp])
        k1=f(y);k2=f(y+0.5*h*d*k1);k3=f(y+0.5*h*d*k2);k4=f(y+h*d*k3)
        y=y+(h*d/6)*(k1+2*k2+2*k3+k4)
        N,A,om,V=y;x+=h*d
        if not np.all(np.isfinite(y)) or abs(V)>=1.0 or A<=0 or A>50: break
        traj.append((x,N,A,om,V))
    return traj[-1] if traj else None

for V0 in np.linspace(-0.5,0.5,21):
    N0,A0,om0=sonic_state(V0)
    if A0 is None: continue
    sls=sonic_slopes(V0,N0,A0,om0)
    print(f"V0={V0:+.3f} N0={N0:.3f} A0={A0:.3f} om0={om0:.4f}  slopes={[(round(a,3),round(b,3)) for a,b in sls]}")
    for j,sl in enumerate(sls):
        ri=integrate(V0,sl,-1,-7.0)   # inward -> center
        ro=integrate(V0,sl,+1,+7.0)   # outward -> infinity
        si=f"in:x={ri[0]:.1f},A={ri[2]:.3f},V={ri[3]:.3f}" if ri else "in:blow"
        so=f"out:x={ro[0]:.1f},A={ro[2]:.3f},V={ro[3]:.3f},om={ro[3]:.3f}" if ro else "out:blow"
        print(f"    slope#{j} {si}  {so}")
