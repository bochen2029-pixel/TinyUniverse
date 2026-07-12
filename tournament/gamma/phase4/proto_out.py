# Validate the OUTWARD solution against Evans-Coleman infinity asymptotics:
#   A(x->+inf) -> a^2 ~ 1.145 (a~1.07),  omega -> ~0.0109,  V -> V_inf>0 (outgoing).
# Scan V0 (sonic velocity) and the sonic eigendirection; integrate outward far; record asymptote.
import numpy as np
from scipy.optimize import fsolve

def branch1_N(V): return (-2*V+np.sqrt(3)*(V*V-1))/(3*V*V-1)
def branch0_N(V): return (-2*V-np.sqrt(3)*(V*V-1))/(3*V*V-1)
def metric_slopes(N,A,om,V):
    Ap=A*(1-A+(2*om/(1-V*V))*(1+V*V/3.0)); Np=N*(-2+A-2.0*om/3.0); return Ap,Np
def full_slopes(Y):
    N,A,om,V=Y; Ap,Np=metric_slopes(N,A,om,V); omV2=1-V*V
    m11=(1+N*V);m12=4*(N+V)/(3*omV2);m21=(4*V+N+3*N*V*V);m22=4*(1+V*V+2*N*V)/omV2
    rhs1=-(-N*V*Ap/(3*A)+4*V*Np/3.0+2*N*(1+4*om/(9*omV2)))
    rhs2=-( N*omV2*Ap/A+4*(1+V*V)*Np+2*N*(1+3*V*V))
    det=m11*m22-m12*m21
    return np.array([Np,Ap,(rhs1*m22-m12*rhs2)/det*om,(m11*rhs2-rhs1*m21)/det])

def sonic_state(V0,Nb):
    N0=Nb(V0)
    def res(p):
        A0,om0=p; Ap,Np=metric_slopes(N0,A0,om0,V0); omV2=1-V0*V0
        m11=(1+N0*V0);m12=4*(N0+V0)/(3*omV2);m21=(4*V0+N0+3*N0*V0*V0);m22=4*(1+V0*V0+2*N0*V0)/omV2
        rhs1=-(-N0*V0*Ap/(3*A0)+4*V0*Np/3.0+2*N0*(1+4*om0/(9*omV2)))
        rhs2=-( N0*omV2*Ap/A0+4*(1+V0*V0)*Np+2*N0*(1+3*V0*V0))
        return [rhs1*m22-m12*rhs2, m11*rhs2-rhs1*m21]
    sols=[]
    for g in [[1.1,0.01],[1.15,0.011],[1.0,0.02],[1.2,0.05],[1.5,0.5],[2,1],[1.05,0.005]]:
        s=fsolve(res,g,full_output=True)
        if s[2]==1 and abs(res(s[0])[0])+abs(res(s[0])[1])<1e-11 and s[0][0]>0 and s[0][1]>0:
            if not any(abs(s[0][0]-q[0])<1e-5 for q in sols): sols.append(tuple(s[0]))
    return N0,sols

def sonic_slopes(V0,N0,A0,om0):
    eps=1e-7
    def detf(N,A,om,V):
        omV2=1-V*V; return (1+N*V)*(4*(1+V*V+2*N*V)/omV2)-(4*(N+V)/(3*omV2))*(4*V+N+3*N*V*V)
    def numUf(N,A,om,V):
        Ap,Np=metric_slopes(N,A,om,V);omV2=1-V*V;m12=4*(N+V)/(3*omV2);m22=4*(1+V*V+2*N*V)/omV2
        rhs1=-(-N*V*Ap/(3*A)+4*V*Np/3.0+2*N*(1+4*om/(9*omV2)));rhs2=-(N*omV2*Ap/A+4*(1+V*V)*Np+2*N*(1+3*V*V))
        return rhs1*m22-m12*rhs2
    def numVf(N,A,om,V):
        Ap,Np=metric_slopes(N,A,om,V);omV2=1-V*V;m11=(1+N*V);m21=(4*V+N+3*N*V*V)
        rhs1=-(-N*V*Ap/(3*A)+4*V*Np/3.0+2*N*(1+4*om/(9*omV2)));rhs2=-(N*omV2*Ap/A+4*(1+V*V)*Np+2*N*(1+3*V*V))
        return m11*rhs2-rhs1*m21
    Ap,Np=metric_slopes(N0,A0,om0,V0)
    def g(f): return ((f(N0+eps,A0,om0,V0)-f(N0-eps,A0,om0,V0))/(2*eps),(f(N0,A0+eps,om0,V0)-f(N0,A0-eps,om0,V0))/(2*eps),(f(N0,A0,om0+eps,V0)-f(N0,A0,om0-eps,V0))/(2*eps),(f(N0,A0,om0,V0+eps)-f(N0,A0,om0,V0-eps))/(2*eps))
    gd=g(detf);gu=g(numUf);gv=g(numVf)
    def eqs(p):
        s,w=p; dd=gd[0]*Np+gd[1]*Ap+gd[2]*w+gd[3]*s;du=gu[0]*Np+gu[1]*Ap+gu[2]*w+gu[3]*s;dv=gv[0]*Np+gv[1]*Ap+gv[2]*w+gv[3]*s
        return [s*dd-dv,(w/om0)*dd-du]
    out=[]
    for gg in [[-1,-1],[1,1],[-0.5,0.5],[0.5,-0.5],[0,-2],[2,2],[-2,2],[3,-1]]:
        r=fsolve(eqs,gg,full_output=True)
        if r[2]==1 and abs(eqs(r[0])[0])+abs(eqs(r[0])[1])<1e-6 and not any(abs(r[0][0]-o[0])<1e-4 for o in out): out.append(tuple(r[0]))
    return out

def integ_out(V0,A0,om0,N0,slope,xfar=12.0,h=1e-3):
    Np0,Ap0=metric_slopes(N0,A0,om0,V0); sV,som=slope; h0=1e-3
    Y=np.array([N0+Np0*h0,A0+Ap0*h0,om0+som*h0,V0+sV*h0]); x=h0
    n=int((xfar-x)/h); last=None
    for i in range(n):
        k1=full_slopes(Y);k2=full_slopes(Y+0.5*h*k1);k3=full_slopes(Y+0.5*h*k2);k4=full_slopes(Y+h*k3)
        Y=Y+(h/6)*(k1+2*k2+2*k3+k4);x+=h
        if not np.all(np.isfinite(Y)) or abs(Y[3])>=1 or Y[1]<=0 or Y[1]>10: return None
        last=(x,Y[0],Y[1],Y[2],Y[3])
    return last

print("Target (EC): A_inf~1.145, omega_inf~0.0109, V_inf>0 outgoing")
print("V0     branch  (N0,A0,om0)          slope       -> (x,N,A,om,V)_far")
for Nb,nm in [(branch1_N,'b1'),(branch0_N,'b0')]:
    for V0 in np.linspace(-0.55,0.55,23):
        N0,sols=sonic_state(V0,Nb)
        for (A0,om0) in sols:
            sls=sonic_slopes(V0,N0,A0,om0)
            for sl in sls:
                r=integ_out(V0,A0,om0,N0,sl)
                if r and 0.9<r[2]<1.4 and r[3]>1e-4:   # A in plausible range, omega>0 (not decayed)
                    print(f"{V0:+.3f} {nm}  ({N0:.3f},{A0:.3f},{om0:.4f})  ({sl[0]:+.2f},{sl[1]:+.2f}) -> x={r[0]:.1f} A={r[2]:.4f} om={r[3]:.5f} V={r[4]:.4f}")
