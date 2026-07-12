# Integrate the separatrices from the V0=0 sonic point (N0=sqrt3, A0=1.75, om0=0.375) and
# verify the center<-sonic->dispersal structure with the derived equations.
import sympy as sp, pickle, numpy as np
from scipy.optimize import fsolve
d=pickle.load(open("rhs_exprs.pkl","rb"))
A0,N0,o0,V0=sp.symbols('A0 N0 om0 V0')
omx=sp.sympify(d['omx']); Vxx=sp.sympify(d['Vxx']); det=sp.sympify(d['det'])
nV=sp.sympify(d['nV']); nU=sp.sympify(d['nU'])
f_om=sp.lambdify((A0,N0,o0,V0),omx,'numpy',cse=True)
f_V =sp.lambdify((A0,N0,o0,V0),Vxx,'numpy',cse=True)
f_det=sp.lambdify((A0,N0,o0,V0),det,'numpy',cse=True)
f_nV=sp.lambdify((A0,N0,o0,V0),nV,'numpy',cse=True)
f_nU=sp.lambdify((A0,N0,o0,V0),nU,'numpy',cse=True)
def Axf(A,N,om,V): return A*(1-A+(2*om/(1-V*V))*(1+V*V/3.0))
def Nxf(A,N,om,V): return N*(-2+A-2.0*om/3.0)
def slopes(Y):
    N,A,om,V=Y; return np.array([Nxf(A,N,om,V),Axf(A,N,om,V),float(f_om(A,N,om,V)),float(f_V(A,N,om,V))])

# sonic point state at V0=0
V0v=0.0; N0v=np.sqrt(3.0)
def res(p):
    A,om=p
    return [float(f_nV(A,N0v,om,V0v)), float(f_nU(A,N0v,om,V0v))]
sol=fsolve(res,[1.75,0.375],full_output=True)
A0v,om0v=sol[0]
print(f"sonic point: V0={V0v} N0={N0v:.6f} A0={A0v:.6f} om0={om0v:.6f} resid={max(abs(np.array(res(sol[0])))):.2e}")

# analytic slopes (finite) at sonic via L'Hopital
eps=1e-6
def g(f):
    return np.array([(f(A0v,N0v+eps,om0v,V0v)-f(A0v,N0v-eps,om0v,V0v))/(2*eps),
                     (f(A0v+eps,N0v,om0v,V0v)-f(A0v-eps,N0v,om0v,V0v))/(2*eps),
                     (f(A0v,N0v,om0v+eps,V0v)-f(A0v,N0v,om0v-eps,V0v))/(2*eps),
                     (f(A0v,N0v,om0v,V0v+eps)-f(A0v,N0v,om0v,V0v-eps))/(2*eps)])
gd=g(f_det); gV=g(f_nV); gU=g(f_nU)
Np=Nxf(A0v,N0v,om0v,V0v); Ap=Axf(A0v,N0v,om0v,V0v)
def eqs(p):
    s,w=p
    dd=gd[0]*Np+gd[1]*Ap+gd[2]*w+gd[3]*s
    dv=gV[0]*Np+gV[1]*Ap+gV[2]*w+gV[3]*s
    du=gU[0]*Np+gU[1]*Ap+gU[2]*w+gU[3]*s
    return [dv-s*dd, du-w*dd]
roots=[]
for gg in [[-1,-1],[1,1],[-.5,.5],[.5,-.5],[-2,2],[2,-2],[0,-2],[-.3,-.3],[3,3],[-3,3]]:
    rr=fsolve(eqs,gg,full_output=True)
    if rr[2]==1 and max(abs(np.array(eqs(rr[0]))))<1e-6 and not any(abs(rr[0][0]-o[0])<1e-4 and abs(rr[0][1]-o[1])<1e-4 for o in roots):
        roots.append(tuple(rr[0]))
print("analytic (dV/dx, dom/dx) eigen-slopes at sonic:", [(round(s,4),round(w,4)) for s,w in roots])
print(f"  (metric slopes there: N'={Np:.4f} A'={Ap:.4f})")

def integ(dY, dirn, xspan=8.0, h=1e-4):
    Y0=np.array([N0v,A0v,om0v,V0v]); Y=Y0+dirn*h*dY; x=dirn*h; traj=[(0,*Y0),(x,*Y)]
    for i in range(int(xspan/h)):
        k1=slopes(Y);k2=slopes(Y+0.5*h*dirn*k1);k3=slopes(Y+0.5*h*dirn*k2);k4=slopes(Y+h*dirn*k3)
        Y=Y+(h*dirn/6)*(k1+2*k2+2*k3+k4);x+=h*dirn
        if not np.all(np.isfinite(Y)) or abs(Y[3])>=1 or Y[1]<=0 or Y[1]>500: break
        traj.append((x,*Y))
    return traj
for j,(s,w) in enumerate(roots):
    dY=np.array([Np,Ap,w,s])
    for dr in [+1,-1]:
        tr=integ(dY,dr)
        e=tr[-1]
        tag=""
        if abs(e[4])<0.02 and abs(e[2]-1)<0.1: tag=" <== CENTER (A->1,V->0)"
        elif abs(e[4])>0.9: tag=" (dispersal V->1)"
        print(f"  eig{j} slope(V',om')=({s:+.3f},{w:+.3f}) dir{dr:+d}: end x={e[0]:+.2f} N={e[1]:.3f} A={e[2]:.4f} om={e[3]:.5f} V={e[4]:+.4f}{tag}")
