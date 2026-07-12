import sympy as sp, pickle, numpy as np
from scipy.optimize import fsolve, brentq
d=pickle.load(open("rhs_exprs.pkl","rb"))
A0,N0,o0,V0=sp.symbols('A0 N0 om0 V0')
omx=sp.sympify(d['omx']); Vxx=sp.sympify(d['Vxx']); nV=sp.sympify(d['nV']); nU=sp.sympify(d['nU']); det=sp.sympify(d['det'])
f_om=sp.lambdify((A0,N0,o0,V0),omx,'numpy',cse=True)
f_V =sp.lambdify((A0,N0,o0,V0),Vxx,'numpy',cse=True)
f_nV=sp.lambdify((A0,N0,o0,V0),nV,'numpy',cse=True)
f_nU=sp.lambdify((A0,N0,o0,V0),nU,'numpy',cse=True)
def Axf(A,N,om,V): return A*(1-A+(2*om/(1-V*V))*(1+V*V/3.0))
def Nxf(A,N,om,V): return N*(-2+A-2.0*om/3.0)
def slopes(Y):
    N,A,om,V=Y; return np.array([Nxf(A,N,om,V),Axf(A,N,om,V),float(f_om(A,N,om,V)),float(f_V(A,N,om,V))])
def branch1_N(V): return (-2*V+np.sqrt(3)*(V*V-1))/(3*V*V-1)
def Dloc(N,V): return 3*N*N*V*V-N*N+4*N*V-V*V+3

def sonic_state(V0v):
    N0v=branch1_N(V0v)
    def res(p):
        A,om=p
        if A<=0 or om<=0: return [1e3,1e3]
        return [float(f_nV(A,N0v,om,V0v)), float(f_nU(A,N0v,om,V0v))]
    for g in [[1.75,0.375],[1.5,0.3],[1.1,0.05],[2,0.5],[1.3,0.1],[1.9,0.4],[1.6,0.2]]:
        s=fsolve(res,g,full_output=True)
        if s[2]==1 and max(abs(np.array(res(s[0]))))<1e-8 and s[0][0]>0 and s[0][1]>0:
            return N0v,s[0][0],s[0][1]
    return N0v,None,None

# analytic slope at sonic (finite dY/dx). Near sonic det->0; slopes() would blow up. Get the
# finite limit by L'Hopital: dV/dx there solves quadratic. Estimate numerically: evaluate
# slopes() at two nearby points ON the analytic branch (approach sonic from both sides via a
# tiny parametric move that keeps nV,nU~0). Use central slope from f_V at V0+-delta after
# re-solving (A,om) regularity at those V -> traces the analytic separatrix.
def analytic_step(V0v, sgn, dv=1e-3):
    # move sonic-velocity by dv, get the regular state there (this IS the separatrix, since the
    # 1-param family of sonic-regular states parametrized by V0 traces the analytic solution
    # through the sonic point). This gives the neighboring point on the physical solution!
    Vn=V0v+sgn*dv
    Nn,An,omn=sonic_state(Vn)
    if An is None: return None
    return np.array([Nn,An,omn,Vn])

# THE INSIGHT: the set of sonic-regular states {(N0(V),A0(V),om0(V),V)} as V varies is itself a
# curve; the physical solution passes through ONE of them. But actually each V0 is a DIFFERENT
# candidate solution. So we can't trace along V0. We must integrate the ODE from a fixed sonic
# point. Do that: fix V0, step off along the analytic eigendirection (finite slope), integrate.
def sonic_slope_finite(V0v,N0v,A0v,om0v):
    # finite dY/dx via L'Hopital: differentiate det, nV, nU along solution and solve.
    eps=1e-6
    def g(f):
        return np.array([(f(A0v,N0v+eps,om0v,V0v)-f(A0v,N0v-eps,om0v,V0v))/(2*eps),
                         (f(A0v+eps,N0v,om0v,V0v)-f(A0v-eps,N0v,om0v,V0v))/(2*eps),
                         (f(A0v,N0v,om0v+eps,V0v)-f(A0v,N0v,om0v-eps,V0v))/(2*eps),
                         (f(A0v,N0v,om0v,V0v+eps)-f(A0v,N0v,om0v,V0v-eps))/(2*eps)])
    fdet=sp.lambdify((A0,N0,o0,V0),det,'numpy',cse=True)
    gd=g(fdet); gV=g(f_nV)
    Np=Nxf(A0v,N0v,om0v,V0v); Ap=Axf(A0v,N0v,om0v,V0v)
    # dV/dx=s, dom/dx=w. det=0, nV=0 at sonic. d(nV)/dx = s*d(det)/dx (from V'=nV/det, L'Hopital)
    # => (gV.[Np,Ap,w,s]) = s*(gd.[Np,Ap,w,s]). Also need om': from nU similarly.
    gU=g(f_nU)
    def eqs(p):
        s,w=p
        dd=gd[0]*Np+gd[1]*Ap+gd[2]*w+gd[3]*s
        dv=gV[0]*Np+gV[1]*Ap+gV[2]*w+gV[3]*s
        du=gU[0]*Np+gU[1]*Ap+gU[2]*w+gU[3]*s
        return [dv-s*dd, du-w*dd]
    outs=[]
    for gg in [[-1,-1],[1,1],[-0.5,0.5],[0.5,-0.5],[-2,2],[2,-2],[0,-2],[-0.3,-0.3]]:
        rr=fsolve(eqs,gg,full_output=True)
        if rr[2]==1 and max(abs(np.array(eqs(rr[0]))))<1e-6 and not any(abs(rr[0][0]-o[0])<1e-4 for o in outs):
            outs.append(tuple(rr[0]))
    return Np,Ap,outs

def integ(Y0,dY0,dirn,xspan,h=2e-4):
    Y=Y0+dirn*h*dY0; x=dirn*h; traj=[(x,*Y)]
    n=int(abs(xspan/h))
    for i in range(n):
        k1=slopes(Y);k2=slopes(Y+0.5*h*dirn*k1);k3=slopes(Y+0.5*h*dirn*k2);k4=slopes(Y+h*dirn*k3)
        Y=Y+(h*dirn/6)*(k1+2*k2+2*k3+k4);x+=h*dirn
        if not np.all(np.isfinite(Y)) or abs(Y[3])>=1 or Y[1]<=0 or Y[1]>200: break
        traj.append((x,*Y))
    return traj

print("For each sonic V0: eigendirections and where each separatrix goes:")
for V0v in np.linspace(0.1,0.5,17):
    N0v,A0v,om0v=sonic_state(V0v)
    if A0v is None: print(f" V0={V0v:.3f}: no sonic state"); continue
    Np,Ap,slps=sonic_slope_finite(V0v,N0v,A0v,om0v)
    dirs=[np.array([Np,Ap,w,s]) for (s,w) in slps]
    msg=f" V0={V0v:.3f} N0={N0v:.3f} A0={A0v:.3f} om0={om0v:.4f}: "
    for j,dY in enumerate(dirs):
        for dr in [+1,-1]:
            tr=integ(np.array([N0v,A0v,om0v,V0v]),dY,dr,7.0)
            e=tr[-1]
            if abs(e[4])<0.03 and abs(e[2]-1)<0.15:  # reached center-like
                msg+=f"[eig{j} dir{dr:+d}->CENTER x={e[0]:.1f} A={e[2]:.3f} V={e[4]:.3f}] "
    print(msg)
