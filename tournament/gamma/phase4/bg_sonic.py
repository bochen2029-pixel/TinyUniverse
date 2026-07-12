# Background via sonic-point saddle separatrices, using derived (correct) Euler ODEs.
# At the sonic point det=0: it's a saddle. The two analytic solutions leave along the two
# eigendirections of the local flow. Integrate each separatrix; one -> regular center, the
# other -> dispersal. The sonic point is parametrized by V0 (with N0,A0,om0 from regularity).
# The physical (critical) solution is the V0 whose center-separatrix hits A->1,V->0.
import sympy as sp, pickle, numpy as np
from scipy.optimize import fsolve

en_s,mo_s=pickle.load(open("fluid_manual.pkl","rb"))
en=sp.sympify(en_s); mo=sp.sympify(mo_s)
A0s,N0s,o0s,V0s,Ax,Nx,ox,Vx=sp.symbols('A0 N0 om0 V0 A_x N_x om_x V_x')
Axv=A0s*(1-A0s+(2*o0s/(1-V0s**2))*(1+V0s**2/3)); Nxv=N0s*(-2+A0s-sp.Rational(2,3)*o0s)
e1=en.subs({Ax:Axv,Nx:Nxv}); e2=mo.subs({Ax:Axv,Nx:Nxv})
sol=sp.solve([e1,e2],[ox,Vx],dict=True)[0]
omx=sp.simplify(sol[ox]); Vxx=sp.simplify(sol[Vx])
f_om=sp.lambdify((A0s,N0s,o0s,V0s),omx,'numpy',cse=True)
f_V =sp.lambdify((A0s,N0s,o0s,V0s),Vxx,'numpy',cse=True)
# numerators for regularity at sonic (both must vanish for analytic passage):
# my 2x2: [[a11,a12],[a21,a22]] (ox,Vx)=(b1,b2). numV = det[[a11,b1],[a21,b2]].
e1e=sp.expand(e1); e2e=sp.expand(e2)
a11=e1e.coeff(ox);a12=e1e.coeff(Vx);a21=e2e.coeff(ox);a22=e2e.coeff(Vx)
b1=-(e1e.subs({ox:0,Vx:0}));b2=-(e2e.subs({ox:0,Vx:0}))
numV=sp.simplify(a11*b2-b1*a21); numU=sp.simplify(b1*a22-a12*b2)
f_numV=sp.lambdify((A0s,N0s,o0s,V0s),numV,'numpy',cse=True)
f_numU=sp.lambdify((A0s,N0s,o0s,V0s),numU,'numpy',cse=True)

def Axf(A,N,om,V): return A*(1-A+(2*om/(1-V*V))*(1+V*V/3.0))
def Nxf(A,N,om,V): return N*(-2+A-2.0*om/3.0)
def slopes(Y):
    N,A,om,V=Y; return np.array([Nxf(A,N,om,V),Axf(A,N,om,V),float(f_om(A,N,om,V)),float(f_V(A,N,om,V))])
def branch1_N(V): return (-2*V+np.sqrt(3)*(V*V-1))/(3*V*V-1)
def Dloc(N,V): return 3*N*N*V*V-N*N+4*N*V-V*V+3

def sonic_state(V0):
    # regularity at sonic: numV=numU=0 (both), solve (A0,om0) given V0,N0=branch1.
    N0=branch1_N(V0)
    def res(p):
        A,om=p; return [float(f_numV(A,N0,om,V0)), float(f_numU(A,N0,om,V0))]
    for g in [[1.7,0.37],[1.5,0.3],[1.1,0.01],[2,0.5],[1.2,0.05],[1.0,0.2],[1.8,0.4]]:
        s=fsolve(res,g,full_output=True)
        if s[2]==1 and max(abs(np.array(res(s[0]))))<1e-9 and s[0][0]>0 and s[0][1]>0:
            return N0,s[0][0],s[0][1]
    return N0,None,None

def sonic_eigendirs(V0,N0,A0,om0):
    # analytic slope dY/dx at sonic = L'Hopital. Since ox,Vx=0/0, use: near sonic the flow is
    # dY/dx with om',V' the finite limit. Compute by evaluating slopes at (state + eps along a
    # small displacement that stays ~on the analytic branch). Use Jacobian of (Nx,Ax,numV/det..).
    # Simpler robust: get the two slopes from L'Hopital quadratic. Numerically differentiate
    # numV,numU,det along candidate dY and solve. We'll estimate via small finite step both ways.
    eps=1e-6
    def det(N,A,om,V):
        return a11f(N,A,om,V)*a22f(N,A,om,V)-a12f(N,A,om,V)*a21f(N,A,om,V)
    return None

# Since both numerators vanish at sonic (regularity), dV/dx there is 0/0. Step off with the
# metric slopes exact and om',V' from a one-sided L'Hopital: evaluate slopes at V0+-tiny using
# the ACTUAL f_V (which is numV/det -> finite limit by construction since numV,det both ~ linear).
def step_off(V0,N0,A0,om0,sgn,h0=1e-4):
    # move a hair along V then let f_V,f_om give the (finite) slopes
    Vt=V0+sgn*1e-7
    # near sonic f_V should be finite (0/0 resolved). evaluate:
    Y=np.array([N0,A0,om0,V0])
    s=slopes(Y+np.array([0,0,0,sgn*1e-7]))
    Y2=Y+sgn*h0*s
    return Y2, s

def integrate(Y0, s0, direction, xspan, h=5e-4):
    Y=Y0.copy(); x=0.0; traj=[]
    n=int(abs(xspan/h)); d=np.sign(direction)
    for i in range(n):
        k1=slopes(Y);k2=slopes(Y+0.5*h*d*k1);k3=slopes(Y+0.5*h*d*k2);k4=slopes(Y+h*d*k3)
        Y=Y+(h*d/6)*(k1+2*k2+2*k3+k4);x+=h*d
        if not np.all(np.isfinite(Y)) or abs(Y[3])>=1 or Y[1]<=0 or Y[1]>100: break
        traj.append((x,*Y))
    return traj

print("sonic-state scan (my eqs), then integrate both directions:")
print("target center: A->1, V->0")
for V0 in np.linspace(0.05,0.55,26):
    N0,A0,om0=sonic_state(V0)
    if A0 is None: continue
    # step off both signs, integrate toward center (expect decreasing x) and dispersal
    for sgn in [+1,-1]:
        Y2,s0=step_off(V0,N0,A0,om0,sgn)
        # integrate in -x (toward center) and +x
        tin=integrate(Y2,s0,-1,8.0); tout=integrate(Y2,s0,+1,8.0)
        ei=tin[-1] if tin else None; eo=tout[-1] if tout else None
        def fmt(e): return f"x={e[0]:.1f},A={e[2]:.3f},V={e[4]:.3f},om={e[3]:.4f}" if e else "blow"
        # only print if one side approaches center
        if (ei and abs(ei[4])<0.05 and abs(ei[2]-1)<0.1) or (eo and abs(eo[4])<0.05 and abs(eo[2]-1)<0.1):
            print(f"  V0={V0:.3f} sgn={sgn:+d} N0={N0:.3f} A0={A0:.3f} om0={om0:.4f}: in[{fmt(ei)}] out[{fmt(eo)}]")
