# Trace the separatrix through the V0=0 sonic point robustly. Since V'=nV/det is 0/0 there,
# step off using a series: to first order use the metric slopes for N,A and a SMALL guess for
# the (V,om) step, then let the adaptive integrator recover the analytic branch a hair away.
# We approach from BOTH sides in x and glue. Use very fine RK4 and start the step-off at a point
# a small delta AWAY where det is small-but-nonzero, reading f_V,f_om (which are finite ratios).
import sympy as sp, pickle, numpy as np
d=pickle.load(open("rhs_exprs.pkl","rb"))
A0,N0,o0,V0=sp.symbols('A0 N0 om0 V0')
f_om=sp.lambdify((A0,N0,o0,V0),sp.sympify(d['omx']),'numpy',cse=True)
f_V =sp.lambdify((A0,N0,o0,V0),sp.sympify(d['Vxx']),'numpy',cse=True)
f_det=sp.lambdify((A0,N0,o0,V0),sp.sympify(d['det']),'numpy',cse=True)
def Axf(A,N,om,V): return A*(1-A+(2*om/(1-V*V))*(1+V*V/3.0))
def Nxf(A,N,om,V): return N*(-2+A-2.0*om/3.0)
def slopes(Y):
    N,A,om,V=Y
    return np.array([Nxf(A,N,om,V),Axf(A,N,om,V),float(f_om(A,N,om,V)),float(f_V(A,N,om,V))])

# sonic point
Ns=np.sqrt(3.0); As=1.75; oms=0.375; Vs=0.0
# The physical solution near the sonic point: expand to 2nd order. dY/dx finite. We estimate
# the analytic slope by evaluating f_V,f_om at points slightly displaced ALONG a trial dY and
# requiring self-consistency (fixed point of the slope map). Iterate:
Np=Nxf(As,Ns,oms,Vs); Ap=Axf(As,Ns,oms,Vs)
# trial: approach the sonic point along x by delta; the analytic soln has Y(x)=Ys+dY*x+...
# Evaluate slope at Ys+eps*dY_guess and iterate dY=(Np,Ap,f_om,f_V) there.
def refine_slope(sign, delta=1e-3, iters=60):
    dY=np.array([Np,Ap,0.0,0.0])  # start with metric slopes, zero V',om'
    for it in range(iters):
        Yt=np.array([Ns,As,oms,Vs])+sign*delta*dY
        st=slopes(Yt)
        dY=0.5*dY+0.5*np.array([Np,Ap,st[2],st[3]])   # damped fixed point
    return dY
for sign in [+1,-1]:
    dY=refine_slope(sign)
    print(f"sign={sign:+d}: refined slope dY=(N',A',om',V')=({dY[0]:.4f},{dY[1]:.4f},{dY[2]:.4f},{dY[3]:.4f})")
    # integrate away from sonic in +x and -x from the stepped-off point
    for dr in [+1,-1]:
        h=1e-4; delta=1e-3
        Y=np.array([Ns,As,oms,Vs])+sign*delta*dY  # step off
        x=sign*delta
        traj=[(x,*Y)]
        for i in range(int(8.0/h)):
            k1=slopes(Y);k2=slopes(Y+0.5*h*dr*k1);k3=slopes(Y+0.5*h*dr*k2);k4=slopes(Y+h*dr*k3)
            Y=Y+(h*dr/6)*(k1+2*k2+2*k3+k4);x+=h*dr
            if not np.all(np.isfinite(Y)) or abs(Y[3])>=1 or Y[1]<=0 or Y[1]>500: break
            traj.append((x,*Y))
        e=traj[-1]
        tag=""
        if abs(e[4])<0.02 and abs(e[2]-1)<0.15: tag=" <==CENTER"
        elif abs(e[4])>0.85: tag=" dispersal"
        print(f"    step sign{sign:+d} then integrate dir{dr:+d}: end x={e[0]:+.2f} N={e[1]:.3f} A={e[2]:.4f} om={e[3]:.5f} V={e[4]:+.4f}{tag}")
