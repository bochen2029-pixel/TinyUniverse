# Launch off the sonic point along the desingularized real eigenspace, then integrate the RESOLVED
# flow dY/dx=rhs(Y) with RK4 in both +/-x. Identify which side/eps-sign heads to the regular center
# (A->1, V->0, N->+inf as x->-inf).
import numpy as np, math
import hka_desing as D
import hka_background as HB

def resolved_rk4(Y, h):
    k1=D.resolved_rhs(Y); k2=D.resolved_rhs(Y+0.5*h*k1); k3=D.resolved_rhs(Y+0.5*h*k2); k4=D.resolved_rhs(Y+h*k3)
    return Y+(h/6)*(k1+2*k2+2*k3+k4)

def run(V0, eps, xdir, h=2e-4, nmax=200000, xspan=16.0):
    Y0,J,w,Vc = D.sonic_jacobian(V0)
    order=np.argsort(-w.real); vec=Vc[:,order[0]]
    vr=np.real(vec); vr/=np.linalg.norm(vr)
    Y = Y0 + eps*vr
    x = 0.0; step = xdir*abs(h)
    traj=[Y.copy()]; xs=[0.0]
    for i in range(nmax):
        Y = resolved_rk4(Y, step); x += step
        if not np.all(np.isfinite(Y)): break
        if abs(Y[1])>1e7 or abs(Y[0])>50:
            traj.append(Y.copy()); xs.append(x); break
        if i%100==0: traj.append(Y.copy()); xs.append(x)
        if abs(x)>xspan: break
    return np.array(xs), np.array(traj)

if __name__=="__main__":
    import sys
    V0=-0.25
    if len(sys.argv)>1: V0=float(sys.argv[1])
    Y0=D.sonic_point(V0)
    print(f"V0={V0} sonic=(A,N,om,V)={Y0.round(5)}")
    for eps in (+1e-5,-1e-5):
        for xdir in (-1,+1):
            xs,tr = run(V0,eps,xdir)
            if len(tr)<2: print(f" eps={eps:+.0e} xdir={xdir:+d}: died"); continue
            A,N,om,V=tr[-1]
            vz=int(np.sum(np.diff(np.sign(tr[:,3][tr[:,3]!=0]))!=0))
            print(f" eps={eps:+.0e} xdir={xdir:+d}: x_end={xs[-1]:+.3f} A={A:.4f} N={N:.4g} om={om:.5f} V={V:.4f}  Amin={tr[:,0].min():.4f} Vzeros={vz}")
