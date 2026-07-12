# Scan V0: launch inward (toward center) off the sonic point and record the state where the
# trajectory either (a) reaches the regular center (A->1,V->0,N large) or (b) diverges. Look for
# the sign change / bracket that isolates the Evans-Coleman critical V0.
import numpy as np, math
import hka_desing as D

def resolved_rk4(Y,h):
    k1=D.resolved_rhs(Y);k2=D.resolved_rhs(Y+0.5*h*k1);k3=D.resolved_rhs(Y+0.5*h*k2);k4=D.resolved_rhs(Y+h*k3)
    return Y+(h/6)*(k1+2*k2+2*k3+k4)

def inward(V0, eps=-1e-5, h=1e-4, nmax=400000):
    """Launch inward along the desing eigendirection; integrate resolved flow in -x.
    Return diagnostics: (status, x, A, N, om, V, Amin, Amax, Vzeros).
    status: 'center' if reaches A~1,V~0 with N>1e4; 'diverge' otherwise."""
    Y0,J,w,Vc = D.sonic_jacobian(V0)
    order=np.argsort(-w.real); vec=Vc[:,order[0]]; vr=np.real(vec); vr/=np.linalg.norm(vr)
    # choose eps sign so trajectory goes toward center; test both, keep the one with N increasing
    best=None
    for es in (eps, -eps):
        Y=Y0+es*vr; x=0.0; A0hist=[]; ok=True; Vprev=Y[3]; vz=0; Amin=Y[0]; Amax=Y[0]
        started_center=False
        for i in range(nmax):
            Y=resolved_rk4(Y,-abs(h)); x-=abs(h)
            if not np.all(np.isfinite(Y)): ok=False; break
            A,N,om,V=Y
            Amin=min(Amin,A); Amax=max(Amax,A)
            if V*Vprev<0: vz+=1
            Vprev=V
            if abs(N)>1e6 or abs(A)>1e3 or abs(V)>1e3:
                break
            if N>50 and abs(V)<0.05 and abs(A-1)<0.2:
                started_center=True
        A,N,om,V=Y
        # score: reached center if N large and A near 1 and V near 0
        reached = (N>1e3 and abs(A-1)<0.5 and abs(V)<0.3)
        entry=('center' if reached else 'diverge', x, A, N, om, V, Amin, Amax, vz, es)
        # prefer the branch that got N largest with A closest to 1
        key = (reached, -abs(A-1) if reached else -1e9, N)
        if best is None or key>best[0]:
            best=(key, entry)
    return best[1]

if __name__=="__main__":
    print(f"{'V0':>8} {'status':>8} {'x':>8} {'A':>10} {'N':>10} {'om':>9} {'V':>9} {'Amin':>8} {'Vz':>3} {'es':>7}")
    for V0 in np.arange(-0.55,-0.04,0.01):
        try:
            r=inward(round(V0,4))
        except Exception as e:
            print(f"{V0:8.3f}  err {e}"); continue
        st,x,A,N,om,V,Amin,Amax,vz,es=r
        print(f"{V0:8.3f} {st:>8} {x:8.3f} {A:10.4g} {N:10.4g} {om:9.4f} {V:9.4f} {Amin:8.4f} {vz:3d} {es:+.0e}")
