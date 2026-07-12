# Scan V0 with the stable constraint-reduced 3D flow; launch center-ward off the sonic point,
# integrate in -x, and record where the trajectory ends up. Look for the EC critical V0 where it
# lands on the regular center (A->1, V->0, N->+inf, NV->-2/(3g)=-0.5).
import numpy as np, math
import hka_reduced as R

def probe(V0, es, h=1e-4, nmax=300000, launch=1e-6):
    Y0,J,w,Vc=R.eigdir3(V0)
    order=np.argsort(-w.real); vec=Vc[:,order[0]]; vr=np.real(vec); vr/=np.linalg.norm(vr)
    Y=Y0+es*launch*vr; x=0.0
    Vprev=Y[2]; vz=0; Amin=R.A_of(*Y); Nmax=Y[0]; Vmin=Y[2]; Vmax=Y[2]
    # record A and V at the moment N first crosses thresholds, to define a monotonic miss functional
    A_at=[];
    for i in range(nmax):
        Y=R.rk4_3(Y,-abs(h)); x-=abs(h)
        if not np.all(np.isfinite(Y)): break
        N,om,V=Y; A=R.A_of(N,om,V)
        Amin=min(Amin,A); Nmax=max(Nmax,N); Vmin=min(Vmin,V); Vmax=max(Vmax,V)
        if V*Vprev<0: vz+=1
        Vprev=V
        if N>1e6 or abs(A)>1e3 or abs(V)>10 or om>1e2 or om<-1e2:
            break
        if x<-22: break
    N,om,V=Y; A=R.A_of(N,om,V)
    return dict(x=x,A=A,N=N,om=om,V=V,vz=vz,Amin=Amin,Nmax=Nmax,Vmin=Vmin,Vmax=Vmax,M=N*V)

if __name__=="__main__":
    import sys,time
    es=+1
    if len(sys.argv)>1: es=int(sys.argv[1])
    t0=time.time()
    print(f"launch es={es}")
    print(f"{'V0':>8} {'x_end':>8} {'A':>10} {'N':>10} {'om':>9} {'V':>9} {'Amin':>8} {'Nmax':>9} {'Vz':>3} {'M':>9}")
    for V0 in np.arange(-0.55,-0.04,0.01):
        r=probe(round(V0,4),es)
        print(f"{V0:8.3f} {r['x']:8.3f} {r['A']:10.4g} {r['N']:10.4g} {r['om']:9.4f} {r['V']:9.4f} "
              f"{r['Amin']:8.4f} {r['Nmax']:9.4g} {r['vz']:3d} {r['M']:9.3g}")
    print(f"# {time.time()-t0:.1f}s")
