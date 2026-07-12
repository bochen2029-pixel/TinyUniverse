# Examine the trajectory's approach to the sonic locus near oi*~0.375. Does it graze (analytic) or
# cross transversally? Use Dson(x) profile and a robust shooting functional = signed MIN of Dson
# (closest approach). For the EC solution the trajectory is TANGENT to Dson=0 (min Dson -> 0 from one
# side). This min-Dson functional is far less noisy than crossing-point lhop.
import numpy as np, math
import hka_center as C
G=C.G; GM1=C.GM1; CS=C.CS

def traj_minD(N_inf, oi, x0=-12.0, h=5e-5, nmax=400000):
    """Return the signed extremum of Dson nearest to 0 along the outward trajectory (before break)."""
    z0=math.exp(x0); Y4=C.center_seed(N_inf,oi,z0); Y3=Y4[[1,2,3]]
    x=x0; D=C.Dson(Y3[0],Y3[2])
    # Dson at center is large NEGATIVE (N huge). It rises toward 0. Track its max (closest to 0 from below)
    Dmax=D; xatmax=x; Yatmax=Y3.copy()
    crossed=False
    for i in range(nmax):
        Yn=C.rk4_3(Y3,abs(h)); xn=x+abs(h)
        if not np.all(np.isfinite(Yn)): break
        Dn=C.Dson(Yn[0],Yn[2])
        if Dn>Dmax:
            Dmax=Dn; xatmax=xn; Yatmax=Yn.copy()
        if Dn>0:            # crossed to positive side (transversal crossing)
            crossed=True
            Y3=Yn; x=xn; D=Dn
            break
        Y3=Yn; x=xn; D=Dn
        if abs(Y3[2])>2 or Y3[1]<-1 or Y3[1]>60 or Y3[0]<1e-5: break
        if x>2: break
    N,om,V=Yatmax; A=C.A_of(N,om,V)
    return dict(Dmax=Dmax, x=xatmax, N=N, om=om, V=V, A=A, crossed=crossed)

if __name__=="__main__":
    import sys,time
    N_inf=1.0; t0=time.time()
    print("signed closest-approach of Dson (from below) vs oi; sign(Dmax) flips at EC solution (tangency):")
    print(f"{'oi':>8} {'Dmax':>12} {'x@max':>8} {'V@max':>9} {'N@max':>8} {'om@max':>8} {'A@max':>8} {'crossed':>7}")
    for oi in np.arange(0.30,0.46,0.005):
        r=traj_minD(N_inf,round(oi,4))
        print(f"{oi:8.4f} {r['Dmax']:12.4e} {r['x']:8.3f} {r['V']:9.5f} {r['N']:8.5f} {r['om']:8.5f} {r['A']:8.5f} {str(r['crossed']):>7}")
    print(f"# {time.time()-t0:.1f}s")
