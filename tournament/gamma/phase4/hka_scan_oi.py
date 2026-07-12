# Scan the central-density parameter oi (N_inf = gauge); launch from the regular center, integrate
# outward, and see how the trajectory approaches the sonic locus Dson=0. The EC value of oi is where
# the trajectory reaches the sonic point ANALYTICALLY (numerators numOm,numV -> 0 together with Dson,
# so om',V' stay finite). We record the closest approach to Dson=0 (min|D|) and the value of V there;
# the L'Hopital residual (a*RHS_e - c*RHS_d) at that point should -> 0 at the critical oi.
import numpy as np, math
import hka_center as C

G=C.G; GM1=C.GM1
def lhop_resid(A,N,om,V):
    """Row-proportionality residual (4.6): a*RHS_e - c*RHS_d, with a=(1+NV)/om, c=(g-1)(N+V)/om.
    Times om: (1+NV)*RHS_e - (g-1)(N+V)*RHS_d. Vanishes at an analytic sonic passage."""
    g=G
    RHS_d=(3*(2-g)/2)*N*V-((2+g)/2)*A*N*V+(2-g)*N*V*om
    RHS_e=(2-g)*(g-1)*N*om+((7*g-6)/2)*N+((2-3*g)/2)*A*N
    return (1+N*V)*RHS_e-(g-1)*(N+V)*RHS_d

def probe(N_inf, oi, x0=-12.0, h=2e-4, nmax=200000):
    z0=math.exp(x0)
    Y4=C.center_seed(N_inf,oi,z0); Y3=Y4[[1,2,3]]
    x=x0; Dprev=C.Dson(Y3[0],Y3[2])
    minabsD=abs(Dprev); Y_at_min=Y3.copy(); x_at_min=x
    ncross=0; Vcross=None; xcross=None
    for i in range(nmax):
        Yn=C.rk4_3(Y3,abs(h)); xn=x+abs(h)
        if not np.all(np.isfinite(Yn)):
            break
        Dn=C.Dson(Yn[0],Yn[2])
        if abs(Dn)<minabsD:
            minabsD=abs(Dn); Y_at_min=Yn.copy(); x_at_min=xn
        if Dprev*Dn<0:
            ncross+=1; frac=Dprev/(Dprev-Dn); Yc=Y3+frac*(Yn-Y3)
            if ncross==1: Vcross=Yc[2]; xcross=x+frac*abs(h); Ycross=Yc.copy()
        Y3=Yn; x=xn; Dprev=Dn
        if abs(Y3[2])>3 or Y3[1]<-2 or Y3[1]>80 or Y3[0]>1e9 or Y3[0]<1e-5:
            break
        if x>3: break
    N,om,V=Y_at_min; A=C.A_of(N,om,V)
    lh=lhop_resid(A,N,om,V)
    out=dict(x_end=x,minabsD=minabsD,x_at_min=x_at_min,Vmin=V,Nmin=N,ommin=om,Amin=A,lhop=lh,ncross=ncross,Yend=Y3)
    if ncross>=1: out.update(Vcross=Vcross,xcross=xcross)
    return out

if __name__=="__main__":
    import sys,time
    N_inf=1.0
    if len(sys.argv)>1: N_inf=float(sys.argv[1])
    lo,hi,step=(0.1,6.0,0.1)
    if len(sys.argv)>4: lo,hi,step=float(sys.argv[2]),float(sys.argv[3]),float(sys.argv[4])
    t0=time.time()
    print(f"N_inf={N_inf} Vi={C.M_CENTER/N_inf:.4f}. Closest approach to sonic locus + L'Hopital resid:")
    print(f"{'oi':>8} {'min|D|':>10} {'x@min':>8} {'V@min':>9} {'N@min':>8} {'om@min':>8} {'A@min':>8} {'lhop':>11} {'ncr':>3} {'xend':>7}")
    for oi in np.arange(lo,hi+1e-9,step):
        r=probe(N_inf,round(oi,4))
        print(f"{oi:8.4f} {r['minabsD']:10.3e} {r['x_at_min']:8.3f} {r['Vmin']:9.4f} {r['Nmin']:8.4f} {r['ommin']:8.4f} {r['Amin']:8.4f} {r['lhop']:11.3e} {r['ncross']:3d} {r['x_end']:7.2f}")
    print(f"# {time.time()-t0:.1f}s")
