# hka_shoot_oi.py — shoot the central density oi so the center-launched trajectory reaches the sonic
# locus (Dson=0) ANALYTICALLY: the L'Hopital residual (4.6) must vanish AT the first sonic crossing.
# The coarse scan (hka_scan_oi.py) shows lhop(cross) flips sign between oi~0.3 and ~0.4 -> bracket it.
import numpy as np, math
import hka_center as C
G=C.G

def lhop(A,N,om,V):
    g=G
    RHS_d=(3*(2-g)/2)*N*V-((2+g)/2)*A*N*V+(2-g)*N*V*om
    RHS_e=(2-g)*(g-1)*N*om+((7*g-6)/2)*N+((2-3*g)/2)*A*N
    return (1+N*V)*RHS_e-(g-1)*(N+V)*RHS_d

def first_cross(N_inf, oi, x0=-12.0, h=1e-4, nmax=300000):
    """Integrate center->outward; at the FIRST Dson=0 crossing return (x,A,N,om,V,lhop). If it breaks
    before crossing, return status='break'."""
    z0=math.exp(x0); Y4=C.center_seed(N_inf,oi,z0); Y3=Y4[[1,2,3]]
    x=x0; Dprev=C.Dson(Y3[0],Y3[2])
    for i in range(nmax):
        Yn=C.rk4_3(Y3,abs(h)); xn=x+abs(h)
        if not np.all(np.isfinite(Yn)):
            return dict(status='nan',x=x,Y3=Y3)
        Dn=C.Dson(Yn[0],Yn[2])
        if Dprev<0 and Dn>=0 or (Dprev>0 and Dn<=0):
            frac=Dprev/(Dprev-Dn); Yc=Y3+frac*(Yn-Y3); xc=x+frac*abs(h)
            N,om,V=Yc; A=C.A_of(N,om,V)
            return dict(status='cross',x=xc,A=A,N=N,om=om,V=V,lhop=lhop(A,N,om,V),Dfrom=Dprev)
        Y3=Yn; x=xn; Dprev=Dn
        if abs(Y3[2])>2 or Y3[1]<-1 or Y3[1]>60 or Y3[0]>1e9 or Y3[0]<1e-5:
            return dict(status='break',x=x,Y3=Y3)
        if x>3: return dict(status='nocross',x=x,Y3=Y3)
    return dict(status='maxstep',x=x,Y3=Y3)

def F(oi, N_inf=1.0, **kw):
    r=first_cross(N_inf, oi, **kw)
    if r['status']!='cross':
        return None, r
    return r['lhop'], r

if __name__=="__main__":
    import sys,time
    N_inf=1.0
    t0=time.time()
    print("fine scan of lhop at first sonic crossing (bracket the analytic passage):")
    print(f"{'oi':>9} {'status':>7} {'xcross':>8} {'V':>9} {'N':>8} {'om':>8} {'A':>8} {'lhop':>12}")
    prev=None; bracket=[]
    for oi in np.arange(0.20,0.60,0.01):
        val,r=F(round(oi,4),N_inf)
        st=r['status']
        if st=='cross':
            print(f"{oi:9.4f} {st:>7} {r['x']:8.3f} {r['V']:9.4f} {r['N']:8.4f} {r['om']:8.4f} {r['A']:8.4f} {r['lhop']:12.4e}")
            if prev is not None and prev[1] is not None and val is not None and prev[1]*val<0:
                bracket.append((prev[0],oi))
            prev=(oi,val)
        else:
            print(f"{oi:9.4f} {st:>7} {r['x']:8.3f}")
            prev=(oi,None)
    print("brackets (lhop sign change):", bracket)
    print(f"# {time.time()-t0:.1f}s")
