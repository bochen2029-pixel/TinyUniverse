# hka_beta3.py — Stage-B eigenvalue via a clean 3x3 solvability determinant (sonic->center).
#
# The 3 sonic ANALYTIC modes (ker(R)) are integrated to the center. A physical eigenperturbation is a
# combination c=(c0,c1,c2) of them satisfying the THREE physical requirements:
#   (I)   gauge Nbar_p(x_s)=0            -> row [Nbar comp of a_0 for each mode]
#   (II)  no e^{-2x} growing mode at center  -> row [e^{-2x} amplitude of each integrated mode]
#   (III) no e^{-x}  growing mode at center  -> row [e^{-x}  amplitude of each integrated mode]
# A nontrivial c exists iff the 3x3 matrix M(kappa)=[I;II;III] is SINGULAR: det M = 0. This is the
# eigenvalue condition. Gauge modes kappa~0.357 and kappa=1 are DISCARDED (md fn15).
#
# Rows are scaled to O(1) (each integrated mode normalized before reading amplitudes) so det M is a
# clean, well-conditioned function of kappa. HONEST: no target used; scan prints |det M|.
# Eq #s per HKA_beta_equations.md.
import numpy as np, math
import hka_pert_core as PC
import hka_frobenius as FR
import hka_pert_sonic as PS
import hka_ec as E

_BGSER=None;_XS=None
def bg():
    global _BGSER,_XS
    if _BGSER is None:
        _BGSER,_,_=PS.bg_series_near_sonic(order=8,dt=2e-3,npts=30); _XS=_BGSER['xs']
    return _BGSER

_BGPATH=None
def bg_path(N_inf=1.0, oi=0.375, xc=-16.0):
    global _BGPATH
    if _BGPATH is None:
        from scipy.integrate import solve_ivp
        z0=math.exp(xc); Y0=E.center_state3(N_inf,oi,z0); xs=bg()['xs']
        sol=solve_ivp(E.rhs3,[xc,xs-1e-7],Y0,method='DOP853',rtol=1e-12,atol=1e-14,dense_output=True)
        _BGPATH=(sol,xs,xc)
    return _BGPATH

def bg_state(x):
    sol,xs,xc=bg_path(); N,om,V=sol.sol(x); return PC.bg_fields(N,om,V)
def L_of_x(x,kappa): return PC.Lnum(bg_state(x),complex(kappa))

def eval_series(a,t):
    s=np.zeros(4,complex); tp=1.0
    for n in range(len(a)): s+=a[n]*tp; tp*=t
    return s

def integ(vec, kappa, x0, xc, h=2e-3):
    Y=vec.astype(complex); x=x0; n=int(round((x0-xc)/h)); hh=-(x0-xc)/n
    for i in range(n):
        k1=L_of_x(x,kappa).dot(Y); k2=L_of_x(x+0.5*hh,kappa).dot(Y+0.5*hh*k1)
        k3=L_of_x(x+0.5*hh,kappa).dot(Y+0.5*hh*k2); k4=L_of_x(x+hh,kappa).dot(Y+hh*k3)
        Y=Y+(hh/6)*(k1+2*k2+2*k3+k4); x+=hh
        if not np.all(np.isfinite(Y)): return None,x
        # rescale to avoid overflow but keep amplitude ratios (record log-scale not needed for det=0)
        nb=np.max(np.abs(Y))
        if nb>1e120: Y=Y/nb
    return Y,x

def center_dual(kappa, xc):
    z=math.exp(xc); N=1.0/z; om=0.375*z*z; V=E.M_CENTER*z
    L=PC.Lnum(PC.bg_fields(N,om,V),complex(kappa))
    ev,Vr=np.linalg.eig(L); Vinv=np.linalg.inv(Vr); idx=np.argsort(ev.real)
    return ev,idx,Vinv

def detM(kappa, x0=None, xc=-14.0, h=2e-3, order=8):
    kappa=complex(kappa); xs=bg()['xs']
    if x0 is None: x0=xs-0.02
    modes,R=FR.analytic_modes(kappa, bg(), order=order)
    if len(modes)<3: return None
    a0=np.array([m[0] for m in modes]); Nbar0=a0[:,1]     # row I (gauge): Nbar comp of each a_0
    t0=x0-xs
    rows_growing=[[],[]]     # [e^{-2x} amps],[e^{-x} amps]
    ev=idx=Vinv=None
    for a in modes:
        v0=eval_series(a,t0)
        Yc,xend=integ(v0,kappa,x0,xc,h=h)
        if Yc is None: return None
        Yc=Yc/np.linalg.norm(Yc)              # normalize -> O(1) amplitudes
        if ev is None: ev,idx,Vinv=center_dual(kappa,xend)
        amp=Vinv.dot(Yc)
        rows_growing[0].append(amp[idx[0]])   # e^{-2x} (lambda=-2)
        rows_growing[1].append(amp[idx[1]])   # e^{-x}  (lambda=-1)
    M=np.array([Nbar0/np.linalg.norm(Nbar0), rows_growing[0], rows_growing[1]])
    return np.linalg.det(M), M, ev[idx]

if __name__=="__main__":
    import sys,time
    t0=time.time(); bg(); bg_path()
    xc=-14.0
    if len(sys.argv)>1: xc=float(sys.argv[1])
    print(f"sonic x_s={_XS:.5f}, xc={xc}. 3x3 solvability det (gauge + kill e^-2x + kill e^-x):")
    print(f"{'kappa':>8} {'Re detM':>13} {'Im detM':>13} {'|detM|':>11}")
    rows=[]
    for kap in np.arange(0.3,6.01,0.1):
        out=detM(kap,xc=xc)
        if out is None: print(f"{kap:8.2f}  (fail)"); continue
        d,M,cev=out; rows.append((kap,d))
        print(f"{kap:8.2f} {d.real:13.4e} {d.imag:13.4e} {abs(d):11.3e}")
    print("\n|detM| local minima (candidate real eigenvalues):")
    for i in range(1,len(rows)-1):
        if abs(rows[i][1])<abs(rows[i-1][1]) and abs(rows[i][1])<abs(rows[i+1][1]):
            print(f"   kappa~{rows[i][0]:.2f}  |detM|={abs(rows[i][1]):.3e}")
    print("Re(detM) sign changes:")
    for i in range(1,len(rows)):
        if rows[i-1][1].real*rows[i][1].real<0:
            print(f"   between {rows[i-1][0]:.2f} and {rows[i][0]:.2f}  (|detM|~{abs(rows[i][1]):.2e})")
    print(f"# {time.time()-t0:.1f}s")
