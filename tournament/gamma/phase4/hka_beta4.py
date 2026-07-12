# hka_beta4.py — Stage-B eigenvalue: two-sided match det[ v_center | s1 | s2 | s3 ] = 0.
#
# (This is the prior-work-verified NON-DEGENERATE formulation, now on the CORRECT HKA background.)
#  - v_center: the PHYSICAL 1D center-regular mode = the lambda=0 (bounded, Abar_p=0) subspace of L_c
#    with the GAUGE mode (5.20) projected OUT. The gauge-mode direction at the center is the
#    background-log-derivative vector g=((lnA)'_x,(lnN)'_x,(lnom)'_x,V'_x) ~ (0,-1,2,0). Integrated
#    center -> match point x_m along the true background.
#  - s1,s2,s3: the 3 sonic ANALYTIC Frobenius modes (ker(R)), evaluated at x_m via their Taylor series.
#  - Delta(kappa) = det[ v_center | s1 | s2 | s3 ], each column UNIT-normalized (scale-invariant).
#    Zero <=> the center-regular mode lies in the sonic-analytic span <=> eigenvalue.
#
# Gauge modes kappa~0.35699 (sonic gauge) and kappa=1 (origin gauge) are DISCARDED (md fn15).
# HONEST: no target used in Delta. Scan prints |Delta|. Eq #s per HKA_beta_equations.md.
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

def gauge_dir(xc, N_inf=1.0, oi=0.375):
    """Gauge-mode (5.20) direction at the center: background log-derivative vector."""
    z=math.exp(xc); N=N_inf/z; om=oi*z*z; V=E.M_CENTER/N_inf*z; A=E.A_of(N,om,V)
    omx,Vx=E._fluid_slopes(A,N,om,V)
    Ax=A*(1-A+2*om*(1+(1/3.0)*V*V)/(1-V*V)); Nx=N*(-2+A-(2-4/3.0)*om)
    g=np.array([Ax/A, Nx/N, omx/om, Vx],complex)
    return g/np.linalg.norm(g)

def v_center(kappa, xc=-13.0):
    """Physical 1D center-regular vector: lambda=0 span minus gauge direction."""
    z=math.exp(xc); N=1.0/z; om=0.375*z*z; V=E.M_CENTER*z
    L=PC.Lnum(PC.bg_fields(N,om,V),complex(kappa))
    ev,Vr=np.linalg.eig(L); idx=np.argsort(np.abs(ev.real))
    span=np.column_stack([Vr[:,idx[0]],Vr[:,idx[1]]])   # 2D lambda~0
    Q,_=np.linalg.qr(span)
    g=gauge_dir(xc)
    gp=Q.dot(Q.conj().T.dot(g))                          # gauge projected into the span
    gp=gp/np.linalg.norm(gp)
    # physical = component of span orthogonal to gp: take Q columns, remove gp
    v=Q[:,0]-np.vdot(gp,Q[:,0])*gp
    if np.linalg.norm(v)<1e-6: v=Q[:,1]-np.vdot(gp,Q[:,1])*gp
    return v/np.linalg.norm(v)

def integ(vec, kappa, x0, xm, h=1e-3):
    Y=vec.astype(complex); x=x0; n=int(round((xm-x0)/h)); hh=(xm-x0)/n
    for i in range(n):
        k1=L_of_x(x,kappa).dot(Y); k2=L_of_x(x+0.5*hh,kappa).dot(Y+0.5*hh*k1)
        k3=L_of_x(x+0.5*hh,kappa).dot(Y+0.5*hh*k2); k4=L_of_x(x+hh,kappa).dot(Y+hh*k3)
        Y=Y+(hh/6)*(k1+2*k2+2*k3+k4); x+=hh
        nb=np.linalg.norm(Y)
        if nb>1e120: Y=Y/nb
        if not np.all(np.isfinite(Y)): return None
    return Y/np.linalg.norm(Y)

def eval_series(a,t):
    s=np.zeros(4,complex); tp=1.0
    for n in range(len(a)): s+=a[n]*tp; tp*=t
    return s

def Delta(kappa, xm=None, xc=-13.0, h=1e-3, order=8, return_rej=True):
    """Match residual. c = unit center-regular vector (integrated to x_m). Q = ORTHONORMAL basis of the
    3D sonic-analytic subspace at x_m (QR of the 3 Frobenius modes -> cures the near-parallel
    degeneracy). Residual = rejection ||c - Q Q^H c|| (scale-invariant, in [0,1]); ALSO det[c|Q]
    (== rejection up to unit phase since Q orthonormal). Zero <=> c in span(Q) <=> eigenvalue."""
    kappa=complex(kappa); xs=bg()['xs']
    if xm is None: xm=xs-0.03
    vc=v_center(kappa,xc)
    c=integ(vc,kappa,xc,xm,h=h)
    if c is None: return None
    modes,R=FR.analytic_modes(kappa,bg(),order=order)
    if len(modes)<3: return None
    tm=xm-xs
    scols=[eval_series(a,tm) for a in modes]
    Q,_=np.linalg.qr(np.column_stack(scols))      # 4x3 orthonormal analytic-subspace basis
    Q=Q[:,:3]
    rej=c-Q.dot(Q.conj().T.dot(c))
    det=np.linalg.det(np.column_stack([c,Q]))     # scale-invariant (c unit, Q orthonormal)
    if return_rej:
        return det, np.linalg.norm(rej)
    return det

if __name__=="__main__":
    import sys,time
    t0=time.time(); bg(); bg_path()
    xm=_XS-0.03
    print(f"sonic x_s={_XS:.5f}, x_m={xm:.5f}. det[v_center|s1|s2|s3] (columns unit-normalized):")
    print(f"{'kappa':>8} {'Re Delta':>13} {'Im Delta':>13} {'|Delta|':>11}")
    rows=[]
    for kap in np.arange(0.3,6.01,0.1):
        d=Delta(kap)
        if d is None: print(f"{kap:8.2f}  (fail)"); continue
        rows.append((kap,d)); print(f"{kap:8.2f} {d.real:13.4e} {d.imag:13.4e} {abs(d):11.3e}")
    print("\n|Delta| local minima:")
    for i in range(1,len(rows)-1):
        if abs(rows[i][1])<abs(rows[i-1][1]) and abs(rows[i][1])<abs(rows[i+1][1]):
            print(f"   kappa~{rows[i][0]:.2f}  |Delta|={abs(rows[i][1]):.3e}")
    print("Re(Delta) sign changes:")
    for i in range(1,len(rows)):
        if rows[i-1][1].real*rows[i][1].real<0:
            print(f"   {rows[i-1][0]:.2f} -> {rows[i][0]:.2f}  |Delta|~{min(abs(rows[i-1][1]),abs(rows[i][1])):.2e}")
    print(f"# {time.time()-t0:.1f}s")
