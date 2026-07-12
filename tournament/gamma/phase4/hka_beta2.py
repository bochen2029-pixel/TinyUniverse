# hka_beta2.py — Stage-B eigenvalue via sonic->center shoot (md's prescribed direction).
#
# Launch the sonic ANALYTIC subspace (Frobenius, 3D ker(R)) reduced by gauge Nbar_p(x_s)=0 -> 2D.
# Integrate sonic -> center (decreasing x). At the center, decompose onto the 4 indicial modes of L_c:
#   lambda=0 (x2, regular), lambda=-1 (e^{-x}), lambda=-2 (e^{-2x}).  [these -> +inf as x->-inf].
# Regularity requires the GROWING modes (e^{-x}, e^{-2x}) to be ABSENT. For the 2D sonic subspace this
# is 2 conditions -> eigenvalues are kappa where the 2x2 matrix
#     B(kappa) = [ [ampl e^{-2x} of mode1, ampl e^{-2x} of mode2],
#                  [ampl e^{-x}  of mode1, ampl e^{-x}  of mode2] ]
# is SINGULAR: det B = 0. (If (5.15) auto-kills e^{-x}, that row ~0 and det~0 trivially -> then use
# only the e^{-2x} row with the gauge giving a 1D physical mode; we test both.)
#
# HONEST: no target used. Print |det B| and the growing-mode amplitudes vs kappa.
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

def sonic_gauge_modes(kappa, order=8):
    """2 analytic modes (Taylor series lists) spanning the gauge-fixed (Nbar_p(x_s)=0) subspace."""
    modes,R=FR.analytic_modes(kappa, bg(), order=order)
    if len(modes)<3: return None
    a0=np.array([m[0] for m in modes]); Nbar0=a0[:,1]
    u=Nbar0/np.linalg.norm(Nbar0); basis=[]
    for e in np.eye(3):
        v=e-np.vdot(u,e)*u
        for b in basis: v=v-np.vdot(b,v)*b
        if np.linalg.norm(v)>1e-8: basis.append(v/np.linalg.norm(v))
        if len(basis)==2: break
    combos=[]
    for c in basis:
        combos.append([sum(c[i]*modes[i][n] for i in range(3)) for n in range(order+1)])
    return combos   # 2 modes, each a list of Taylor coeffs

def eval_series(a, t):
    s=np.zeros(4,complex); tp=1.0
    for n in range(len(a)): s+=a[n]*tp; tp*=t
    return s

def integrate_to_center(vec, kappa, x_start, x_c, h=2e-3):
    """RK4 a single perturbation vector from x_start (near sonic) to x_c (deep center)."""
    Y=vec.astype(complex); x=x_start
    n=int(round((x_start-x_c)/h)); hh=-(x_start-x_c)/n
    for i in range(n):
        k1=L_of_x(x,kappa).dot(Y)
        k2=L_of_x(x+0.5*hh,kappa).dot(Y+0.5*hh*k1)
        k3=L_of_x(x+0.5*hh,kappa).dot(Y+0.5*hh*k2)
        k4=L_of_x(x+hh,kappa).dot(Y+hh*k3)
        Y=Y+(hh/6)*(k1+2*k2+2*k3+k4); x+=hh
        nrm=np.linalg.norm(Y)
        if nrm>1e150: return None,x
        if not np.all(np.isfinite(Y)): return None,x
    return Y,x

def center_left_modes(kappa, xc):
    """Left eigenvectors (dual) of L_c for the e^{-x}(lambda=-1) and e^{-2x}(lambda=-2) modes, to read
    off those amplitudes from a state at xc."""
    z=math.exp(xc); N=1.0/z; om=0.375*z*z; V=E.M_CENTER*z
    L=PC.Lnum(PC.bg_fields(N,om,V),complex(kappa))
    ev,Vr=np.linalg.eig(L)
    # right eigvecs Vr; left = rows of inv(Vr). amplitude of mode i = (Vr^{-1} Y)_i
    Vinv=np.linalg.inv(Vr)
    idx=np.argsort(ev.real)   # -2, -1, 0, 0
    return ev,idx,Vinv

def detB(kappa, x_start=None, xc=-15.0, h=2e-3, order=8):
    """2x2 growing-mode amplitude matrix determinant. Zero at eigenvalues."""
    kappa=complex(kappa); xs=bg()['xs']
    if x_start is None: x_start=xs-0.02
    combos=sonic_gauge_modes(kappa,order=order)
    if combos is None: return None
    t0=x_start-xs
    amps=[]
    ev=None;idx=None;Vinv=None
    cols=[]
    for a in combos:
        v0=eval_series(a,t0)
        Yc,xend=integrate_to_center(v0,kappa,x_start,xc,h=h)
        if Yc is None: return None
        if ev is None:
            ev,idx,Vinv=center_left_modes(kappa,xend)
        amp=Vinv.dot(Yc)            # amplitudes in the L_c eigenbasis at xend
        # growing modes: lambda=-2 (idx[0]) and lambda=-1 (idx[1])
        cols.append([amp[idx[0]], amp[idx[1]]])
    B=np.array(cols).T             # 2x2: rows=(e^{-2x},e^{-x}), cols=modes
    return np.linalg.det(B), B, ev[idx]

if __name__=="__main__":
    import sys,time
    t0=time.time(); bg(); bg_path()
    print(f"sonic x_s={_XS:.5f}. sonic->center shoot; det of 2x2 growing-mode amplitude matrix:")
    print(f"{'kappa':>8} {'Re detB':>13} {'Im detB':>13} {'|detB|':>11}   center eigenvalues")
    rows=[]
    for kap in np.arange(0.3,6.01,0.1):
        out=detB(kap, xc=-15.0)
        if out is None: print(f"{kap:8.2f}  (fail)"); continue
        d,B,cev=out; rows.append((kap,d))
        print(f"{kap:8.2f} {d.real:13.4e} {d.imag:13.4e} {abs(d):11.3e}   {cev.round(2)}")
    print("\n|detB| local minima:")
    for i in range(1,len(rows)-1):
        if abs(rows[i][1])<abs(rows[i-1][1]) and abs(rows[i][1])<abs(rows[i+1][1]):
            print(f"   kappa~{rows[i][0]:.2f}  |detB|={abs(rows[i][1]):.3e}")
    print(f"# {time.time()-t0:.1f}s")
