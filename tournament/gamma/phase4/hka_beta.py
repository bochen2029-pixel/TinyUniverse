# hka_beta.py — Stage-B eigenvalue shoot for the radiation-fluid critical exponent beta=1/Re kappa.
#
# Two-sided linear match on the CORRECT HKA background (Stage A: EC solution, sonic at (3/2,2/sqrt3,
# 3/4,-1/sqrt3), oi*=3/8):
#   - CENTER side: 2D regular subspace = the two lambda=0 modes of L_c (both have Abar_p=0), integrated
#     center -> match point x_m (RK4). (The e^{-x}, e^{-2x} modes are excluded/killed.)
#   - SONIC side: analytic Frobenius modes (ker(R), 3D), REDUCED to 2D by the gauge Nbar_p(x_s)=0
#     (a_0 has zero Nbar component), integrated sonic -> x_m.
#   - EIGENVALUE: kappa where the 4 vectors [c1|c2|s1|s2] at x_m become linearly dependent:
#         Delta(kappa) = det[ c1 | c2 | s1 | s2 ] = 0.
#     Columns are normalized; QR-orthonormalize each 2D block for conditioning.
#
# The gauge modes (md fn15): kappa~0.357 (sonic gauge) and kappa=1 (origin gauge) are DISCARDED.
# The physical relevant mode is expected at kappa~2.81 -> beta~0.3558. We scan complex kappa, find the
# roots of Delta, exclude the gauge modes, and report the physical kappa + convergence.
#
# HONEST: no target is used in Delta. The scan prints |Delta| so structure is visible first.
# Eq #s per HKA_beta_equations.md.

import numpy as np, math
import hka_pert_core as PC
import hka_frobenius as FR
import hka_pert_sonic as PS
import hka_ec as E

G=4.0/3.0

# ---- background near sonic (Taylor series in t=x-x_s), computed once ----
_BGSER=None; _XS=None
def bg():
    global _BGSER,_XS
    if _BGSER is None:
        _BGSER,_,_=PS.bg_series_near_sonic(order=8,dt=2e-3,npts=30)
        _XS=_BGSER['xs']
    return _BGSER

# ---- center regular subspace, integrated to x_m ----
def center_seeds(kappa, xc=-13.0, N_inf=1.0, oi=0.375):
    """Two lambda=0 (regular, Abar_p=0) eigenvectors of L_c at deep center xc."""
    z=math.exp(xc); N=N_inf/z; om=oi*z*z; V=E.M_CENTER/N_inf*z
    fld=PC.bg_fields(N,om,V)
    L=PC.Lnum(fld,complex(kappa))
    ev,Vc=np.linalg.eig(L)
    idx=np.argsort(np.abs(ev.real))         # smallest |Re| = the lambda=0 modes
    seeds=[Vc[:,idx[0]], Vc[:,idx[1]]]
    return [s/np.linalg.norm(s) for s in seeds], xc

# We need the TRUE background (A,N,om,V)(x) along [xc, x_m], not just the center asymptotics.
# Provide it by integrating the EC background once and interpolating.
_BGPATH=None
def bg_path(N_inf=1.0, oi=0.375, xc=-13.0):
    global _BGPATH
    if _BGPATH is None:
        from scipy.integrate import solve_ivp
        z0=math.exp(xc); Y0=E.center_state3(N_inf,oi,z0)
        xs=bg()['xs']
        sol=solve_ivp(E.rhs3,[xc,xs-1e-6],Y0,method='DOP853',rtol=1e-12,atol=1e-14,dense_output=True)
        _BGPATH=(sol,xs)
    return _BGPATH

def bg_state(x):
    sol,xs=bg_path()
    N,om,V=sol.sol(x)
    return PC.bg_fields(N,om,V)

def L_of_x(x, kappa):
    return PC.Lnum(bg_state(x), complex(kappa))

def integrate_center(seeds, kappa, xc, x_m, h=1e-3):
    """RK4 the center-regular subspace (matrix of columns) from xc to x_m along the TRUE background."""
    M=np.column_stack(seeds).astype(complex)
    x=xc; n=int(round((x_m-xc)/h)); hh=(x_m-xc)/n
    for i in range(n):
        k1=L_of_x(x,kappa).dot(M)
        k2=L_of_x(x+0.5*hh,kappa).dot(M+0.5*hh*k1)
        k3=L_of_x(x+0.5*hh,kappa).dot(M+0.5*hh*k2)
        k4=L_of_x(x+hh,kappa).dot(M+hh*k3)
        M=M+(hh/6)*(k1+2*k2+2*k3+k4)
        x+=hh
        # periodic re-orthonormalize for conditioning
        if i%50==49:
            M,_=np.linalg.qr(M)
        if not np.all(np.isfinite(M)): return None
    Q,_=np.linalg.qr(M)
    return Q       # 4x2 orthonormal

def sonic_subspace(kappa, x_m, order=8):
    """2D analytic subspace at x_m: 3 Frobenius modes reduced by gauge Nbar_p(x_s)=0, evaluated at
    t_m=x_m-x_s via the Taylor series (x_m just inside sonic so series converges)."""
    modes,R=FR.analytic_modes(kappa, bg(), order=order)
    if len(modes)<3: return None
    xs=bg()['xs']; t_m=x_m-xs
    # gauge Nbar_p(x_s)=0: a_0 Nbar component = modes[i][0][1]. Build 2 combos with zero Nbar in a_0.
    a0=np.array([m[0] for m in modes])          # 3x4
    Nbar0=a0[:,1]                                # Nbar comp of each a_0
    # null space of Nbar0 (1x3) -> 2D combos
    # combos c (3-vector) with Nbar0.c=0
    u=Nbar0/np.linalg.norm(Nbar0)
    # two orthonormal vectors orthogonal to u in C^3
    basis=[]
    for e in np.eye(3):
        v=e-np.vdot(u,e)*u
        for b in basis: v=v-np.vdot(b,v)*b
        if np.linalg.norm(v)>1e-8:
            basis.append(v/np.linalg.norm(v))
        if len(basis)==2: break
    cols=[]
    for c in basis:
        # combined mode = sum_i c_i modes[i]; evaluate at t_m
        val=np.zeros(4,complex)
        for n in range(order+1):
            an=sum(c[i]*modes[i][n] for i in range(3))
            val+=an*(t_m**n)
        cols.append(val)
    Q,_=np.linalg.qr(np.column_stack(cols))
    return Q      # 4x2

def Delta(kappa, x_m=None, xc=-13.0, h=1e-3):
    """Match determinant det[ Qc | Qs ] (4x4). Zero at eigenvalues."""
    kappa=complex(kappa)
    xs=bg()['xs']
    if x_m is None: x_m=xs-0.03
    seeds,xc=center_seeds(kappa,xc)
    Qc=integrate_center(seeds,kappa,xc,x_m,h=h)
    if Qc is None: return None
    Qs=sonic_subspace(kappa,x_m)
    if Qs is None: return None
    Mmat=np.column_stack([Qc,Qs])
    return np.linalg.det(Mmat)

if __name__=="__main__":
    import sys,time
    t0=time.time()
    bg(); bg_path()
    print(f"sonic x_s={_XS:.5f}, match x_m={_XS-0.03:.5f}. Scanning real kappa (|Delta| structure):")
    print(f"{'kappa':>8} {'Re Delta':>13} {'Im Delta':>13} {'|Delta|':>11}")
    rows=[]
    for kap in np.arange(0.2,6.01,0.1):
        d=Delta(kap)
        if d is None:
            print(f"{kap:8.2f}   (fail)"); continue
        rows.append((kap,d)); print(f"{kap:8.2f} {d.real:13.4e} {d.imag:13.4e} {abs(d):11.3e}")
    print("\n|Delta| local minima (candidate real eigenvalues):")
    for i in range(1,len(rows)-1):
        if abs(rows[i][1])<abs(rows[i-1][1]) and abs(rows[i][1])<abs(rows[i+1][1]):
            print(f"   kappa~{rows[i][0]:.2f}  |Delta|={abs(rows[i][1]):.3e}")
    print("Re(Delta) sign changes:")
    for i in range(1,len(rows)):
        if rows[i-1][1].real*rows[i][1].real<0:
            print(f"   between kappa={rows[i-1][0]:.2f} and {rows[i][0]:.2f}")
    print(f"# {time.time()-t0:.1f}s")
