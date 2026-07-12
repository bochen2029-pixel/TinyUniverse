# hka_beta_solve.py — eigenvalue via the INTERSECTION rank-drop test, on the corrected operator.
#
# Geometry (4D solution space at match point x_mid):
#   U = center-regular subspace  = the 2 lambda=0 modes of L(xc)  [gauge + physical], integrated xc->xmid.
#   W = sonic-analytic subspace  = the 3 analytic Frobenius modes, integrated (xs+tson)->xmid.
# Generic dim(U cap W)=2+3-4=1 (the GAUGE mode, which is both). At an eigenvalue a 2nd (physical) mode
# joins => dim=2 => the stacked 4x5 matrix [U|W] drops rank 4->3 => its 4th singular value sigma4 -> 0.
# No hand gauge-removal. Columns unit-normalized => scale-invariant. HONEST: no target used.
# Eq #s per HKA_beta_equations.md.
import numpy as np, sys
import hka_beta4 as B, hka_pert_core as PC, hka_frobenius as FR
B.bg(); B.bg_path()
xs=B.bg()['xs']; bgser=B.bg()

def L_of_x(x,k): return PC.Lnum(B.bg_state(x),complex(k))

def rk4(Y,x0,x1,k,h=5e-4):
    n=max(1,int(round(abs(x1-x0)/h))); hh=(x1-x0)/n; x=x0; Y=Y.astype(complex)
    for _ in range(n):
        k1=L_of_x(x,k).dot(Y); k2=L_of_x(x+hh/2,k).dot(Y+hh/2*k1)
        k3=L_of_x(x+hh/2,k).dot(Y+hh/2*k2); k4=L_of_x(x+hh,k).dot(Y+hh*k3)
        Y=Y+hh/6*(k1+2*k2+2*k3+k4); x+=hh
        nb=np.linalg.norm(Y)
        if nb>1e150: Y=Y/nb
        if not np.all(np.isfinite(Y)): return None
    return Y

def center_modes(k, xc=-13.0):
    """The 2 regular (lambda~0) modes of L(xc)."""
    L=PC.Lnum(B.bg_state(xc),complex(k))
    ev,Vr=np.linalg.eig(L); idx=np.argsort(np.abs(ev.real))
    return [Vr[:,idx[0]], Vr[:,idx[1]]], np.sort(np.abs(ev.real))

def svals(k, xmid=None, xc=-13.0, tson=-0.01, order=12, h=5e-4):
    k=complex(k)
    if xmid is None: xmid=xs-1.0
    cols=[]
    cm,_=center_modes(k,xc)
    for u in cm:
        c=rk4(u,xc,xmid,k,h)
        if c is None: return None
        cols.append(c/np.linalg.norm(c))
    modes,R=FR.analytic_modes(k,bgser,order=order)
    if len(modes)<3: return None
    for a in modes:
        s=rk4(FR.eval_mode(a,tson),xs+tson,xmid,k,h)
        if s is None: return None
        cols.append(s/np.linalg.norm(s))
    M=np.column_stack(cols)                       # 4x5
    return np.linalg.svd(M,compute_uv=False)      # 4 singular values; sv[3]=sigma4 -> 0 at eigenvalue

if __name__=="__main__":
    grid=[float(a) for a in sys.argv[1:]] or list(np.round(np.arange(2.0,4.01,0.2),3))
    print(f"intersection rank-drop test: xs={xs:.5f}, xmid={xs-1.0:.5f}")
    print(f"{'kappa':>10} {'sigma4':>12} {'sigma3':>12}")
    for k in grid:
        sv=svals(k)
        if sv is None: print(f"{k:10.5f}  fail"); continue
        print(f"{k:10.5f} {sv[3]:12.4e} {sv[2]:12.4e}")
