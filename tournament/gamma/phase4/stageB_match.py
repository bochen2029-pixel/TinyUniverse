# stageB_match.py — robust two-point shoot (no divergent branch).
# LEFT: center-regular seed -> integrate outward to matching point x_m (just before sonic).
# RIGHT: sonic-ANALYTIC seed (Taylor series on the mu=0 subspace, c_na=0 built in) -> integrate
#        INWARD to x_m.  Eigenvalue: the two solution SPACES coincide.
#   Left space = span{ hp_L } (1-D, center-regular).
#   Right space = span of analytic solutions at sonic. dim = # independent analytic Taylor seeds.
# The match condition: hp_L is in the right analytic space. Residual = det[ hp_L | basis_R ]-type,
# i.e. the component of hp_L OUTSIDE the analytic space. Equivalent: project hp_L onto the
# NON-analytic left dual wL_n AFTER transporting; at eigenvalue that component (relative to analytic) =0.
#
# Simplest robust residual: integrate BOTH, and compute the 4x4 determinant of
#   [ hp_L(x_m) , d1, d2, d3 ]   where d1,d2,d3 = the 3 analytic solutions integrated inward.
# If hp_L lies in span{d1,d2,d3} the det=0 -> eigenvalue. This is the Wronskian match. No divergence.
import numpy as np, math, sympy as sp, pickle
import css_core as C, pert_core as P

S3=math.sqrt(3.0); A2STAR=0.250067905275; NC=1.0
_d=pickle.load(open("Lmat.pkl","rb"))
A0,N0,o0,V0=sp.symbols('A0 N0 om0 V0'); Nxs,Axs,oxs,Vxs=sp.symbols('N_x A_x om_x V_x'); ksym=sp.symbols('k')
Lsym=sp.sympify(_d['Lmat'])
_ARGS=(N0,A0,o0,V0,Nxs,Axs,oxs,Vxs,ksym)
fL_ij=[[sp.lambdify(_ARGS,Lsym[i,j],'numpy') for j in range(4)] for i in range(4)]
def Lmat_num(Y,slp,k):
    N,A,om,V=Y; Nb,Ab,ob,Vb=slp
    return np.array([[complex(fL_ij[i][j](N,A,om,V,Nb,Ab,ob,Vb,k)) for j in range(4)] for i in range(4)])

def bg_slopes(Y):
    N,A,om,V=Y; return (C.Nx(A,N,om,V),C.Ax(A,N,om,V),C.omx(A,N,om,V),C.Vx(A,N,om,V))

# ---- analytic Taylor seed at the sonic point: solve the recursion. Near sonic, background is the
# exact series; the perturbation ODE (x-x_s) hp' = [R + (x-x_s)R1 + ...] hp . Analytic solution
# hp = v0 + v1 (x-x_s) + ...  Plug in: matching orders,
#   order 0:  0 = R v0   -> v0 in ker(R) (mu=0 space, 3-dim)
#   order 1:  v1 = R v1 + R1 v0  -> (I-R) v1 = R1 v0 ... careful with the (x-x_s) hp' structure.
# Let t=x-x_s. ODE: hp' = L hp, L = R/t + L0 + L1 t + ...  Then t hp' = (R + L0 t + ...) hp.
# Write hp=sum a_n t^n. LHS t hp' = sum n a_n t^n. RHS = sum_m R? ... = R a_0 + (R a_1 + L0 a_0) t + ...
#   t^0:  0 = R a_0         => a_0 in ker R
#   t^1:  a_1 = R a_1 + L0 a_0  => (I-R) a_1 = L0 a_0
#   t^n:  n a_n = R a_n + sum_{j>=0} L_j a_{n-1-j}  => (nI - R) a_n = sum L_j a_{n-1-j}
# Build L0, L1 from the series of (L - R/t).
def build_analytic_seeds(k, order=6):
    t=sp.symbols('t')
    Ns=2/sp.sqrt(3)+(-2/sp.sqrt(3))*t+sp.Rational(1,2)*(2/sp.sqrt(3))*t**2
    As=sp.Rational(3,2)+3*t+sp.Rational(1,2)*33*t**2
    Os=sp.Rational(3,4)+sp.Rational(9,2)*t+sp.Rational(1,2)*sp.Rational(99,2)*t**2
    Vs=1/sp.sqrt(3)+(2/sp.sqrt(3))*t+sp.Rational(1,2)*(10*sp.sqrt(3)/3)*t**2
    sub={N0:Ns,A0:As,o0:Os,V0:Vs,Nxs:sp.diff(Ns,t),Axs:sp.diff(As,t),oxs:sp.diff(Os,t),Vxs:sp.diff(Vs,t)}
    Lt=Lsym.subs(sub).subs(ksym,k)
    # Laurent: tL = R + L0 t + L1 t^2 + ...
    tL=sp.Matrix(4,4,lambda i,j: sp.series(Lt[i,j]*t, t,0,order+1).removeO())
    def coeffM(p):
        return np.array([[complex(sp.simplify(tL[i,j].coeff(t,p))) for j in range(4)] for i in range(4)])
    R=coeffM(0); Ls=[coeffM(p+1) for p in range(order)]   # Ls[j] multiplies a_{n-1-j}
    # ker R (mu=0): numeric nullspace
    U,S,Vh=np.linalg.svd(R)
    ker=[Vh[i].conj() for i in range(4) if S[i]<1e-9]
    seeds=[]
    for a0 in ker:
        a=[np.array(a0,dtype=complex)]
        ok=True
        for n in range(1,order+1):
            rhs=np.zeros(4,dtype=complex)
            for j in range(min(n,len(Ls))):
                rhs+=Ls[j].dot(a[n-1-j])
            M=n*np.eye(4)-R
            try: an=np.linalg.solve(M,rhs)
            except np.linalg.LinAlgError: ok=False; break
            a.append(an)
        if ok: seeds.append(a)
    return seeds   # list of Taylor coeff lists

def eval_series(a, t):
    return sum(a[n]*(t**n) for n in range(len(a)))

def center_regular_seed(k, z0):
    Y=C.center_seed(NC,A2STAR,z0); slp=bg_slopes(Y)
    L=Lmat_num(Y,slp,k); ev,Vec=np.linalg.eig(L)
    idx=int(np.argmax(ev.real)); vec=Vec[:,idx]
    for c in [1,0,2,3]:
        if abs(vec[c])>1e-9: vec=vec/vec[c]; break
    return vec

def integrate_bg_pert(Y, hp, k, x0, x1, h):
    """RK4 both background and perturbation from x0 to x1 (h signed)."""
    x=x0; n=int(round(abs((x1-x0)/h))); h=math.copysign(h,x1-x0)
    for i in range(n):
        def d(Yr,hpc):
            slp=bg_slopes(Yr); return np.array(slp), Lmat_num(Yr,slp,k).dot(hpc)
        dY1,dh1=d(Y,hp)
        dY2,dh2=d(tuple(Y[j]+0.5*h*dY1[j] for j in range(4)),hp+0.5*h*dh1)
        dY3,dh3=d(tuple(Y[j]+0.5*h*dY2[j] for j in range(4)),hp+0.5*h*dh2)
        dY4,dh4=d(tuple(Y[j]+h*dY3[j] for j in range(4)),hp+h*dh3)
        Y=tuple(Y[j]+(h/6)*(dY1[j]+2*dY2[j]+2*dY3[j]+dY4[j]) for j in range(4))
        hp=hp+(h/6)*(dh1+2*dh2+2*dh3+dh4)
        nrm=np.max(np.abs(hp))
        if nrm>1e10: hp=hp/nrm
        if not np.all(np.isfinite(hp)): return None,None
    return Y,hp

# background sonic location (from bg cache-like quick integ)
def find_sonic(z0,h):
    Y=C.center_seed(NC,A2STAR,z0); x=math.log(z0); prevD=C.Dson(Y[0],Y[3]); prevY=Y; prevx=x
    for i in range(int((1.0-x)/h)):
        Y2=C.rk4_step(Y,h); D=C.Dson(Y2[0],Y2[3])
        if prevD<0 and D>=0:
            fr=prevD/(prevD-D); xh=prevx+fr*h
            Yh=tuple(prevY[j]+fr*(Y2[j]-prevY[j]) for j in range(4)); return xh,Yh
        prevD=D; prevY=Y2; prevx=x+h; Y=Y2; x+=h
    return None,None

def residual(k, dm=0.02, z0=math.exp(-12), h=1e-4):
    k=complex(k)
    xs,Ys=find_sonic(z0,h)
    if xs is None: return None
    x_m = xs - dm    # matching point (t=x-xs=-dm)
    # LEFT: center -> x_m
    hpL0=center_regular_seed(k,z0)
    YL,hpL=integrate_bg_pert(C.center_seed(NC,A2STAR,z0),hpL0,k,math.log(z0),x_m,h)
    if hpL is None: return None
    # RIGHT: analytic seeds at sonic (t small negative), integrate inward to x_m
    seeds=build_analytic_seeds(k)
    t_seed=-dm*0.3   # start integration a bit inside from sonic on the analytic side (t<0)
    cols=[]
    for a in seeds:
        hp0=eval_series(a, t_seed)
        # background at x=xs+t_seed via its own series
        tt=t_seed
        Ns=2/S3+(-2/S3)*tt+0.5*(2/S3)*tt**2; As=1.5+3*tt+0.5*33*tt**2
        Os=0.75+4.5*tt+0.5*49.5*tt**2; Vs=1/S3+(2/S3)*tt+0.5*(10*S3/3)*tt**2
        Yr,hpr=integrate_bg_pert((Ns,As,Os,Vs),np.array(hp0,dtype=complex),k,xs+t_seed,x_m,h)
        if hpr is None: return None
        cols.append(hpr)
    # residual: det[ hpL | cols ] (4x4). Zero when hpL in span(cols).
    Mmat=np.column_stack([hpL]+cols)
    # normalize columns to avoid scale blowup
    for c in range(Mmat.shape[1]):
        nc=np.linalg.norm(Mmat[:,c])
        if nc>0: Mmat[:,c]/=nc
    return np.linalg.det(Mmat)

if __name__=="__main__":
    print("Two-point Wronskian match residual vs k (zeros=eigenvalues):")
    print(f"{'k':>7} {'Re det':>14} {'Im det':>14} {'|det|':>12}")
    rows=[]
    for k in np.arange(0.2,3.61,0.1):
        r=residual(k)
        if r is None: print(f"{k:7.2f}  fail"); continue
        rows.append((k,r)); print(f"{k:7.2f} {r.real:14.5e} {r.imag:14.5e} {abs(r):12.4e}")
    print("\nSign changes in Re(det):")
    for i in range(1,len(rows)):
        if rows[i-1][1].real*rows[i][1].real<0: print(f"   between k={rows[i-1][0]:.2f} and {rows[i][0]:.2f}")
    print("|det| local minima:")
    for i in range(1,len(rows)-1):
        if abs(rows[i][1])<abs(rows[i-1][1]) and abs(rows[i][1])<abs(rows[i+1][1]):
            print(f"   min at k={rows[i][0]:.2f} |det|={abs(rows[i][1]):.3e}")
