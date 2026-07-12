# stageB_match_fast.py — fast two-point Wronskian match using a cached background.
# LEFT: center-regular solution integrated outward (cached bg L-stacks).
# RIGHT: 3 sonic-analytic Taylor solutions integrated INWARD to the matching node.
# residual = det[ hpL | hpR1 | hpR2 | hpR3 ] (column-normalized). Zeros = eigenvalues. No divergent branch.
import numpy as np, math, sympy as sp, pickle, sys
import css_core as C, pert_core as P

S3=math.sqrt(3.0); A2STAR=0.250067905275; NC=1.0
_d=pickle.load(open("Lmat.pkl","rb"))
A0,N0,o0,V0=sp.symbols('A0 N0 om0 V0'); Nxs,Axs,oxs,Vxs=sp.symbols('N_x A_x om_x V_x'); ksym=sp.symbols('k')
Lsym=sp.sympify(_d['Lmat']); _ARGS=(N0,A0,o0,V0,Nxs,Axs,oxs,Vxs,ksym)
fL_ij=[[sp.lambdify(_ARGS,Lsym[i,j],'numpy') for j in range(4)] for i in range(4)]

def Lstack(stage,k):
    N,A,om,V,Nx,Ax,ox,Vx=[stage[:,j] for j in range(8)]; ns=len(N); T=np.empty((4,4,ns),complex)
    for i in range(4):
        for j in range(4):
            m=fL_ij[i][j](N,A,om,V,Nx,Ax,ox,Vx,k); T[i,j,:]=m if np.ndim(m)>0 else np.full(ns,m)
    return T
_LC={}
def Lcache(bg,k):
    key=id(bg),k
    if key not in _LC:
        ns=bg['nsteps']; _LC[key]=(Lstack(bg['S1_full'][:ns],k),Lstack(bg['S2'],k),Lstack(bg['S3'],k),Lstack(bg['S4'],k))
    return _LC[key]

def bg_slopes(Y):
    N,A,om,V=Y; return (C.Nx(A,N,om,V),C.Ax(A,N,om,V),C.omx(A,N,om,V),C.Vx(A,N,om,V))

def center_regular_seed(k):
    Y=C.center_seed(NC,A2STAR,math.exp(-12)); slp=bg_slopes(Y)
    L=P.Lmat(Y[0],Y[1],Y[2],Y[3],*slp,k); ev,Vec=np.linalg.eig(L); idx=int(np.argmax(ev.real)); v=Vec[:,idx]
    for c in [1,0,2,3]:
        if abs(v[c])>1e-9: v=v/v[c]; break
    return v

# --- Precompute the Laurent coefficient matrices M_p(k) of tL = R + M1 t + M2 t^2 + ... ONCE
#     (symbolic in k), lambdified. Per-k call is then pure numpy. ---
_ORDER=6
def _precompute_ML():
    t=sp.symbols('t')
    Ns=2/sp.sqrt(3)+(-2/sp.sqrt(3))*t+sp.Rational(1,2)*(2/sp.sqrt(3))*t**2
    As=sp.Rational(3,2)+3*t+sp.Rational(1,2)*33*t**2
    Os=sp.Rational(3,4)+sp.Rational(9,2)*t+sp.Rational(1,2)*sp.Rational(99,2)*t**2
    Vs=1/sp.sqrt(3)+(2/sp.sqrt(3))*t+sp.Rational(1,2)*(10*sp.sqrt(3)/3)*t**2
    sub={N0:Ns,A0:As,o0:Os,V0:Vs,Nxs:sp.diff(Ns,t),Axs:sp.diff(As,t),oxs:sp.diff(Os,t),Vxs:sp.diff(Vs,t)}
    Lt=Lsym.subs(sub)   # still symbolic in k
    Mfuncs=[]  # Mfuncs[p][i][j] = lambda k -> coeff of t^p in (t*Lt)[i,j]
    for p in range(_ORDER+1):
        row=[]
        for i in range(4):
            r2=[]
            for j in range(4):
                c=sp.series(Lt[i,j]*t,t,0,_ORDER+2).removeO().coeff(t,p)
                r2.append(sp.lambdify(ksym, c, 'numpy'))
            row.append(r2)
        Mfuncs.append(row)
    return Mfuncs
_MFUNCS=_precompute_ML()
def _Mp(p,k):
    return np.array([[complex(_MFUNCS[p][i][j](k)) for j in range(4)] for i in range(4)])
def build_analytic_seeds(k, order=_ORDER):
    R=_Mp(0,k); Ls=[_Mp(p+1,k) for p in range(order)]
    U,S,Vh=np.linalg.svd(R); ker=[Vh[i].conj() for i in range(4) if S[i]<1e-9]
    seeds=[]
    for a0 in ker:
        a=[np.array(a0,complex)]; ok=True
        for n in range(1,order+1):
            rhs=np.zeros(4,complex)
            for j in range(min(n,len(Ls))): rhs+=Ls[j].dot(a[n-1-j])
            try: a.append(np.linalg.solve(n*np.eye(4)-R,rhs))
            except np.linalg.LinAlgError: ok=False; break
        if ok: seeds.append(a)
    return R,seeds
def eval_series(a,t): return sum(a[n]*(t**n) for n in range(len(a)))

def rk4_pert_stack(L1,L2,L3,L4,hp,h,i0,i1):
    """integrate perturbation over cached bg steps [i0,i1) (i1>i0 forward; i1<i0 backward)."""
    hp=hp.copy(); step=1 if i1>i0 else -1; hh=h*step
    idx=range(i0,i1,step)
    for i in idx:
        # for backward, use stage matrices at step i (approx; stages symmetric enough for small h)
        d1=L1[:,:,i].dot(hp); d2=L2[:,:,i].dot(hp+0.5*hh*d1); d3=L3[:,:,i].dot(hp+0.5*hh*d2); d4=L4[:,:,i].dot(hp+hh*d3)
        hp=hp+(hh/6)*(d1+2*d2+2*d3+d4)
        nrm=np.max(np.abs(hp))
        if nrm>1e10: hp=hp/nrm
        if not np.all(np.isfinite(hp)): return None
    return hp

def residual(bg,k,dm=0.02):
    k=complex(k); L1,L2,L3,L4=Lcache(bg,k); Dn=bg['Dnode']; ns=bg['nsteps']; h=bg['h']
    # matching node: where Dson ~ -(4/3)dm  (t=-dm). find node index
    Dtar=-(4.0/3.0)*dm
    im=None
    for i in range(ns):
        if Dn[i]<=Dtar: im=i; break
    if im is None: im=ns-1
    # LEFT: center (i=0) -> im
    hpL=center_regular_seed(k)
    hpL=rk4_pert_stack(L1,L2,L3,L4,hpL,h,0,im)
    if hpL is None: return None
    # RIGHT: analytic seeds at the last cached node (closest to sonic, i=ns-1), integrate inward to im.
    R,seeds=build_analytic_seeds(k)
    if len(seeds)<3: return None
    # seed at node ns-1: t = Dnode[ns-1]/(4/3)
    t_end=Dn[ns-1]/(4.0/3.0)
    cols=[]
    for a in seeds[:3]:
        hp0=np.array(eval_series(a,t_end),complex)
        hpR=rk4_pert_stack(L1,L2,L3,L4,hp0,h,ns-1,im)
        if hpR is None: return None
        cols.append(hpR)
    M=np.column_stack([hpL]+cols)
    for c in range(4):
        n=np.linalg.norm(M[:,c]);
        if n>0: M[:,c]/=n
    return np.linalg.det(M)

if __name__=="__main__":
    bg=pickle.load(open(sys.argv[1] if len(sys.argv)>1 else "bg_h1e4.pkl","rb"))
    dm=float(sys.argv[2]) if len(sys.argv)>2 else 0.03
    print(f"# fast Wronskian match  bg nsteps={bg['nsteps']}  dm={dm}")
    print(f"{'k':>7} {'Re det':>14} {'Im det':>14} {'|det|':>12}")
    rows=[]
    for k in np.arange(0.2,3.61,0.1):
        r=residual(bg,k,dm)
        if r is None: print(f"{k:7.2f} fail"); continue
        rows.append((k,r)); print(f"{k:7.2f} {r.real:14.5e} {r.imag:14.5e} {abs(r):12.4e}")
    print("\nSign changes Re(det):")
    for i in range(1,len(rows)):
        if rows[i-1][1].real*rows[i][1].real<0: print(f"   {rows[i-1][0]:.2f}..{rows[i][0]:.2f}")
    print("|det| minima:")
    for i in range(1,len(rows)-1):
        if abs(rows[i][1])<abs(rows[i-1][1]) and abs(rows[i][1])<abs(rows[i+1][1]): print(f"   k={rows[i][0]:.2f} |det|={abs(rows[i][1]):.3e}")
