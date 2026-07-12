# bg_analyze.py — load cached background; integrate the linear perturbation for a given kappa along
# it; extract the non-analytic amplitude c_na (eigenvalue residual) with the correct normalization,
# and a boundedness profile. Fast: no SymPy in the inner loop (L built vectorized over steps once/k).
import numpy as np, math, sympy as sp, pickle
import css_core as C, pert_core as P

S3 = math.sqrt(3.0)
A2STAR = 0.250067905275
NC = 1.0
DSON_SLOPE = 4.0/3.0        # Dson ~ (4/3)(x-x_s) near sonic (verified diag_residue)

_d = pickle.load(open("Lmat.pkl","rb"))
_A0,_N0,_o0,_V0 = sp.symbols('A0 N0 om0 V0'); _Nx,_Ax,_ox,_Vx = sp.symbols('N_x A_x om_x V_x'); _ksym=sp.symbols('k')
_L = sp.sympify(_d['Lmat'])
_ARGS = (_N0,_A0,_o0,_V0,_Nx,_Ax,_ox,_Vx,_ksym)
_fL_ij = [[sp.lambdify(_ARGS, _L[i,j], 'numpy') for j in range(4)] for i in range(4)]

# pole-normalized left dual (from diag_residue). c_na normalization: wL_n . v_na = 4k-2, divide it out.
wL_n_sym = sp.Matrix([[ -5, -7*sp.sqrt(3)/9, 2*sp.sqrt(3)*(2*_ksym+1)/9, 2*_ksym-3 ]])
f_wL_n = sp.lambdify(_ksym, wL_n_sym, 'numpy')
v_na = np.array([0.0,0.0,3*S3/2.0,1.0], dtype=complex)

def bg_slopes(Y):
    N,A,om,V = Y
    return (C.Nx(A,N,om,V), C.Ax(A,N,om,V), C.omx(A,N,om,V), C.Vx(A,N,om,V))

def center_regular_seed(k, z0):
    Y = C.center_seed(NC, A2STAR, z0)
    Nb,Ab,ob,Vb = bg_slopes(Y)
    L = P.Lmat(Y[0],Y[1],Y[2],Y[3], Nb,Ab,ob,Vb, k)
    ev, Vec = np.linalg.eig(L)
    idx = int(np.argmax(ev.real))
    vec = Vec[:, idx]
    for c in [1,0,2,3]:
        if abs(vec[c]) > 1e-9:
            vec = vec/vec[c]; break
    return vec

def Lstack(stage, k):
    N,A,om,V,Nx,Ax,ox,Vx = [stage[:,j] for j in range(8)]
    ns=len(N); T=np.empty((4,4,ns),dtype=complex)
    for i in range(4):
        for j in range(4):
            mij=_fL_ij[i][j](N,A,om,V,Nx,Ax,ox,Vx,k)
            T[i,j,:]= mij if np.ndim(mij)>0 else np.full(ns,mij)
    return T

_LC={}
def Lcache(bg,k):
    key=id(bg),k
    if key not in _LC:
        ns=bg['nsteps']
        _LC[key]=(Lstack(bg['S1_full'][:ns],k),Lstack(bg['S2'],k),Lstack(bg['S3'],k),Lstack(bg['S4'],k))
    return _LC[key]

def integrate(bg, k, dstop, record_profile=False):
    """Integrate perturbation along cached bg; stop-shell at |Dson|=dstop. Return c_na and profile."""
    L1,L2,L3,L4 = Lcache(bg,k)
    h=bg['h']; ns=bg['nsteps']
    hp=center_regular_seed(k, bg['z0']); logscale=0.0
    Dn=bg['Dnode']
    prof=[]  # (xs_minus_x, log|hp|+logscale)
    i_stop=None; fr=0.0
    prevhp=hp.copy()
    for i in range(ns):
        d1=L1[:,:,i].dot(hp)
        d2=L2[:,:,i].dot(hp+0.5*h*d1)
        d3=L3[:,:,i].dot(hp+0.5*h*d2)
        d4=L4[:,:,i].dot(hp+h*d3)
        prevhp=hp
        hp=hp+(h/6)*(d1+2*d2+2*d3+d4)
        if not np.all(np.isfinite(hp)): return None
        nrm=np.max(np.abs(hp))
        if nrm>1e12: hp=hp/nrm; prevhp=prevhp/nrm; logscale+=math.log(nrm)
        Dcur=Dn[i+1] if i+1<len(Dn) else Dn[-1]
        if record_profile:
            dxm=-Dcur/DSON_SLOPE
            if dxm>0: prof.append((dxm, math.log(nrm)+logscale))
        # stop shell
        if i_stop is None and Dcur<0 and -Dcur<=dstop:
            Dprev=Dn[i] if i<len(Dn) else Dn[-1]
            if -Dprev>dstop:
                fr=(Dprev+dstop)/(Dprev-Dcur); fr=min(max(fr,0.0),1.0)
                hph=prevhp+fr*(hp-prevhp)
                Dh=Dprev+fr*(Dcur-Dprev)
                dxm=-Dh/DSON_SLOPE
                w=np.array(f_wL_n(k),dtype=complex).flatten()
                proj=w.dot(hph)
                power=dxm**(1.0-2.0*k)
                norm=(4.0*k-2.0)         # wL_n . v_na
                cna=(proj/power/norm)*math.exp(logscale)
                i_stop=i
                res=dict(cna=cna, proj=proj, dxm=dxm, logabs=math.log(np.max(np.abs(hph)))+logscale)
                if not record_profile:
                    return res
                res_saved=res
    if record_profile:
        return dict(profile=np.array(prof), stop=(res_saved if i_stop is not None else None))
    return None

def c_na(bg, k, dstop):
    r=integrate(bg,k,dstop,record_profile=False)
    return None if r is None else (r['cna'], r['proj'], r['dxm'], r['logabs'])

if __name__=="__main__":
    import sys
    bgfn=sys.argv[1] if len(sys.argv)>1 else "bgcache.pkl"
    dstop=float(sys.argv[2]) if len(sys.argv)>2 else 0.01
    bg=pickle.load(open(bgfn,"rb"))
    print(f"# loaded {bgfn}: nsteps={bg['nsteps']} h={bg['h']} dmin={bg['dmin']} finalDson={bg['Dnode'][-1]:.2e}")
    print(f"# c_na scan (normalized), dstop={dstop}")
    print(f"{'k':>7} {'Re c_na':>14} {'Im c_na':>14} {'|c_na|':>13} {'log|hp|':>11}")
    rows=[]
    for k in np.arange(0.1,3.61,0.1):
        out=c_na(bg,complex(k),dstop)
        if out is None: print(f"{k:7.2f} (fail)"); continue
        cna,proj,dxm,logabs=out
        rows.append((k,cna,logabs)); print(f"{k:7.2f} {cna.real:14.5e} {cna.imag:14.5e} {abs(cna):13.4e} {logabs:11.3f}")
    print("\nSign changes in Re(c_na):")
    for i in range(1,len(rows)):
        if rows[i-1][1].real*rows[i][1].real<0:
            print(f"   between k={rows[i-1][0]:.2f} and k={rows[i][0]:.2f}")
