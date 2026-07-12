# stageB_fast.py — corrected AND fast eigenvalue extraction for fluid-CSS beta.
# Speed: the background trajectory is kappa-INDEPENDENT. Precompute it ONCE (all RK4 stage-states),
# then for each kappa build L on the whole cached trajectory in ONE vectorized lambdify call and run
# only the cheap linear perturbation RK4 in the Python loop (no SymPy in the inner loop).
#
# Physics identical to stageB_fix.py (verified diag_residue.py):
#   Fuchsian at sonic, exponents {0,0,0,1-2k}, v_na=[0,0,3sqrt3/2,1], pole-normalized left dual
#   wL_n=(2k-3)wL. Eigenvalue condition: c_na(k)= [wL_n.hp]/(x_s-x)^(1-2k) = 0.
import numpy as np, math, sympy as sp, pickle
import css_core as C, pert_core as P

S3 = math.sqrt(3.0)
A2STAR = 0.250067905275
NC = 1.0
DSON_SLOPE = 4.0/3.0

# --- vectorized L builder: lambdify EACH entry separately (avoids inhomogeneous-array assembly when
#     some entries are constants and others are arrays). Each _fL_ij(args...) -> scalar or array. ---
_d = pickle.load(open("Lmat.pkl","rb"))
_A0,_N0,_o0,_V0 = sp.symbols('A0 N0 om0 V0'); _Nx,_Ax,_ox,_Vx = sp.symbols('N_x A_x om_x V_x'); _ksym=sp.symbols('k')
_L = sp.sympify(_d['Lmat'])
_ARGS = (_N0,_A0,_o0,_V0,_Nx,_Ax,_ox,_Vx,_ksym)
_fL_ij = [[sp.lambdify(_ARGS, _L[i,j], 'numpy') for j in range(4)] for i in range(4)]

# pole-normalized left dual + non-analytic right vector (exact, from diag_residue)
wL_n_sym = sp.Matrix([[ -5, -7*sp.sqrt(3)/9, 2*sp.sqrt(3)*(2*_ksym+1)/9, 2*_ksym-3 ]])
f_wL_n = sp.lambdify(_ksym, wL_n_sym, 'numpy')

def bg_slopes(Y):
    N,A,om,V = Y
    return (C.Nx(A,N,om,V), C.Ax(A,N,om,V), C.omx(A,N,om,V), C.Vx(A,N,om,V))

def build_background(z0, h, dstop):
    """Integrate background center->just-before-sonic ONCE. Cache RK4 stage-states.
    Returns arrays over steps of the 4 stage-states' (N,A,om,V,Nx,Ax,ox,Vx), plus the stop index and
    the interp fraction & (x_s-x) at the stop shell. Stage layout per step matches classic RK4:
      s1=Y ; s2=Y+0.5h k1 ; s3=Y+0.5h k2 ; s4=Y+h k3  (k_i are BACKGROUND slopes)."""
    Y = C.center_seed(NC, A2STAR, z0); x = math.log(z0)
    S1=[];S2=[];S3=[];S4=[]; Dser=[]; Yser=[]
    prevD = C.Dson(Y[0], Y[3])
    n = int((3.0 - x)/h)
    stop = None
    for i in range(n):
        k1 = bg_slopes(Y)
        Y2 = tuple(Y[j]+0.5*h*k1[j] for j in range(4)); k2 = bg_slopes(Y2)
        Y3 = tuple(Y[j]+0.5*h*k2[j] for j in range(4)); k3 = bg_slopes(Y3)
        Y4 = tuple(Y[j]+h*k3[j] for j in range(4));      k4 = bg_slopes(Y4)
        # store stage states with their slopes (slopes at each stage)
        S1.append(Y + k1); S2.append(Y2 + k2); S3.append(Y3 + k3); S4.append(Y4 + k4)
        Ynext = tuple(Y[j]+(h/6)*(k1[j]+2*k2[j]+2*k3[j]+k4[j]) for j in range(4))
        Dnext = C.Dson(Ynext[0], Ynext[3])
        Yser.append(Y); Dser.append(prevD)
        if Dnext < 0 and -Dnext <= dstop and -prevD > dstop:
            fr = (prevD + dstop)/(prevD - Dnext)
            fr = min(max(fr,0.0),1.0)
            Yh = tuple(Y[j] + fr*(Ynext[j]-Y[j]) for j in range(4))
            Dh = C.Dson(Yh[0], Yh[3])
            xs_minus_x = -Dh/DSON_SLOPE
            stop = (i, fr, xs_minus_x, Dh)
            # include this final step then break
            Y = Ynext; prevD = Dnext
            break
        Y = Ynext; prevD = Dnext
    if stop is None: return None
    A1=np.array(S1); A2=np.array(S2); A3=np.array(S3); A4=np.array(S4)   # (nsteps,8)
    return dict(S1=A1,S2=A2,S3=A3,S4=A4,stop=stop,h=h,z0=z0)

_BGCACHE = {}
def get_bg(z0,h,dstop):
    key=(round(math.log(z0),6),h,dstop)
    if key not in _BGCACHE:
        _BGCACHE[key]=build_background(z0,h,dstop)
    return _BGCACHE[key]

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
    """Build L (4x4) at every step for a given RK4 stage array -> tensor (4,4,nsteps)."""
    N,A,om,V,Nx,Ax,ox,Vx = [stage[:,j] for j in range(8)]
    ns = len(N)
    T = np.empty((4,4,ns), dtype=complex)
    for i in range(4):
        for j in range(4):
            mij = _fL_ij[i][j](N,A,om,V,Nx,Ax,ox,Vx,k)
            T[i,j,:] = mij if np.ndim(mij)>0 else np.full(ns, mij)
    return T

def c_na(k, dstop=0.02, z0=math.exp(-12), h=5e-5):
    bg = get_bg(z0,h,dstop)
    if bg is None: return None
    i_stop = bg['stop'][0]; fr = bg['stop'][1]; dxm = bg['stop'][2]
    ns = i_stop+1
    # precompute L tensors for the 4 stages over steps [0..i_stop]
    L1 = Lstack(bg['S1'][:ns], k); L2 = Lstack(bg['S2'][:ns], k)
    L3 = Lstack(bg['S3'][:ns], k); L4 = Lstack(bg['S4'][:ns], k)
    hp = center_regular_seed(k, z0)
    h = bg['h']; logscale = 0.0
    prevhp = hp.copy()
    for i in range(ns):
        d1 = L1[:,:,i].dot(hp)
        d2 = L2[:,:,i].dot(hp+0.5*h*d1)
        d3 = L3[:,:,i].dot(hp+0.5*h*d2)
        d4 = L4[:,:,i].dot(hp+h*d3)
        prevhp = hp
        hp = hp + (h/6)*(d1+2*d2+2*d3+d4)
        nrm = np.max(np.abs(hp))
        if not np.all(np.isfinite(hp)): return None
        if nrm > 1e12:
            hp = hp/nrm; prevhp = prevhp/nrm; logscale += math.log(nrm)
    # interpolate to the stop shell within the last step
    hph = prevhp + fr*(hp-prevhp)
    w = np.array(f_wL_n(k), dtype=complex).flatten()
    proj = w.dot(hph)
    power = dxm**(1.0-2.0*k)
    cna = (proj/power)*math.exp(logscale)
    return cna, proj, dxm, math.log(np.max(np.abs(hph)))+logscale

if __name__=="__main__":
    import sys
    dstop = float(sys.argv[1]) if len(sys.argv)>1 else 0.02
    h = float(sys.argv[2]) if len(sys.argv)>2 else 5e-5
    print(f"# corrected c_na scan  dstop={dstop} h={h}  z0=e^-12")
    print(f"{'k':>7} {'log|hp|':>12} {'|c_na|':>13} {'Re c_na':>13} {'Im c_na':>13}")
    rows=[]
    for k in np.arange(0.1, 3.61, 0.1):
        out = c_na(complex(k), dstop=dstop, h=h)
        if out is None:
            print(f"{k:7.2f}   (fail)"); continue
        cna,proj,dxm,logabs = out
        rows.append((k,logabs,cna))
        print(f"{k:7.2f} {logabs:12.3f} {abs(cna):13.4e} {cna.real:13.4e} {cna.imag:13.4e}")
    print("\nSign changes in Re(c_na):")
    for i in range(1,len(rows)):
        if rows[i-1][2].real*rows[i][2].real < 0:
            print(f"   between k={rows[i-1][0]:.2f} and k={rows[i][0]:.2f}")
    print("Boundedness (log|hp|) local minima:")
    for i in range(1,len(rows)-1):
        if rows[i][1]<rows[i-1][1] and rows[i][1]<rows[i+1][1]:
            print(f"   min at k={rows[i][0]:.2f}  log|hp|={rows[i][1]:.3f}")
