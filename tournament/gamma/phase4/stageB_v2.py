# stageB_v2.py — robust eigenvalue residual via the sonic Fuchsian structure.
# Indicial exponents at sonic: mu=0 (x3, analytic), mu=1-2k (x1, the non-analytic branch for Re k>1/2).
# Physical/gauge eigenvalues: k such that the center-REGULAR solution has ZERO coefficient on the
# mu=(1-2k) eigendirection of the residue matrix R at the sonic point.
#
# Residual R(k) = (left-eigenvector of R for eigenvalue 1-2k) . hp(x_m)   [x_m just before sonic],
# normalized by |hp|. Zeros = eigenvalues. Never integrate into Dson=0.
import numpy as np, math, sympy as sp, pickle, css_core as C, pert_core as P
from scipy.optimize import brentq

S3 = math.sqrt(3.0)
A2STAR = 0.250067905275
NC = 1.0

# --- symbolic residue matrix R(k) and its left null-vector for eigenvalue (1-2k) ---
_d = pickle.load(open("Lmat.pkl","rb"))
_A0,_N0,_o0,_V0 = sp.symbols('A0 N0 om0 V0'); _Nx,_Ax,_ox,_Vx = sp.symbols('N_x A_x om_x V_x'); _k=sp.symbols('k')
_L = sp.sympify(_d['Lmat']); _x = sp.symbols('x')
_Ns=2/sp.sqrt(3) + (-2/sp.sqrt(3))*_x + sp.Rational(1,2)*(2/sp.sqrt(3))*_x**2
_As=sp.Rational(3,2)+3*_x+sp.Rational(1,2)*33*_x**2
_Os=sp.Rational(3,4)+sp.Rational(9,2)*_x+sp.Rational(1,2)*sp.Rational(99,2)*_x**2
_Vs=1/sp.sqrt(3)+(2/sp.sqrt(3))*_x+sp.Rational(1,2)*(10*sp.sqrt(3)/3)*_x**2
_sub={_N0:_Ns,_A0:_As,_o0:_Os,_V0:_Vs,_Nx:sp.diff(_Ns,_x),_Ax:sp.diff(_As,_x),_ox:sp.diff(_Os,_x),_Vx:sp.diff(_Vs,_x)}
_Lx=_L.subs(_sub)
_R=sp.zeros(4,4)
for i in range(4):
    for j in range(4):
        s=sp.series(_Lx[i,j]*_x,_x,0,1).removeO()
        _R[i,j]=sp.simplify(s.subs(_x,0)) if s!=0 else sp.Integer(0)
# left eigenvector wL: wL (R - (1-2k) I) = 0
mu = 1-2*_k
M = (_R - mu*sp.eye(4)).T          # right-null of M^T = left-null of (R-mu I)
ns = M.nullspace()
assert ns, "no left-eigenvector for mu=1-2k"
wL = ns[0].T                        # 1x4
f_wL = sp.lambdify(_k, wL, 'numpy')
print("left-eigenvector wL(k) for exponent (1-2k):", sp.simplify(wL))

def bg_slopes(Y):
    N,A,om,V=Y; return (C.Nx(A,N,om,V),C.Ax(A,N,om,V),C.omx(A,N,om,V),C.Vx(A,N,om,V))

def center_regular_seed(k, z0):
    Y=C.center_seed(NC,A2STAR,z0); Nb,Ab,ob,Vb=bg_slopes(Y)
    L=P.Lmat(Y[0],Y[1],Y[2],Y[3],Nb,Ab,ob,Vb,k)
    ev,Vec=np.linalg.eig(L)
    idx=int(np.argmax(ev.real))          # lambda=+2 regular mode
    vec=Vec[:,idx]
    for c in [1,0,2,3]:
        if abs(vec[c])>1e-9: vec=vec/vec[c]; break
    return Y, vec

def rk4_joint(Y,hp,k,h):
    def deriv(Yr,hpc):
        Nb,Ab,ob,Vb=bg_slopes(Yr); dY=(Nb,Ab,ob,Vb)
        L=P.Lmat(Yr[0],Yr[1],Yr[2],Yr[3],Nb,Ab,ob,Vb,k); return np.array(dY,float),L.dot(hpc)
    dY1,dh1=deriv(Y,hp)
    dY2,dh2=deriv(tuple(Y[i]+0.5*h*dY1[i] for i in range(4)),hp+0.5*h*dh1)
    dY3,dh3=deriv(tuple(Y[i]+0.5*h*dY2[i] for i in range(4)),hp+0.5*h*dh2)
    dY4,dh4=deriv(tuple(Y[i]+h*dY3[i] for i in range(4)),hp+h*dh3)
    Yn=tuple(Y[i]+(h/6)*(dY1[i]+2*dY2[i]+2*dY3[i]+dY4[i]) for i in range(4))
    return Yn, hp+(h/6)*(dh1+2*dh2+2*dh3+dh4)

def residual(k, dstop=0.06, z0=math.exp(-11), h=1e-4):
    """Project hp onto the singular left-eigenvector at Dson=dstop before sonic. Zero at eigenvalues."""
    Y,hp=center_regular_seed(k,z0); x=math.log(z0)
    prevD=C.Dson(Y[0],Y[3])
    n=int((2.0-x)/h)
    for i in range(n):
        Y,hp=rk4_joint(Y,hp,k,h); x+=h
        if not (all(math.isfinite(v) for v in Y) and np.all(np.isfinite(hp))): return None
        nrm=np.max(np.abs(hp))
        if nrm>1e8: hp=hp/nrm
        D=C.Dson(Y[0],Y[3])
        if D<dstop and Y[3]>0.3 and prevD>=dstop:
            # reached the stop shell just before sonic
            w=np.array(f_wL(k),dtype=complex).flatten()
            proj=w.dot(hp)/(np.max(np.abs(hp))+1e-30)
            return complex(proj)
        prevD=D
    return None

if __name__=="__main__":
    print("\nProjection residual onto singular direction vs real k:")
    print(f"{'k':>7} {'|proj|':>12} {'Re proj':>12} {'Im proj':>12}")
    for k in np.arange(0.1, 5.01, 0.2):
        r=residual(complex(k))
        if r is None: print(f"{k:7.2f}   none"); continue
        print(f"{k:7.2f} {abs(r):12.4e} {r.real:12.4e} {r.imag:12.4e}")
