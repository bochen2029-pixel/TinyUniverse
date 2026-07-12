# kha_fast_struct.py — FAST: KHA-verbatim system sonic separatrix slopes (numeric) + perturbation
# operator residue (DFT Laurent) at the ingoing sonic point V=-1/sqrt3. Determines whether the KHA
# operator supports the physical mode kappa~2.81 (indicial exponents = eig of residue R).
import sympy as sp, numpy as np, math, pickle
S3 = math.sqrt(3.0)

# ---- KHA resolved fluid slopes w'(N,A,w,V), V'(N,A,w,V) lambdified (fast numeric) ----
N, A, w, V = sp.symbols('N A w V', real=True)
Ap = A*(1 - A + (2*w/(1-V**2))*(1+V**2/3)); Np = N*(-2 + A - sp.Rational(2,3)*w)
g, u = sp.symbols('g u')
f1 = (1+N*V)*g + 4*(N+V)*u/(3*(1-V**2)) - N*V*Ap/(3*A) + 4*V*Np/3 + 2*N*(1+4*w/(9*(1-V**2)))
f2 = (4*V+N+3*N*V**2)*g + 4*(1+V**2+2*N*V)*u/(1-V**2) + N*(1-V**2)*Ap/A + 4*(1+V**2)*Np + 2*N*(1+3*V**2)
sol = sp.solve([f1,f2],[g,u],dict=True)[0]
wp_e = sp.cancel(w*sol[g]); Vp_e = sp.cancel(sol[u])
f_Np = sp.lambdify((N,A,w,V), Np, 'numpy'); f_Ap = sp.lambdify((N,A,w,V), Ap, 'numpy')
f_wp = sp.lambdify((N,A,w,V), wp_e, 'numpy'); f_Vp = sp.lambdify((N,A,w,V), Vp_e, 'numpy')

# sonic point (ingoing): N0=2/S3, A0=3/2, w0=3/4, V0=-1/S3.  metric slopes A'=3, N'=-2/S3.
N0, A0, w0, V0 = 2/S3, 1.5, 0.75, -1/S3
Ap0, Np0 = 3.0, -2/S3
print("KHA sonic (ingoing V=-1/sqrt3): N0=%.5f A0=%.3f w0=%.3f V0=%.5f, A'=%.3f N'=%.5f" % (N0,A0,w0,V0,Ap0,Np0))

# separatrix slopes (w',V') via L'Hopital: along Y'=(N',A',w',V'), numerators/det -> ratios of
# directional derivatives. Build gradients numerically by finite difference of numG,numU,det.
# Use the Cramer numerators symbolically-lambdified.
c11=(1+N*V); c12=4*(N+V)/(3*(1-V**2)); c21=(4*V+N+3*N*V**2); c22=4*(1+V**2+2*N*V)/(1-V**2)
b1=-(-N*V*Ap/(3*A)+4*V*Np/3+2*N*(1+4*w/(9*(1-V**2))))
b2=-(N*(1-V**2)*Ap/A+4*(1+V**2)*Np+2*N*(1+3*V**2))
numG = sp.cancel(b1*c22 - c12*b2); numU = sp.cancel(c11*b2 - b1*c21); det_full = sp.cancel(c11*c22 - c12*c21)
f_numG=sp.lambdify((N,A,w,V),numG,'numpy'); f_numU=sp.lambdify((N,A,w,V),numU,'numpy'); f_det=sp.lambdify((N,A,w,V),det_full,'numpy')
def grad(f,pt,e=1e-6):
    out=[]
    for i in range(4):
        p2=list(pt); p2[i]+=e; pm=list(pt); pm[i]-=e
        out.append((f(*p2)-f(*pm))/(2*e))
    return np.array(out)
pt0=[N0,A0,w0,V0]
gG=grad(f_numG,pt0); gU=grad(f_numU,pt0); gD=grad(f_det,pt0)
# unknowns W1=w', U1=V'. Yp=(Np0,Ap0,W1,U1). W1*(gD.Yp)=w0*(gG.Yp)? No: g=numG/det, W1=w*g so
# W1 = w0*numG/det; L'Hopital W1 = w0 * (gG.Yp)/(gD.Yp). Similarly U1=(gU.Yp)/(gD.Yp). Two eqns:
from sympy import symbols as _s
W1,U1=sp.symbols('W1 U1')
Yp=[Np0,Ap0,W1,U1]
DG=sum(gG[i]*Yp[i] for i in range(4)); DU=sum(gU[i]*Yp[i] for i in range(4)); DD=sum(gD[i]*Yp[i] for i in range(4))
eqW=sp.expand(W1*DD - w0*DG); eqU=sp.expand(U1*DD - DU)
slopes=sp.solve([eqW,eqU],[W1,U1],dict=True)
print("\nKHA separatrix slopes (w',V') [numeric-gradient L'Hopital]:")
branches=[]
for s in slopes:
    wv=complex(s[W1]); vv=complex(s[U1])
    real = abs(wv.imag)<1e-6 and abs(vv.imag)<1e-6
    print("   w'=%.5f%+.5fj , V'=%.5f%+.5fj  (%s)" % (wv.real,wv.imag,vv.real,vv.imag,'REAL' if real else 'complex'))
    if real: branches.append((wv.real,vv.real))

# 2nd order coeffs for a branch, numerically: enforce N'',A'' from d/dx of metric slopes, and
# w'',V'' from d/dx of the fluid slope-implicit relations. Simpler: fit the branch series by
# integrating the resolved ODE a tiny step from the sonic point using a high-order local solve.
# For the residue we only need L ~ R/t; R depends on Y0,Y'(0). L may need Y'' too? R = lim t*L; if L
# has a simple pole with residue built from Y0 and Y'(0) only (the pole comes from det~t), R needs the
# 1st-order series. So build background series to 1st order (Y0 + Y'(0) t) — sufficient for R.
d=pickle.load(open('Lmat_kha.pkl','rb')); L=sp.sympify(d['Lmat'])
N0s,A0s,w0s,V0s=sp.symbols('N0 A0 w0 V0'); Nxs,Axs,wxs,Vxs=sp.symbols('N_x A_x w_x V_x'); ksym=sp.symbols('k')
# lambdify t*L entries with background series (1st order) substituted, functions of (t,k)
tsym=sp.symbols('t')
def residue_for_branch(w1,v1,order1=True):
    # 1st-order background series
    Ns=N0+Np0*tsym; As=A0+Ap0*tsym; Ws=w0+w1*tsym; Vs=V0+v1*tsym
    sub={N0s:Ns,A0s:As,w0s:Ws,V0s:Vs,Nxs:Np0,Axs:Ap0,wxs:w1,Vxs:v1}
    Lt=L.subs(sub)   # functions of t,k (affine in k)
    # residue R = lim_{t->0} t*Lt, via DFT on small circle in t (k kept symbolic-> evaluate at test k)
    R_of_k=[]
    for ktest in [2.81055255, 0.35699, 1.0, 0.0]:
        Ltk=[[sp.lambdify(tsym, (tsym*Lt[i,j]).subs(ksym,ktest), 'numpy') for j in range(4)] for i in range(4)]
        M=64; eps=1e-4; wud=np.exp(2j*np.pi*np.arange(M)/M); tt=eps*wud
        R=np.zeros((4,4),complex)
        for i in range(4):
            for j in range(4):
                gvals=np.array([complex(Ltk[i][j](tv)) for tv in tt])
                R[i,j]=np.fft.ifft(gvals)[0]   # constant term = lim t*L
        R_of_k.append((ktest,R))
    return R_of_k

print("\n=== KHA operator residue R and indicial exps (eig R) per branch ===")
for bi,(w1,v1) in enumerate(branches):
    print(f"\n branch {bi}: w'={w1:.4f} V'={v1:.4f}")
    for ktest,R in residue_for_branch(w1,v1):
        ev=np.sort_complex(np.linalg.eigvals(R))
        print(f"   k={ktest:9.5f}: eig(R)={np.round(ev,5)}")
