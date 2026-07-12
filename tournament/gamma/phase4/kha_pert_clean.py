# kha_pert_clean.py — re-derive the KHA perturbation operator hp'=L(x;k)hp FROM SCRATCH, carefully,
# and VERIFY its Fuchsian residue at the ingoing sonic point has the textbook indicial exponents
# {0,0,0,1-2k}. If clean, dump Lmat_kha_clean.pkl (affine in k) for the eigenvalue match.
#
# KHA autonomous EOM (verbatim, gr-qc/9503007 eq.18, s-derivs kept for the perturbation):
#   metric:  A_,x/A = 1-A+(2w/(1-V^2))(1+V^2/3)        [M_A := A_,x - A*(...) = 0]
#            N_,x/N = -2+A-2w/3                          [M_N := N_,x - N*(...) = 0]
#   energy:  (w_,s+(1+NV)w_,x)/w + 4{V V_,s+(N+V)V_,x}/(3(1-V^2)) - NV A_,x/(3A)
#            + 4V N_,x/3 + 2N(1+4w/(9(1-V^2))) = 0
#   moment:  (4V w_,s+(4V+N+3NV^2)w_,x)/w + 4{(1+V^2)V_,s+(1+V^2+2NV)V_,x}/(1-V^2)
#            + N(1-V^2)A_,x/A + 4(1+V^2)N_,x + 2N(1+3V^2) = 0
# Perturb F=H(x)+eps hp(x)e^{ks}: F_,x -> H'+eps hp' e^{ks}; F_,s -> eps k hp e^{ks}. O(eps):
#   metric gives Ap'=dH[A(...)]·hp, Np'=dH[N(...)]·hp (linear in hp).
#   fluid gives 2 eqns: cX·hp' + k cS·hp + J·hp = 0 (energy,moment), cX,cS,J on the background.
#   Solve 4x4 for hp' = L hp.
import sympy as sp, numpy as np, math, pickle
S3=math.sqrt(3.0)

N,A,w,V = sp.symbols('N A w V', real=True)
Nx,Ax,wx,Vx = sp.symbols('Nx Ax wx Vx', real=True)     # x-derivs (background slopes as symbols)
Ns_,As_,ws_,Vs_ = sp.symbols('Ns As ws Vs', real=True) # s-derivs
k = sp.symbols('k')

# metric slope RHS
ArhsN = A*(1 - A + (2*w/(1-V**2))*(1+V**2/3))
NrhsN = N*(-2 + A - sp.Rational(2,3)*w)
# fluid residuals E1(energy), E2(moment) as expressions in fields & their x,s derivs (=0 on soln)
E1 = (ws_+(1+N*V)*wx)/w + 4*(V*Vs_+(N+V)*Vx)/(3*(1-V**2)) - N*V*Ax/(3*A) + 4*V*Nx/3 + 2*N*(1+4*w/(9*(1-V**2)))
E2 = (4*V*ws_+(4*V+N+3*N*V**2)*wx)/w + 4*((1+V**2)*Vs_+(1+V**2+2*N*V)*Vx)/(1-V**2) + N*(1-V**2)*Ax/A + 4*(1+V**2)*Nx + 2*N*(1+3*V**2)

# perturbation symbols
Np_,Ap_,wp_,Vp_ = sp.symbols('Np_ Ap_ wp_ Vp_')          # hp
NpP,ApP,wpP,VpP = sp.symbols('NpP ApP wpP VpP')          # hp'
eps = sp.symbols('eps')
flds=[N,A,w,V]; hp=[Np_,Ap_,wp_,Vp_]; hpP=[NpP,ApP,wpP,VpP]
# substitution F -> F+eps hp; F_x -> F_x + eps hp'; F_s -> 0 + eps k hp  (background F_s=0)
def perturb(expr):
    s={ N:N+eps*Np_, A:A+eps*Ap_, w:w+eps*wp_, V:V+eps*Vp_,
        Nx:Nx+eps*NpP, Ax:Ax+eps*ApP, wx:wx+eps*wpP, Vx:Vx+eps*VpP,
        Ns_:eps*k*Np_, As_:eps*k*Ap_, ws_:eps*k*wp_, Vs_:eps*k*Vp_ }
    return sp.diff(expr.subs(s), eps).subs(eps,0)
# metric perturbation slopes
ApP_e = sp.cancel(sum(sp.diff(ArhsN,f)*h for f,h in zip(flds,hp)))   # = Ap'
NpP_e = sp.cancel(sum(sp.diff(NrhsN,f)*h for f,h in zip(flds,hp)))   # = Np'
# linearize fluid
L1 = sp.expand(perturb(E1)); L2 = sp.expand(perturb(E2))
# substitute NpP,ApP (metric-determined)
L1 = L1.subs({NpP:NpP_e, ApP:ApP_e}); L2 = L2.subs({NpP:NpP_e, ApP:ApP_e})
# collect: a11 wpP + a12 VpP + (rest) = 0
a11=sp.cancel(L1.coeff(wpP,1)); a12=sp.cancel(L1.coeff(VpP,1))
a21=sp.cancel(L2.coeff(wpP,1)); a22=sp.cancel(L2.coeff(VpP,1))
r1=sp.cancel(-(L1.coeff(wpP,0).coeff(VpP,0)))
r2=sp.cancel(-(L2.coeff(wpP,0).coeff(VpP,0)))
detf=sp.cancel(a11*a22-a12*a21)
wpP_e=sp.cancel((r1*a22-a12*r2)/detf)
VpP_e=sp.cancel((a11*r2-r1*a21)/detf)
rows=[NpP_e,ApP_e,wpP_e,VpP_e]
Lmat=sp.Matrix(4,4, lambda i,j: sp.cancel(sp.diff(rows[i],hp[j])))
print("KHA operator (clean) assembled. Checking affine in k:", all(sp.simplify(sp.diff(Lmat[i,j],k,2))==0 for i in range(4) for j in range(4)), flush=True)
print("detf (fluid denominator) factored:", sp.factor(sp.numer(sp.together(detf))))
# sonic locus from detf:
pt_p={N:2/sp.sqrt(3),A:sp.Rational(3,2),w:sp.Rational(3,4),V:sp.Rational(1,1)/sp.sqrt(3)}
pt_m={N:2/sp.sqrt(3),A:sp.Rational(3,2),w:sp.Rational(3,4),V:-1/sp.sqrt(3)}
print("detf at V=+1/sqrt3:", sp.simplify(detf.subs(pt_p)))
print("detf at V=-1/sqrt3:", sp.simplify(detf.subs(pt_m)))

# ---- residue at the ingoing sonic point (V=-1/sqrt3) using 1st-order background series ----
# separatrix slopes numeric (reuse the resolved slopes on the background)
g,u=sp.symbols('g u')
f1b=(1+N*V)*g+4*(N+V)*u/(3*(1-V**2))-N*V*ArhsN/(3*A)+4*V*NrhsN/3+2*N*(1+4*w/(9*(1-V**2)))
f2b=(4*V+N+3*N*V**2)*g+4*(1+V**2+2*N*V)*u/(1-V**2)+N*(1-V**2)*ArhsN/A+4*(1+V**2)*NrhsN+2*N*(1+3*V**2)
solb=sp.solve([f1b,f2b],[g,u],dict=True)[0]
f_wp=sp.lambdify((N,A,w,V), sp.cancel(w*solb[g]),'numpy'); f_Vp=sp.lambdify((N,A,w,V), sp.cancel(solb[u]),'numpy')
f_numG=sp.lambdify((N,A,w,V), sp.cancel(sp.numer(sp.together(w*solb[g]))),'numpy')
# L'Hopital slopes via numeric gradient of numerators/det
c11n=(1+N*V); c12n=4*(N+V)/(3*(1-V**2)); c21n=(4*V+N+3*N*V**2); c22n=4*(1+V**2+2*N*V)/(1-V**2)
b1n=-(-N*V*ArhsN/(3*A)+4*V*NrhsN/3+2*N*(1+4*w/(9*(1-V**2))))
b2n=-(N*(1-V**2)*ArhsN/A+4*(1+V**2)*NrhsN+2*N*(1+3*V**2))
numGn=sp.cancel(b1n*c22n-c12n*b2n); numUn=sp.cancel(c11n*b2n-b1n*c21n); detn=sp.cancel(c11n*c22n-c12n*c21n)
fG=sp.lambdify((N,A,w,V),numGn,'numpy'); fU=sp.lambdify((N,A,w,V),numUn,'numpy'); fD=sp.lambdify((N,A,w,V),detn,'numpy')
def grad(f,p,e=1e-6):
    out=[]
    for i in range(4):
        a=list(p);a[i]+=e;b=list(p);b[i]-=e;out.append((f(*a)-f(*b))/(2*e))
    return np.array(out)
N0,A0,w0,V0=2/S3,1.5,0.75,-1/S3; Np0,Ap0=-2/S3,3.0
p0=[N0,A0,w0,V0]; gG=grad(fG,p0);gU=grad(fU,p0);gD=grad(fD,p0)
W1,U1=sp.symbols('W1 U1'); Yp=[Np0,Ap0,W1,U1]
DG=sum(gG[i]*Yp[i] for i in range(4));DU=sum(gU[i]*Yp[i] for i in range(4));DD=sum(gD[i]*Yp[i] for i in range(4))
slp=sp.solve([sp.expand(W1*DD-w0*DG),sp.expand(U1*DD-DU)],[W1,U1],dict=True)
branches=[]
for s in slp:
    wv=complex(s[W1]);vv=complex(s[U1])
    if abs(wv.imag)<1e-6 and abs(vv.imag)<1e-6: branches.append((wv.real,vv.real))
print("\nreal separatrix branches (w',V') at V=-1/sqrt3:",[(round(a,4),round(b,4)) for a,b in branches])

N0s,A0s,w0s,V0s=N,A,w,V; tsym=sp.symbols('t')
for bi,(w1,v1) in enumerate(branches):
    Nser=N0+Np0*tsym;Aser=A0+Ap0*tsym;Wser=w0+w1*tsym;Vser=V0+v1*tsym
    sub={N:Nser,A:Aser,w:Wser,V:Vser,Nx:Np0,Ax:Ap0,wx:w1,Vx:v1}
    Lt=Lmat.subs(sub)
    for ktest in [2.81055255,0.35699,1.0]:
        Ltk=[[sp.lambdify(tsym,(tsym*Lt[i,j]).subs(k,ktest),'numpy') for j in range(4)] for i in range(4)]
        M=64;eps2=1e-4;wud=np.exp(2j*np.pi*np.arange(M)/M);tt=eps2*wud
        R=np.zeros((4,4),complex)
        for i in range(4):
            for j in range(4):
                gv=np.array([complex(Ltk[i][j](tv)) for tv in tt]);R[i,j]=np.fft.ifft(gv)[0]
        ev=np.sort_complex(np.linalg.eigvals(R))
        print(f"  branch{bi} k={ktest:.5f}: eig(R)={np.round(ev,5)}  (want {{0,0,0,{1-2*ktest:.4f}}})")

pickle.dump({'Lmat':sp.srepr(Lmat.subs(sp.Symbol('r',positive=True),1) if False else Lmat)}, open('Lmat_kha_clean.pkl','wb'))
print("\nwrote Lmat_kha_clean.pkl")
