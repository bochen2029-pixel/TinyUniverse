# verify_divT.py — independently verify the fluid equations by computing nabla_mu T^{mu nu}=0 with an
# EXPLICIT, from-scratch covariant divergence  nabla_mu T^{mu b} = (1/sqrt(-g)) d_mu(sqrt(-g) T^{mu b})
#   + Gamma^b_{mu c} T^{mu c}, and compare to css_symbols' en, mo (which used the other Christoffel form).
# If they agree, the css fluid system is trustworthy (independent of KHA-printed).
import sympy as sp, pickle
tau,r,th,ph=sp.symbols('tau r theta phi', positive=True)
A0,N0,o0,V0=sp.symbols('A0 N0 om0 V0'); Ax,Nx,ox,Vx=sp.symbols('A_x N_x om_x V_x')
E=r/tau; a=sp.sqrt(A0); alpha=N0*a*E; W=1/sp.sqrt(1-V0**2)
rho=o0/(4*sp.pi*r**2*A0); p=rho/3
coords=[tau,r,th,ph]
g=sp.diag(-alpha**2,a**2,r**2,r**2*sp.sin(th)**2); gi=g.inv()
detg=sp.simplify(g.det()); sqrtmg=sp.sqrt(-detg)

def Dk(expr,kk):
    e=sp.diff(expr,coords[kk])
    fac=(-1/tau) if kk==0 else ((1/r) if kk==1 else 0)
    e+=(sp.diff(expr,A0)*Ax+sp.diff(expr,N0)*Nx+sp.diff(expr,o0)*ox+sp.diff(expr,V0)*Vx)*fac
    return e
# u^a upper
uu=[W/alpha, W*V0/a, 0, 0]
Tuu=[[(rho+p)*uu[i]*uu[j]+p*gi[i,j] for j in range(4)] for i in range(4)]

def christ(l,i,j):
    s=sp.S.Zero
    for dd in range(4):
        if gi[l,dd]==0: continue
        s+=gi[l,dd]*(Dk(g[dd,i],j)+Dk(g[dd,j],i)-Dk(g[i,j],dd))
    return s/2

# form 1: (1/sqrt-g) d_mu(sqrt-g T^{mu b}) + Gamma^b_{mu c} T^{mu c}
def divT_form1(b):
    s=sp.S.Zero
    for mu in range(4):
        s+=Dk(sqrtmg*Tuu[mu][b],mu)
    s=s/sqrtmg
    for mu in range(4):
        for c in range(4):
            G=christ(b,mu,c)
            if G!=0: s+=G*Tuu[mu][c]
    return s

Es=sp.symbols('E',positive=True)
def fin(D):
    D=sp.simplify(D*4*sp.pi*r**2); D=D.subs(tau,r/Es); return sp.simplify(D.subs(Es,1))

print("computing div(T) form1 (energy, momentum) with sqrt(-g) formula...", flush=True)
en1=fin(divT_form1(0)); mo1=fin(divT_form1(1))

# load css en/mo? css_symbols didn't save en/mo directly; recompute css form (form2) here to compare.
def divT_form2(b):
    ss=sp.S.Zero
    for ai in range(4):
        ss+=Dk(Tuu[ai][b],ai)
        for c in range(4):
            g1=christ(ai,ai,c); ss+= g1*Tuu[c][b] if g1!=0 else 0
            g2=christ(b,ai,c);  ss+= g2*Tuu[ai][c] if g2!=0 else 0
    return ss
en2=fin(divT_form2(0)); mo2=fin(divT_form2(1))
print("energy form1==form2 ?", sp.simplify(en1-en2)==0)
print("moment form1==form2 ?", sp.simplify(mo1-mo2)==0)

# Now reduce to resolved (om',V') and compare to css_syms.pkl
Axv=A0*(1-A0+(2*o0/(1-V0**2))*(1+V0**2/3)); Nxv=N0*(-2+A0-sp.Rational(2,3)*o0)
e1=sp.expand(en1.subs({Ax:Axv,Nx:Nxv})); e2=sp.expand(mo1.subs({Ax:Axv,Nx:Nxv}))
a11=e1.coeff(ox);a12=e1.coeff(Vx);a21=e2.coeff(ox);a22=e2.coeff(Vx)
b1=-(e1.subs({ox:0,Vx:0}));b2=-(e2.subs({ox:0,Vx:0}))
det=a11*a22-a12*a21; Vxx=sp.simplify((a11*b2-b1*a21)/det); omx=sp.simplify((b1*a22-a12*b2)/det)
d=pickle.load(open("css_syms.pkl","rb")); rr=sp.Symbol('r',positive=True)
print("\ncss Vxx match ?", sp.simplify(Vxx.subs(r,1)-sp.sympify(d['Vxx']).subs(rr,1))==0)
print("css omx match ?", sp.simplify(omx.subs(r,1)-sp.sympify(d['omx']).subs(rr,1))==0)
# sonic locus from this det
detnum=sp.factor(sp.numer(sp.together(det.subs(r,1))))
print("num(det) factored =", detnum)
