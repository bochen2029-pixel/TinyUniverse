# pert_rebuild.py — INDEPENDENT re-derivation of the linear perturbation operator L(x;k), fully
# consistent with the verified background (css_symbols machinery). Then dump Lmat2.pkl and print the
# residue matrix + sonic indicial exponents for cross-check against diag_residue.py.
#
# Ansatz: F(s,x) = H(x) + eps * hp(x) * e^{k s}, F in {N,A,om,V}.  s=-ln(-t), x=ln(-r/t).
# Total time derivative on a field: d/dt f = (-1/t)(f_s + f_x).  For the perturbation piece
# f = eps hp(x) e^{k s}: f_s = k f, f_x = eps hp'(x) e^{k s}. Background piece: f_s=0.
# We build nabla_a T^{ab}=0 and the metric constraints WITH s-derivatives, linearize, divide e^{ks}.
import sympy as sp, pickle

tau, r, th = sp.symbols('tau r theta', positive=True)
eps, k = sp.symbols('epsilon k')
xx = sp.symbols('x')   # not used directly; fields are symbols

# Background field symbols + their x-derivatives (H, H')
N0,A0,o0,V0 = sp.symbols('N0 A0 om0 V0')
Nx,Ax,ox,Vx = sp.symbols('N_x A_x om_x V_x')
# perturbation amplitudes + their x-derivatives
Np,Ap,op,Vp = sp.symbols('Np Ap op Vp')
Npx,Apx,opx,Vpx = sp.symbols('Npx Apx opx Vpx')

# A field and its total value F = H + eps hp e^{k s}. We encode each field as a small object carrying
# (background_value, background_xderiv, pert_amp, pert_xderiv). Represent the full field via a symbol
# and provide its x-derivative and s-derivative explicitly through substitution dictionaries.
# Simplest robust approach: introduce full-field symbols Ff and their x/s derivatives, build E(Ff, Ff_x, Ff_s),
# then substitute the linearization.
NF,AF,oF,VF = sp.symbols('NF AF oF VF')                 # full fields
NFx,AFx,oFx,VFx = sp.symbols('NFx AFx oFx VFx')          # full x-derivs
NFs,AFs,oFs,VFs = sp.symbols('NFs AFs oFs VFs')          # full s-derivs
flds_f  = [NF,AF,oF,VF]; flds_fx=[NFx,AFx,oFx,VFx]; flds_fs=[NFs,AFs,oFs,VFs]

E = r/tau
a = sp.sqrt(AF); alpha = NF*a*E; W = 1/sp.sqrt(1-VF**2)
rho = oF/(4*sp.pi*r**2*AF); p = rho/3
g = sp.diag(-alpha**2, a**2, r**2, r**2*sp.sin(th)**2)
gi = g.inv()
coords=[tau,r,th,sp.symbols('phi')]

def Dk(expr, kk):
    """total derivative wrt coordinate kk, fields depend on (s,x): d/dt=(-1/tau)(d_s+d_x), d/dr=(1/r)d_x."""
    e = sp.diff(expr, coords[kk])
    if kk==0:   # time: bring in BOTH s and x field-derivs
        for F,Fx,Fs in zip(flds_f,flds_fx,flds_fs):
            e += sp.diff(expr,F)*(Fx+Fs)*(-1/tau)
    elif kk==1: # radius: x field-deriv only
        for F,Fx in zip(flds_f,flds_fx):
            e += sp.diff(expr,F)*Fx*(1/r)
    return e

def christ(l,i,j):
    ss=sp.S.Zero
    for dd in range(4):
        if gi[l,dd]==0: continue
        ss+=gi[l,dd]*(Dk(g[dd,i],j)+Dk(g[dd,j],i)-Dk(g[i,j],dd))
    return ss/2
uu=[W/alpha, W*VF/a, sp.S.Zero, sp.S.Zero]
T=[[(rho+p)*uu[i]*uu[j]+p*gi[i,j] for j in range(4)] for i in range(4)]
def divT(b):
    ss=sp.S.Zero
    for ai in range(4):
        ss+=Dk(T[ai][b],ai)
        for c in range(4):
            g1=christ(ai,ai,c); ss+= g1*T[c][b] if g1!=0 else 0
            g2=christ(b,ai,c);  ss+= g2*T[ai][c] if g2!=0 else 0
    return ss

print("building full-field divergence (energy, momentum) with s-derivs...", flush=True)
Dtau=divT(0); Dr=divT(1)
def fin(D):
    D=sp.expand(D*4*sp.pi*r**2); D=D.subs(tau,r); return D   # E=r/tau=1 on tau=r slice; overall E scaling drops
en=fin(Dtau); mo=fin(Dr)
# metric constraints (polar-radial), full-field form (no s-deriv): M_A = A_x - A(1-A+2om/(1-V^2)(1+V^2/3)); M_N = N_x - N(-2+A-2om/3)
MA = AFx - AF*(1-AF+(2*oF/(1-VF**2))*(1+VF**2/3))
MN = NFx - NF*(-2+AF-sp.Rational(2,3)*oF)

# ---- linearize: F = H + eps hp e^{ks}. Substitute:
#   F -> H + eps*hp ;  F_x -> H' + eps*hp' ;  F_s -> 0 + eps*k*hp  (background s-deriv=0)
# then take d/d eps at eps=0.
subs_lin = {
 NF:N0+eps*Np, AF:A0+eps*Ap, oF:o0+eps*op, VF:V0+eps*Vp,
 NFx:Nx+eps*Npx, AFx:Ax+eps*Apx, oFx:ox+eps*opx, VFx:Vx+eps*Vpx,
 NFs:eps*k*Np, AFs:eps*k*Ap, oFs:eps*k*op, VFs:eps*k*Vp,
}
def linearize(expr):
    e=expr.subs(subs_lin)
    return sp.diff(e, eps).subs(eps,0)
print("linearizing...", flush=True)
Len = sp.expand(linearize(en)); Lmo=sp.expand(linearize(mo))
LMA = sp.expand(linearize(MA)); LMN=sp.expand(linearize(MN))
# metric perturbation slopes: LMA=0 -> Apx = ... ; LMN=0 -> Npx = ...
Apx_sol = sp.solve(LMA, Apx)[0]
Npx_sol = sp.solve(LMN, Npx)[0]
# substitute into fluid perturbations, then solve the 2x2 for (opx,Vpx)
Len2 = Len.subs({Apx:Apx_sol, Npx:Npx_sol})
Lmo2 = Lmo.subs({Apx:Apx_sol, Npx:Npx_sol})
a11=Len2.coeff(opx,1); a12=Len2.coeff(Vpx,1)
a21=Lmo2.coeff(opx,1); a22=Lmo2.coeff(Vpx,1)
b1=-(Len2.coeff(opx,0).coeff(Vpx,0)); b2=-(Lmo2.coeff(opx,0).coeff(Vpx,0))
detf=sp.cancel(a11*a22-a12*a21)
opx_sol=sp.cancel((b1*a22-a12*b2)/detf)
Vpx_sol=sp.cancel((a11*b2-b1*a21)/detf)
# assemble L: hp'=(Npx,Apx,opx,Vpx) each linear in hp=(Np,Ap,op,Vp)
rows=[sp.cancel(Npx_sol), sp.cancel(Apx_sol), opx_sol, Vpx_sol]
hp=[Np,Ap,op,Vp]
Lmat=sp.Matrix(4,4, lambda i,j: sp.cancel(sp.diff(rows[i],hp[j])))
Lmat=Lmat.subs(r,1)
pickle.dump({'Lmat':sp.srepr(Lmat)}, open("Lmat2.pkl","wb"))
print("wrote Lmat2.pkl", flush=True)

# ---- cross-check: residue matrix at the exact sonic point (branch-1 series) ----
xs=sp.symbols('xs')
Ns=2/sp.sqrt(3)+(-2/sp.sqrt(3))*xs+sp.Rational(1,2)*(2/sp.sqrt(3))*xs**2
As=sp.Rational(3,2)+3*xs+sp.Rational(1,2)*33*xs**2
Os=sp.Rational(3,4)+sp.Rational(9,2)*xs+sp.Rational(1,2)*sp.Rational(99,2)*xs**2
Vs=1/sp.sqrt(3)+(2/sp.sqrt(3))*xs+sp.Rational(1,2)*(10*sp.sqrt(3)/3)*xs**2
sub={N0:Ns,A0:As,o0:Os,V0:Vs,Nx:sp.diff(Ns,xs),Ax:sp.diff(As,xs),ox:sp.diff(Os,xs),Vx:sp.diff(Vs,xs)}
Lx=Lmat.subs(sub)
R=sp.zeros(4,4)
for i in range(4):
    for j in range(4):
        ser=sp.series(Lx[i,j]*xs,xs,0,1).removeO()
        R[i,j]=sp.nsimplify(sp.simplify(ser.subs(xs,0))) if ser!=0 else sp.Integer(0)
print("\n=== REBUILT residue matrix R ===")
sp.pprint(sp.simplify(R))
lam=sp.symbols('lam')
print("charpoly(lam) =", sp.factor(sp.simplify(R.charpoly(lam).as_expr())))
print("eig(R):", {sp.simplify(e):m for e,m in R.eigenvals().items()})
