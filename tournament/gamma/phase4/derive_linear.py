# derive_linear.py — derive the fluid equations WITH the s-derivative terms, then linearize about
# the background critical solution for the perturbation ansatz h(s,x)=H(x)+eps h_p(x) e^{k s}.
#
# Strategy: the reduced fields depend on (s,x). The autonomous reduction set d/ds=0. To restore the
# s-derivatives we note the SAME nabla_a T^{ab}=0, but now with fields functions of (s,x) where
# s=-ln(-t), x=ln(-r/t). Then d/dt at fixed r pulls in BOTH ds and dx pieces:
#   ds/dt = -1/t ,  dx/dt = -1/t   (since s=-ln(-t): ds/dt=-1/t; x=ln(-r/t): dx/dt = -1/t).
#   dr at fixed t: dx/dr = 1/r ; ds/dr=0.
# So d/dtau f = f_s ds/dt + f_x dx/dt = (-1/t)(f_s + f_x) ;  d/dr f = (1/r) f_x.
# This is the ONLY change vs css_symbols.py (there we had d/dtau f -> (-1/tau) f_x only, i.e. f_s=0).
#
# We keep s-derivatives as new symbols A_s,N_s,om_s,V_s and produce the fluid pair with them, then
# linearize.  We reuse the verified full-4D T^{ab} construction.
import sympy as sp, pickle
tau, r, th = sp.symbols('tau r theta', positive=True)
A0,N0,o0,V0 = sp.symbols('A0 N0 om0 V0')
Ax,Nx,ox,Vx = sp.symbols('A_x N_x om_x V_x')
As,Ns,os,Vs = sp.symbols('A_s N_s om_s V_s')

E = r/tau
a = sp.sqrt(A0); alpha = N0*a*E; W = 1/sp.sqrt(1-V0**2)
rho = o0/(4*sp.pi*r**2*A0); p = rho/3
g = sp.diag(-alpha**2, a**2, r**2, r**2*sp.sin(th)**2)
gi = g.inv()
coords=[tau,r,th,sp.symbols('phi')]

def Dk(expr,k):
    e = sp.diff(expr, coords[k])
    if k==0:      # d/dtau: (-1/tau)(f_x + f_s)
        e += (sp.diff(expr,A0)*(Ax+As)+sp.diff(expr,N0)*(Nx+Ns)
              +sp.diff(expr,o0)*(ox+os)+sp.diff(expr,V0)*(Vx+Vs))*(-1/tau)
    elif k==1:    # d/dr: (1/r) f_x
        e += (sp.diff(expr,A0)*Ax+sp.diff(expr,N0)*Nx
              +sp.diff(expr,o0)*ox+sp.diff(expr,V0)*Vx)*(1/r)
    return e

def christ(l,i,j):
    ss=sp.S.Zero
    for dd in range(4):
        if gi[l,dd]==0: continue
        ss+=gi[l,dd]*(Dk(g[dd,i],j)+Dk(g[dd,j],i)-Dk(g[i,j],dd))
    return sp.simplify(ss/2)
uu=[W/alpha, W*V0/a, sp.S.Zero, sp.S.Zero]
T=[[(rho+p)*uu[i]*uu[j]+p*gi[i,j] for j in range(4)] for i in range(4)]
def divT(b):
    ss=sp.S.Zero
    for ai in range(4):
        ss+=Dk(T[ai][b],ai)
        for c in range(4):
            g1=christ(ai,ai,c)
            if g1!=0: ss+=g1*T[c][b]
            g2=christ(b,ai,c)
            if g2!=0: ss+=g2*T[ai][c]
    return ss
print("energy (b=tau) with s-derivs...", flush=True)
Dtau=divT(0)
print("momentum (b=r) with s-derivs...", flush=True)
Dr=divT(1)
Es=sp.symbols('E',positive=True)
def fin(D):
    D=sp.simplify(D*4*sp.pi*r**2); D=D.subs(tau,r/Es); return sp.simplify(D.subs(Es,1))
en=fin(Dtau); mo=fin(Dr)
# metric constraints unchanged (no s-deriv): reuse
Axv=A0*(1-A0+(2*o0/(1-V0**2))*(1+V0**2/3)); Nxv=N0*(-2+A0-sp.Rational(2,3)*o0)
pickle.dump({'en':sp.srepr(en),'mo':sp.srepr(mo),'Axv':sp.srepr(Axv),'Nxv':sp.srepr(Nxv)},
            open("fluid_sderiv.pkl","wb"))
# sanity: with s-derivs=0 must reproduce css_syms fluid pair (same numV etc.)
en0=en.subs({As:0,Ns:0,os:0,Vs:0}); mo0=mo.subs({As:0,Ns:0,os:0,Vs:0})
d=pickle.load(open("css_syms.pkl","rb")); rr=sp.Symbol('r',positive=True)
# rebuild css a11..: from css we only saved coeffs; just check the resolved slopes match
e1=sp.expand(en0.subs({Ax:Axv,Nx:Nxv})); e2=sp.expand(mo0.subs({Ax:Axv,Nx:Nxv}))
a11=e1.coeff(ox);a12=e1.coeff(Vx);a21=e2.coeff(ox);a22=e2.coeff(Vx)
b1=-(e1.subs({ox:0,Vx:0}));b2=-(e2.subs({ox:0,Vx:0}))
det=sp.simplify(a11*a22-a12*a21); numV=sp.simplify(a11*b2-b1*a21)
Vxx=sp.simplify(numV/det)
Vxx_css=sp.sympify(d['Vxx']).subs(rr,1)
print("s-deriv=0 reproduces background Vxx?", sp.simplify(Vxx.subs(rr,1)-Vxx_css)==0, flush=True)
print("wrote fluid_sderiv.pkl", flush=True)
