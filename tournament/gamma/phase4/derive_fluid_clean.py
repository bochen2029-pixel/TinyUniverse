# Clean derivation of the CSS fluid ODEs from nabla_a T^{ab}=0, avoiding complex logs.
# Use coordinates (tau, r) with tau=-t>0 so the center is r->0, similarity var e^x=r/tau (x real).
# metric ds^2 = -alpha^2 dtau^2 + a^2 dr^2 + r^2 dOmega^2 (dt^2=dtau^2).
# CSS ansatz: alpha=alpha(x), a=a(x), V=V(x), and rho= omega(x)/(4 pi r^2 a^2), with
#   N = alpha a^{-1} e^{-x} = alpha/(a) * (tau/r).  x=ln(r/tau).
# We will (1) compute nabla_a T^{ab}=0 in (tau,r,theta,phi), (2) substitute the ansatz where
# every field is F(x(tau,r)) with x=ln(r/tau) (REAL), (3) reduce to ODEs in x, (4) rewrite in
# KHA variables N,A=a^2,omega,V and COMPARE to KHA eqs (18) fluid part.
import sympy as sp
tau,r,th,ph=sp.symbols('tau r theta phi', positive=True)
x=sp.log(r/tau)   # real
A=sp.Function('A'); N=sp.Function('N'); om=sp.Function('om'); V=sp.Function('V')
Af=A(x); Nf=N(x); omf=om(x); Vf=V(x)
a_=sp.sqrt(Af)
alpha_=Nf*a_*sp.exp(x)          # alpha = N a e^{x}  (since N=alpha a^{-1} e^{-x})
rho_=omf/(4*sp.pi*r**2*Af)
p_=rho_/3
coords=[tau,r,th,ph]
g=sp.diag(-alpha_**2, a_**2, r**2, r**2*sp.sin(th)**2)
gi=g.inv()
n=4
def christ():
    G=[[[sp.S.Zero]*n for _ in range(n)] for _ in range(n)]
    for l in range(n):
        for m in range(n):
            for k in range(m,n):
                ss=sp.S.Zero
                for d in range(n):
                    gld=gi[l,d]
                    if gld==0: continue
                    ss+=gld*(sp.diff(g[d,m],coords[k])+sp.diff(g[d,k],coords[m])-sp.diff(g[m,k],coords[d]))
                v=sp.simplify(ss/2); G[l][m][k]=v; G[l][k][m]=v
    return G
print("christoffels..."); G=christ()
W=1/sp.sqrt(1-Vf**2)
# u^a: u^tau=W/alpha, u^r=W V/a
u=[W/alpha_, W*Vf/a_, sp.S.Zero, sp.S.Zero]
T=[[ (rho_+p_)*u[i]*u[j]+p_*gi[i,j] for j in range(n)] for i in range(n)]
def divT(b):
    ss=sp.S.Zero
    for ai in range(n):
        ss+=sp.diff(T[ai][b],coords[ai])
        for c in range(n):
            if G[ai][ai][c]!=0: ss+=G[ai][ai][c]*T[c][b]
            if G[b][ai][c]!=0: ss+=G[b][ai][c]*T[ai][c]
    return ss
print("div tau..."); Dtau=divT(0)
print("div r...");   Dr=divT(1)
# reduce: x=ln(r/tau); dx/dtau=-1/tau, dx/dr=1/r. Replace F(x) derivatives.
Ax,Nx,ox,Vx,A0,N0,o0,V0=sp.symbols('A_x N_x om_x V_x A0 N0 om0 V0')
dxdtau=sp.Rational(-1,1)/tau; dxdr=sp.Rational(1,1)/r
def red(E):
    E=E.doit()
    subs={sp.Derivative(Af,tau):Ax*dxdtau, sp.Derivative(Nf,tau):Nx*dxdtau, sp.Derivative(omf,tau):ox*dxdtau, sp.Derivative(Vf,tau):Vx*dxdtau,
          sp.Derivative(Af,r):Ax*dxdr, sp.Derivative(Nf,r):Nx*dxdr, sp.Derivative(omf,r):ox*dxdr, sp.Derivative(Vf,r):Vx*dxdr}
    E=E.subs(subs).subs({Af:A0,Nf:N0,omf:o0,Vf:V0})
    return sp.simplify(E)
print("reduce...")
Dtau_r=red(Dtau); Dr_r=red(Dr)
# multiply by convenient factors to clear r,tau. e^x=r/tau. Introduce E=e^x=r/tau. Replace r/tau->E.
E=sp.symbols('E', positive=True)
def clean(D):
    D=sp.simplify(D*4*sp.pi*r**2)      # clear the rho denominator
    D=sp.simplify(D)
    # replace tau = r/E
    D=D.subs(tau, r/E)
    D=sp.simplify(D)
    return D
Dtau_c=clean(Dtau_r); Dr_c=clean(Dr_r)
print("\n=== reduced conservation (energy, tau-component) *4pi r^2, in terms of A0,N0,om0,V0,derivs,E ===")
sp.pprint(sp.collect(sp.expand(Dtau_c),[Ax,Nx,ox,Vx]))
print("\n=== reduced conservation (momentum, r-component) ===")
sp.pprint(sp.collect(sp.expand(Dr_c),[Ax,Nx,ox,Vx]))
import pickle
pickle.dump((sp.srepr(Dtau_c),sp.srepr(Dr_c)), open("fluid_clean.pkl","wb"))
