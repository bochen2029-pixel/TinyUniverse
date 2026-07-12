# Manual-chain-rule derivation: never let sympy see log. Each field F(x) is represented by a
# symbol F0 and its x-derivative Fx; spatial/time derivatives use dx/dtau=-1/tau, dx/dr=1/r.
# We build the relativistic Euler equations directly:
#   energy:   u^a d_a rho + (rho+p) theta = 0,  theta=div u = (1/sqrt(-g)) d_a(sqrt(-g) u^a)
#   momentum: (rho+p) u^a d_a u^b + h^{ab} d_a p = 0, h^{ab}=g^{ab}+u^a u^b
# with rho=om/(4 pi r^2 A), p=rho/3, u^tau=W/alpha, u^r=W V/a, alpha=N a e^x, a=sqrt(A), e^x=r/tau.
import sympy as sp
tau,r=sp.symbols('tau r', positive=True)
A0,N0,o0,V0=sp.symbols('A0 N0 om0 V0', positive=False)
Ax,Nx,ox,Vx=sp.symbols('A_x N_x om_x V_x')
E=r/tau   # e^x
# fields and their partials via chain rule. F(x): d/dtau = Fx*(-1/tau), d/dr=Fx*(1/r)
class Fld:
    def __init__(s,val,dv): s.v=val; s.dv=dv    # value, x-derivative (symbol)
def dtau(v,dvdx): return dvdx*(-1/tau)
def dr(v,dvdx): return dvdx*(1/r)
A=A0; a=sp.sqrt(A0); alpha=N0*a*E; W=1/sp.sqrt(1-V0**2)
rho=o0/(4*sp.pi*r**2*A0); p=rho/3
# derivatives of composite quantities via product/chain rule, expressed with Ax,Nx,ox,Vx.
# We'll compute needed d_tau, d_r of: rho, p, alpha, a, u^tau, u^r.
def Dtau(expr):
    # total tau-derivative: explicit tau/r partials + field x-derivs
    e=sp.diff(expr,tau)  # explicit (treats A0 etc constant, tau,r explicit)
    e+= sp.diff(expr,A0)*Ax*(-1/tau)+sp.diff(expr,N0)*Nx*(-1/tau)+sp.diff(expr,o0)*ox*(-1/tau)+sp.diff(expr,V0)*Vx*(-1/tau)
    return e
def Dr(expr):
    e=sp.diff(expr,r)
    e+= sp.diff(expr,A0)*Ax*(1/r)+sp.diff(expr,N0)*Nx*(1/r)+sp.diff(expr,o0)*ox*(1/r)+sp.diff(expr,V0)*Vx*(1/r)
    return e
uT=W/alpha; uR=W*V0/a
# sqrt(-g)=alpha*a*r^2*sin(theta); drop sin(theta) (theta-independent). SG=alpha*a*r^2
SG=alpha*a*r**2
theta=(Dtau(SG*uT)+Dr(SG*uR))/SG          # divergence of u
# energy eq:
energy=uT*Dtau(rho)+uR*Dr(rho)+(rho+p)*theta
# momentum (b=r): (rho+p) u^a d_a u^r + h^{r a} d_a p + (rho+p)Gamma^r_{ab}u^a u^b ... need
# the geodesic-like term with Christoffels for u^a nabla_a u^r. Use full:
#   u^a nabla_a u^b = u^a d_a u^b + Gamma^b_{ac} u^a u^c
# Build Christoffel^r and ^tau. metric diag(-alpha^2,a^2,r^2,r^2 sin^2).
g=sp.diag(-alpha**2,a**2,r**2,r**2)  # sin^2 dropped (equatorial; fine for radial eq)
coords=[tau,r]
# only need Gamma^r_{tau tau},^r_{rr},^r_{tau r} and Gamma^tau_{..} for b=r and energy check.
def Gam(up,i,j):
    # Christoffel with two lower i,j, upper 'up'; use 2D (tau,r) sector (theta,phi decouple for radial)
    gi=g.inv()
    ss=sp.S.Zero
    for d in range(2):
        if gi[up,d]==0: continue
        ss+=gi[up,d]*(Dcoord(g[d,i],j)+Dcoord(g[d,j],i)-Dcoord(g[i,j],d))
    return sp.simplify(ss/2)
def Dcoord(expr,k):  # derivative wrt coord k (0=tau,1=r) with chain rule
    return Dtau(expr) if k==0 else Dr(expr)
uu=[uT,uR]
# u^a nabla_a u^r:
acc_r=sum(uu[aidx]*(Dcoord(uR,aidx)) for aidx in range(2))
for aidx in range(2):
    for cidx in range(2):
        acc_r+=Gam(1,aidx,cidx)*uu[aidx]*uu[cidx]
gi=g.inv()
hRa_dp=sum((gi[1,aidx]+uR*uu[aidx])*Dcoord(p,aidx) for aidx in range(2))
momentum=(rho+p)*acc_r+hRa_dp
# Now reduce: substitute tau=r/E, simplify, express in E. Multiply by clearing factors.
def fin(expr,mult):
    ex=sp.simplify(expr*mult)
    ex=ex.subs(tau,r/E)
    return sp.simplify(ex)
Es=sp.symbols('E',positive=True)
en=fin(energy, 4*sp.pi*r**2).subs(E,Es)
mo=fin(momentum, 4*sp.pi*r**2).subs(E,Es)
en=sp.simplify(en); mo=sp.simplify(mo)
print("ENERGY eq (=0), collected in field x-derivs:")
sp.pprint(sp.collect(sp.expand(en),[Ax,Nx,ox,Vx]))
print("\nMOMENTUM(r) eq (=0):")
sp.pprint(sp.collect(sp.expand(mo),[Ax,Nx,ox,Vx]))
import pickle; pickle.dump((sp.srepr(en),sp.srepr(mo)),open("fluid_manual.pkl","wb"))
