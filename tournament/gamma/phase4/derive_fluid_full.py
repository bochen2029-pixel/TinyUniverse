# FULL 4D derivation of nabla_a T^{ab}=0 including the angular (theta,phi) pressure terms that
# contribute to the radial momentum equation via Gamma^r_{theta theta} T^{theta theta} etc.
# This fixes the earlier 2D-sector omission. Manual chain rule (no logs).
import sympy as sp, pickle
tau,r,th=sp.symbols('tau r theta', positive=True)
A0,N0,o0,V0=sp.symbols('A0 N0 om0 V0')
Ax,Nx,ox,Vx=sp.symbols('A_x N_x om_x V_x')
E=r/tau
a=sp.sqrt(A0); alpha=N0*a*E; W=1/sp.sqrt(1-V0**2)
rho=o0/(4*sp.pi*r**2*A0); p=rho/3
# full 4D metric
g=sp.diag(-alpha**2, a**2, r**2, r**2*sp.sin(th)**2)
gi=g.inv()
coords=[tau,r,th,sp.symbols('phi')]
# chain-rule derivative wrt coord index k
def Dk(expr,k):
    e=sp.diff(expr,coords[k])
    fac = (-1/tau) if k==0 else ((1/r) if k==1 else 0)  # dx/dtau=-1/tau, dx/dr=1/r, dx/dth=0
    e+= (sp.diff(expr,A0)*Ax+sp.diff(expr,N0)*Nx+sp.diff(expr,o0)*ox+sp.diff(expr,V0)*Vx)*fac
    return e
def christ(l,i,j):
    ss=sp.S.Zero
    for dd in range(4):
        if gi[l,dd]==0: continue
        ss+=gi[l,dd]*(Dk(g[dd,i],j)+Dk(g[dd,j],i)-Dk(g[i,j],dd))
    return sp.simplify(ss/2)
uu=[W/alpha, W*V0/a, sp.S.Zero, sp.S.Zero]  # u^a, correct normalization
T=[[ (rho+p)*uu[i]*uu[j]+p*gi[i,j] for j in range(4)] for i in range(4)]
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
print("computing full 4D divergence (b=tau, energy)...")
Dtau=divT(0)
print("computing full 4D divergence (b=r, radial momentum, WITH angular pressure terms)...")
Dr=divT(1)
# reduce: tau=r/E, express via E, clear rho denom
Es=sp.symbols('E',positive=True)
def fin(D):
    D=sp.simplify(D*4*sp.pi*r**2)
    D=D.subs(tau,r/Es)
    return sp.simplify(D)
en=fin(Dtau); mo=fin(Dr)
# check en,mo depend on E? (should be overall factor)
print("energy has E?",en.has(Es)," momentum has E?",mo.has(Es))
pickle.dump((sp.srepr(en),sp.srepr(mo)),open("fluid_full.pkl","wb"))
# quick: solve for (ox,Vx) with metric slopes substituted, check CENTER dV/dx.
Axv=A0*(1-A0+(2*o0/(1-V0**2))*(1+V0**2/3)); Nxv=N0*(-2+A0-sp.Rational(2,3)*o0)
sol=sp.solve([en.subs({Ax:Axv,Nx:Nxv}),mo.subs({Ax:Axv,Nx:Nxv})],[ox,Vx],dict=True)
if sol:
    Vxx=sp.simplify(sol[0][Vx].subs(Es,1)); omx=sp.simplify(sol[0][ox].subs(Es,1))
    fV=sp.lambdify((A0,N0,o0,V0),Vxx,'numpy'); fom=sp.lambdify((A0,N0,o0,V0),omx,'numpy')
    import numpy as np
    print("CENTER dV/dx (A=1,V=0,om->0,N=1e5):", float(fV(1.0,1e5,1e-9,0.0)))
    print("CENTER dom/dx:", float(fom(1.0,1e5,1e-9,0.0)))
    # sonic locus
    e1=sp.expand(en.subs({Ax:Axv,Nx:Nxv})); e2=sp.expand(mo.subs({Ax:Axv,Nx:Nxv}))
    det=sp.simplify((e1.coeff(ox)*e2.coeff(Vx)-e1.coeff(Vx)*e2.coeff(ox)))
    print("sonic det factor:", sp.factor(sp.numer(sp.together(det.subs(Es,1)))))
    pickle.dump({'omx':sp.srepr(omx),'Vxx':sp.srepr(Vxx)},open("rhs_full.pkl","wb"))
else:
    print("no solution")
