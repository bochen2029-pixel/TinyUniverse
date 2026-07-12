# Derive the regular-center expansion (x -> -inf) for the KHA autonomous system.
# At the center r=0: A=a^2=1 (regular), V=0 (no flow at center), and the density is finite
# so rho(center) finite => omega = 4 pi r^2 a^2 rho -> 0 like r^2 -> 0. Also N = alpha/a * e^{-x}.
# Since r = e^{x-s}, at fixed s, x->-inf is r->0. omega ~ 4pi r^2 rho_c ~ C e^{2(x-s)} but we've
# set _,s=0 so treat omega=omega(x): omega ~ omega2 e^{2x} near center (leading).
#
# Let's posit near-center (xi=x->-inf) expansions in e^{x} (call it z=e^{x}, z->0+):
#   A = 1 + a2 z^2 + ...
#   V = v1 z + ...        (velocity ~ radius, ingoing => v1<0)
#   omega = w2 z^2 + ...  (density finite => omega~r^2)
#   N = n0 + n1 z + ...   ?  but N ~ alpha/a e^{-x} = alpha/a / z ... that DIVERGES as z->0.
# So N ~ N_{-1}/z? Let's check with the ODE N'/N = -2 + A - 2omega/3. As z->0: A->1, omega->0,
# so N'/N -> -1. But d/dx of ln N: if N ~ c z^p = c e^{px}, then N'/N = p. So p=-1 => N ~ c e^{-x}
# = c/z -> DIVERGES. Good, consistent: N ~ N_c e^{-x} near center, i.e. alpha/a -> N_c const.
# So define n_{-1}: N = n_{-1} e^{-x}(1 + ...). Let's do the expansion with sympy.

import sympy as sp
z=sp.symbols('z', positive=True)   # z=e^x -> 0 at center
# ansatz
nm1,n1,n2 = sp.symbols('nm1 n1 n2')
a2,a3,a4  = sp.symbols('a2 a3 a4')
v1,v2,v3  = sp.symbols('v1 v2 v3')
w2,w3,w4  = sp.symbols('w2 w3 w4')
Nser = nm1/z*(1 + n1*z + n2*z**2)
Aser = 1 + a2*z**2 + a3*z**3
Vser = v1*z + v2*z**2 + v3*z**3
Wser = w2*z**2 + w3*z**3 + w4*z**4    # omega

# ODEs in x; note d/dx = z d/dz. So F' = z dF/dz.
def dx(F): return z*sp.diff(F,z)

Ap_ode = Aser*(1 - Aser + (2*Wser/(1-Vser**2))*(1 + Vser**2/3))
Np_ode = Nser*(-2 + Aser - sp.Rational(2,3)*Wser)

# fluid 2x2 -> omega', V'. Build and solve.
N,A,om,V=Nser,Aser,Wser,Vser
omV2=1-V**2
m11=(1+N*V); m12=4*(N+V)/(3*omV2); m21=(4*V+N+3*N*V**2); m22=4*(1+V**2+2*N*V)/omV2
Ap=Ap_ode; Np=Np_ode
rhs1=-(-N*V*Ap/(3*A)+sp.Rational(4,3)*V*Np+2*N*(1+4*om/(9*omV2)))
rhs2=-( N*omV2*Ap/A+4*(1+V**2)*Np+2*N*(1+3*V**2))
det=m11*m22-m12*m21
u=(rhs1*m22-m12*rhs2)/det   # = omega'/omega
Vp=(m11*rhs2-rhs1*m21)/det
omp=u*om

# match orders: dx(A) == Ap ; dx(N)==Np ; dx(omega)==omp ; dx(V)==Vp, order by order in z.
eqA = sp.series(dx(Aser)-Ap_ode, z, 0, 5).removeO()
eqN = sp.series(dx(Nser)-Np_ode, z, 0, 4).removeO()
# fluid: multiply by det to avoid series of ratio blowups; match numerator
eqV = sp.series(dx(Vser)*det-(m11*rhs2-rhs1*m21), z, 0, 5).removeO()
eqW = sp.series(dx(Wser)*det-(u*om*det), z, 0, 6).removeO()

# collect coefficients and solve sequentially
unk=[nm1,n1,n2,a2,a3,v1,v2,v3,w2,w3,w4]
eqs=[]
for e in [eqA,eqN,eqV,eqW]:
    p=sp.Poly(sp.expand(e),z)
    for c in p.all_coeffs():
        if c!=0: eqs.append(c)
sol=sp.solve(eqs, unk, dict=True)
print("center-expansion solutions (parameterized):")
for s in sol[:4]:
    print({str(k):sp.simplify(v) for k,v in s.items()})
print("\nnumber of solution branches:", len(sol))
print("free symbols remaining (the shoot parameter):")
for s in sol[:2]:
    used=set(s.keys()); free=[str(x) for x in unk if x not in used]
    print("  branch:", free)
