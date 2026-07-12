# Center expansion with M = N*z (z=e^x), so M is regular (finite, nonzero) at center.
# N = M/z, N' = z dN/dz. In terms of M: N = M/z, dN/dz = (M' z - M)/z^2 where M'=dM/dz.
# N'_x = z dN/dz = (M' z - M)/z = M' - M/z. And N'/N = (M'-M/z)/(M/z) = z M'/M - 1.
# The N ODE  N'/N = -2 + A - 2omega/3  becomes  z M'/M - 1 = -2 + A - 2omega/3
#   => z M'/M = -1 + A - 2omega/3.  (regular: at z=0, A=1, omega=0 -> zM'/M=0 -> M'~0, M->const. good)
import sympy as sp
z=sp.symbols('z', positive=True)
m0,m1,m2 = sp.symbols('m0 m1 m2')          # M = m0 + m1 z + m2 z^2
a2,a3    = sp.symbols('a2 a3')             # A = 1 + a2 z^2 + a3 z^3 (a1=0 by regularity/parity)
a1       = sp.symbols('a1')
v1,v2,v3 = sp.symbols('v1 v2 v3')          # V = v1 z + v2 z^2 + v3 z^3
w2,w3    = sp.symbols('w2 w3')             # omega = w2 z^2 + w3 z^3
Mser=m0+m1*z+m2*z**2
Aser=1+a1*z+a2*z**2+a3*z**3
Vser=v1*z+v2*z**2+v3*z**3
Wser=w2*z**2+w3*z**3
Nser=Mser/z

def dx(F): return z*sp.diff(F,z)   # d/dx = z d/dz

# A' ODE
Ap_ode = Aser*(1 - Aser + (2*Wser/(1-Vser**2))*(1 + Vser**2/3))
eqA = dx(Aser) - Ap_ode
# N' ODE via M: dx(N) with N=M/z: dx(M/z)= z d/dz(M/z)= z*(M' /z - M/z^2)= M' - M/z... but as series:
# easier: dx(Nser) directly (sympy handles M/z derivative), then eqN = dx(Nser)-Nser*(-2+Aser-2/3 Wser)
eqN = dx(Nser) - Nser*(-2+Aser-sp.Rational(2,3)*Wser)
eqN = sp.simplify(eqN*z)  # clear 1/z

# fluid 2x2
N,A,om,V=Nser,Aser,Wser,Vser
omV2=1-V**2
m11=(1+N*V); m12=4*(N+V)/(3*omV2); m21=(4*V+N+3*N*V**2); m22=4*(1+V**2+2*N*V)/omV2
Ap=Ap_ode; Np=Nser*(-2+Aser-sp.Rational(2,3)*Wser)
rhs1=-(-N*V*Ap/(3*A)+sp.Rational(4,3)*V*Np+2*N*(1+4*om/(9*omV2)))
rhs2=-( N*omV2*Ap/A+4*(1+V**2)*Np+2*N*(1+3*V**2))
det=m11*m22-m12*m21
numV=(m11*rhs2-rhs1*m21)
numU=(rhs1*m22-m12*rhs2)   # = u*det, u=omega'/omega
# equations: dx(V)*det = numV ; dx(omega)= (numU/det)*omega -> dx(omega)*det = numU*omega
eqV = sp.simplify((dx(Vser)*det - numV)*z**0)
eqW = sp.simplify((dx(Wser)*det - numU*om))

def coeffs(e, upto):
    e2=sp.series(e, z, 0, upto).removeO()
    e2=sp.expand(e2)
    p=sp.Poly(e2, z)
    return [c for c in p.all_coeffs() if c!=0]

eqs=[]
eqs += coeffs(eqA,6)
eqs += coeffs(eqN,6)
eqs += coeffs(eqV,7)
eqs += coeffs(eqW,8)

unk=[m0,m1,m2,a1,a2,a3,v1,v2,v3,w2,w3]
sol=sp.solve(eqs, unk, dict=True)
print("branches:", len(sol))
for s in sol:
    d={str(k):sp.nsimplify(sp.simplify(v)) for k,v in s.items()}
    used=set(s.keys()); free=[str(x) for x in unk if x not in used]
    print("free:", free, " ->")
    for k in ['m0','m1','m2','a1','a2','a3','v1','v2','v3','w2','w3']:
        if k in d: print(f"    {k} = {d[k]}")
    print()
