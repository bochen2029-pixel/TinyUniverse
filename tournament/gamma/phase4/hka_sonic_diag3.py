# Exactly characterize the order-1 fluid quadratic at the sonic point.
import numpy as np, sympy as sp
import hka_background as H

V0 = -0.25
A0,N0,om0,_ = H.sonic_data(V0)
g = sp.Rational(4,3)
o1,v1,a1,n1 = sp.symbols('o1 v1 a1 n1', real=True)
x = sp.symbols('x')
# series to order 1
A = sp.Float(A0,20) + a1*x; N = sp.Float(N0,20) + n1*x
O = sp.Float(om0,20) + o1*x; V = sp.Float(V0,20) + v1*x
Ax=sp.diff(A,x); Nx=sp.diff(N,x); Ox=sp.diff(O,x); Vx=sp.diff(V,x)
oV2 = 1-V**2
R_a = sp.expand((Ax - A*(1-A))*oV2 - A*2*O*(1+(g-1)*V**2))
R_b = sp.expand(Nx - N*(-2+A-(2-g)*O))
RHS_d = (3*(2-g)/2)*N*V-((2+g)/2)*A*N*V+(2-g)*N*V*O
R_d = sp.expand((1+N*V)*Ox*oV2 + g*(N+V)*Vx*O - RHS_d*O*oV2)
RHS_e = (2-g)*(g-1)*N*O+((7*g-6)/2)*N+((2-3*g)/2)*A*N
R_e = sp.expand((g-1)*(N+V)*Ox*oV2 + g*(1+N*V)*Vx*O - RHS_e*O*oV2)

# order 0 metric -> a1, n1
ra0 = R_a.coeff(x,0); rb0 = R_b.coeff(x,0)
solAN = sp.solve([ra0,rb0],[a1,n1],dict=True)[0]
print("a1,n1 =", {k:float(v) for k,v in solAN.items()})
# order 0 fluid, with a1,n1 substituted -> quadratic in o1,v1
rd0 = sp.expand(R_d.coeff(x,0).subs(solAN))
re0 = sp.expand(R_e.coeff(x,0).subs(solAN))
print("\nR_d[x^0] =", rd0)
print("R_e[x^0] =", re0)
print("\ndeg in o1:", sp.Poly(rd0,o1).total_degree() if rd0.has(o1) else 0,
      " deg in v1:", sp.Poly(rd0,v1).total_degree() if rd0.has(v1) else 0)
# solve exactly
sols = sp.solve([rd0,re0],[o1,v1],dict=True)
print(f"\n{len(sols)} EXACT roots (o1,v1):")
for s in sols:
    oo=complex(s[o1]); vv=complex(s[v1])
    print(f"  o1={oo}  v1={vv}")
