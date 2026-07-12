# hka_sonic_diag.py — diagnose the order structure of the fluid pair at the sonic point.
import sympy as sp

g = sp.Rational(4, 3); cs = sp.sqrt(g - 1)
x = sp.symbols('x', real=True)
V0val = sp.Rational(-1, 4)
N0 = sp.nsimplify((1 - cs*V0val)/(cs - V0val))
A0 = sp.nsimplify((g**2+4*g-4+8*(g-1)**sp.Rational(3,2)*V0val-(3*g-2)*(2-g)*V0val**2)/(g**2*(1-V0val**2)))
om0 = sp.nsimplify(2*cs*(cs-V0val)*(1+cs*V0val)/(g**2*(1-V0val**2)))
V0 = sp.nsimplify(V0val)
print("A0,N0,om0,V0 =", float(A0), float(N0), float(om0), float(V0))

order = 3
A1,A2,A3 = sp.symbols('A1 A2 A3'); N1,N2,N3 = sp.symbols('N1 N2 N3')
O1,O2,O3 = sp.symbols('O1 O2 O3'); V1,V2,V3 = sp.symbols('V1 V2 V3')
Aser = A0+A1*x+A2*x**2+A3*x**3; Nser = N0+N1*x+N2*x**2+N3*x**3
Oser = om0+O1*x+O2*x**2+O3*x**3; Vser = V0+V1*x+V2*x**2+V3*x**3
Axs=sp.diff(Aser,x); Nxs=sp.diff(Nser,x); Oxs=sp.diff(Oser,x); Vxs=sp.diff(Vser,x)
omV2 = 1 - Vser**2

R_a = sp.expand((Axs - Aser*(1-Aser))*omV2 - Aser*2*Oser*(1+(g-1)*Vser**2))
R_b = sp.expand(Nxs - Nser*(-2+Aser-(2-g)*Oser))
RHS_d = (3*(2-g)/2)*Nser*Vser-((2+g)/2)*Aser*Nser*Vser+(2-g)*Nser*Vser*Oser
R_d = sp.expand((1+Nser*Vser)*Oxs*omV2 + g*(Nser+Vser)*Vxs*Oser - RHS_d*Oser*omV2)
RHS_e = (2-g)*(g-1)*Nser*Oser+((7*g-6)/2)*Nser+((2-3*g)/2)*Aser*Nser
R_e = sp.expand((g-1)*(Nser+Vser)*Oxs*omV2 + g*(1+Nser*Vser)*Vxs*Oser - RHS_e*Oser*omV2)

# metric eqs first, order 0: gives A1, N1
print("\n-- order 0 --")
for nm,R in [('R_a',R_a),('R_b',R_b),('R_d',R_d),('R_e',R_e)]:
    print(nm, "[x^0]:", sp.simplify(R.coeff(x,0)))
# solve A1,N1 from R_a,R_b at order 0
solAN = sp.solve([R_a.coeff(x,0), R_b.coeff(x,0)], [A1,N1], dict=True)[0]
print("A1,N1 =", {k:sp.simplify(v) for k,v in solAN.items()})
# Now fluid order 0 with A1,N1 substituted:
Rd0 = sp.simplify(R_d.coeff(x,0).subs(solAN)); Re0 = sp.simplify(R_e.coeff(x,0).subs(solAN))
print("R_d[x^0] after AN sub:", Rd0)
print("R_e[x^0] after AN sub:", Re0)
# these should be identically 0 (the L'Hopital condition already guarantees consistency at x=0).
# The (O1,V1) quadratic comes from order-1 fluid eqs. Check order 1:
print("\n-- order 1 (fluid), after A1,N1 sub --")
Rd1 = sp.expand(R_d.coeff(x,1).subs(solAN)); Re1 = sp.expand(R_e.coeff(x,1).subs(solAN))
print("R_d[x^1]:", Rd1)
print("R_e[x^1]:", Re1)
# Are these quadratic in O1,V1? show degrees
print("R_d[x^1] as poly in V1:", sp.Poly(Rd1, V1).degree() if Rd1.has(V1) else 0,
      " in O1:", sp.Poly(Rd1, O1).degree() if Rd1.has(O1) else 0)
# Solve the two order-1 fluid eqs for (O1,V1):
sol1 = sp.solve([Rd1, Re1], [O1, V1], dict=True)
print(f"\n{len(sol1)} solution(s) for (O1,V1):")
for si,s in enumerate(sol1):
    print(f"  sol {si}: O1={sp.simplify(s.get(O1))}={float(s.get(O1)):+.5f}  V1={sp.simplify(s.get(V1))}={float(s.get(V1)):+.5f}")
    print(f"           A1={float(solAN[A1].subs(s) if solAN[A1].has(O1,V1) else solAN[A1]):+.5f}  N1={float(solAN[N1].subs(s) if solAN[N1].has(O1,V1) else solAN[N1]):+.5f}")
