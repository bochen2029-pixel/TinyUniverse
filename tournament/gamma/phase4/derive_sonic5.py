# Part 5: correct treatment. At the sonic point x=0, det D=0. Ori-Piran: the solution
# through a regular singular point is analytic iff the numerators vanish there too, BUT
# num_u=0 and num_v=0 are NOT both independent from D=0 -- one of them is implied. The
# right statement: on the locus D=0, the two numerators are proportional (rank drops).
# So imposing num_v=0 (say) is ONE condition; together with D=0 that's 2 conditions on
# the 3 remaining unknowns (A,omega,V) after N=N(V) -> a 1-parameter family. GOOD.
#
# Let's TEST that claim: substitute N=branch into num_u,num_v; are they proportional
# (as functions of A,omega) with a V-dependent ratio? If yes, imposing one suffices.

import sympy as sp
N,A,om,V,u,v = sp.symbols('N A omega V u v', real=True)
Ax = A*(1 - A + (2*om/(1-V**2))*(1 + V**2/3))
Nx = N*(-2 + A - sp.Rational(2,3)*om)
eq1 = ((1+N*V)*u + (4*(N+V)/(3*(1-V**2)))*v - (N*V/(3*A))*Ax + (sp.Rational(4,3)*V)*Nx + 2*N*(1 + 4*om/(9*(1-V**2))))
eq2 = ((4*V+N+3*N*V**2)*u + (4*(1+V**2+2*N*V)/(1-V**2))*v + (N*(1-V**2)/A)*Ax + 4*(1+V**2)*Nx + 2*N*(1+3*V**2))
a11=eq1.coeff(u);a12=eq1.coeff(v);a21=eq2.coeff(u);a22=eq2.coeff(v)
b1=-(eq1.subs({u:0,v:0}));b2=-(eq2.subs({u:0,v:0}))
num_u = sp.Matrix([[b1,a12],[b2,a22]]).det()
num_v = sp.Matrix([[a11,b1],[a21,b2]]).det()

Nbr1 = (-2*V + sp.sqrt(3)*(V**2-1))/(3*V**2-1)
Nbr0 = (-2*V - sp.sqrt(3)*(V**2-1))/(3*V**2-1)

for name,Nbr in [('branch1',Nbr1),('branch0',Nbr0)]:
    nu = sp.simplify(num_u.subs(N,Nbr))
    nv = sp.simplify(num_v.subs(N,Nbr))
    # ratio nu/nv should be independent of A,omega if they're proportional
    ratio = sp.simplify(nu/nv)
    has_A = ratio.has(A); has_om = ratio.has(om)
    print(f"{name}: num_u/num_v depends on A? {has_A}  on omega? {has_om}")
    if not has_A and not has_om:
        print("   -> PROPORTIONAL on the sonic locus. ratio(V) =", sp.simplify(ratio))
    # Solve num_v=0 for omega given (A,V,N=branch): linear in omega
    sol_om = sp.solve(sp.numer(sp.together(nv)), om)
    print(f"   omega from num_v=0 (function of A,V): {[sp.simplify(s) for s in sol_om]}")
    print()
