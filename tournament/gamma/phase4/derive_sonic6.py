# Part 6: CORRECT regularity via Fredholm solvability.
# At det(M)=0, M is rank 1 (generically). M (u,v)^T = b has finite solutions iff b lies
# in Im(M), i.e. l . b = 0 where l is the LEFT null vector (l M = 0). That is the SINGLE
# regularity/solvability condition at the sonic point. Combined with D=0 (the locus) it
# leaves a 1-parameter family. Then the analytic slope is fixed by L'Hopital.

import sympy as sp
N,A,om,V,u,v = sp.symbols('N A omega V u v', real=True)
Ax = A*(1 - A + (2*om/(1-V**2))*(1 + V**2/3))
Nx = N*(-2 + A - sp.Rational(2,3)*om)
eq1 = ((1+N*V)*u + (4*(N+V)/(3*(1-V**2)))*v - (N*V/(3*A))*Ax + (sp.Rational(4,3)*V)*Nx + 2*N*(1 + 4*om/(9*(1-V**2))))
eq2 = ((4*V+N+3*N*V**2)*u + (4*(1+V**2+2*N*V)/(1-V**2))*v + (N*(1-V**2)/A)*Ax + 4*(1+V**2)*Nx + 2*N*(1+3*V**2))
a11=eq1.coeff(u);a12=eq1.coeff(v);a21=eq2.coeff(u);a22=eq2.coeff(v)
b1=-(eq1.subs({u:0,v:0}));b2=-(eq2.subs({u:0,v:0}))
M = sp.Matrix([[a11,a12],[a21,a22]])
bvec = sp.Matrix([b1,b2])

# left null vector l=(l1,l2): l1 a11 + l2 a21 = 0, l1 a12 + l2 a22 = 0.
# take l = (a21, -a11) ... check: l.col? For rank1, l=(-a21,a11) gives l.M=0? l M row =
# (-a21 a11 + a11 a21, -a21 a12 + a11 a22) = (0, det). det=0 so l=( -a21, a11 ) is left null.
# equivalently l=(a22,-a12) works on the other pairing when det=0. Use l=(-a21,a11):
l = sp.Matrix([[-a21, a11]])
solvability = sp.simplify((l*bvec)[0])   # = 0 is the regularity condition
print("Fredholm solvability  l.b  (with l=(-a21,a11)):")
sv = sp.simplify(solvability)
# clear denominators
sv_num = sp.numer(sp.together(sv))
sv_num = sp.expand(sv_num)
print("  numerator degrees: in A", sp.degree(sp.Poly(sv_num,A)), " in omega", sp.degree(sp.Poly(sv_num,om)))

# On the sonic locus D=0 substitute N=branch1 and see the condition on (A,omega,V):
Nbr1 = (-2*V + sp.sqrt(3)*(V**2-1))/(3*V**2-1)
sv_b1 = sp.simplify(sv.subs(N,Nbr1))
print("\n solvability on branch1 (as relation among A,omega,V):")
# solve for omega in terms of A,V
sol_om = sp.solve(sp.numer(sp.together(sv_b1)), om)
print(" omega(A,V) from solvability:", [sp.simplify(s) for s in sol_om])

# BUT we still have TWO constraint eqs (a) that must hold at x=0 relating A,omega to A_,x,N_,x?
# No -- eqs (a) DEFINE A_,x,N_,x, they don't constrain the state. The state (A0,omega0,V0) with
# N0=N(V0) has ONE more relation from solvability -> so free params = V0 alone? Let's count:
#   unknowns at x=0: (N0,A0,omega0,V0) = 4
#   conditions: D=0 (1), solvability (1) => 2. Free = 2. Hmm that's a 2-param family, not 1.
# The resolution: the CENTER BC (A->1,V->0 as x->-infty) is the SECOND shooting condition.
# So at the sonic point we have a 2-parameter family (e.g. V0 and A0), and TWO center
# conditions (A->1 AND V->0) pin both. That's the correct count! KHA phrase it as "one
# parameter Vss(0)" because they ALSO use the gauge/one constraint differently. Let's just
# verify: with solvability giving omega0(A0,V0), we shoot in (A0,V0) to hit (A=1,V=0) center.
#
# Actually simpler & standard (Ori-Piran, Harada): at the sonic point, given V0, the
# regularity (solvability) + the requirement that the SLOPE be the analytic (nodal) one
# fixes A0 and omega0. Let's test whether solvability + a slope condition pins A0.
#
# Cleanest: treat (A0,omega0) both unknown, V0 the shoot param. Solvability is 1 eq.
# The analytic-slope (L'Hopital double root real) is the 2nd selection at the node. Then
# shoot V0 to center. Let's first just get omega0(A0,V0) and move the rest to the C++.
for vv in [-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3,0.4]:
    om_expr = sol_om[0] if sol_om else None
    if om_expr is not None:
        f = sp.lambdify((A,V), om_expr, 'numpy')
        # try A around 1
        import numpy as np
        for Aval in [0.8,1.0,1.2]:
            val = float(f(Aval,vv))
            print(f"   V0={vv:+.2f} A0={Aval:.1f} -> omega0={val:+.5f}")
