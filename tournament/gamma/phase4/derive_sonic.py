# Derivation of the sonic-point regularity condition for the radiation-fluid CSS,
# from the KHA95 (gr-qc/9503007) autonomous ODE system transcribed in equations.md §1.2.
#
# We treat the two fluid equations (b) with all _,s -> 0 as a 2x2 linear system for the
# unknowns (u, v) = (omega_,x / omega, V_,x), substitute the KNOWN metric slopes A_,x, N_,x
# from (a), then:
#   (1) form the coefficient matrix M and RHS vector b  (M [u v]^T = b),
#   (2) det(M) = 0 is the sonic locus  -> compare to the pre-derived 3 - V^2 - N^2 + 4NV + 3N^2V^2,
#   (3) regularity (removable singularity): the numerators from Cramer's rule must also vanish
#       when det=0. That gives V'(0), omega'(0) via L'Hopital.
#
# No fitting, no target injected: pure algebra from the transcribed equations.

import sympy as sp

N, A, om, V = sp.symbols('N A omega V', positive=False, real=True)

# metric slopes (a) -- exact, already in dX/dx form
Ax = A*(1 - A + (2*om/(1-V**2))*(1 + V**2/3))        # A_,x
Nx = N*(-2 + A - sp.Rational(2,3)*om)                # N_,x

# unknowns
u = sp.symbols('u')   # = omega_,x / omega
v = sp.symbols('v')   # = V_,x

# ---- fluid eq 1 (b, first), _,s = 0 ----
# (1+NV) u + (4(N+V)/(3(1-V^2))) v - (N V/(3A)) A_,x + (4V/3) N_,x + 2N(1 + 4om/(9(1-V^2))) = 0
eq1 = ( (1+N*V)*u
        + (4*(N+V)/(3*(1-V**2)))*v
        - (N*V/(3*A))*Ax
        + (sp.Rational(4,3)*V)*Nx
        + 2*N*(1 + 4*om/(9*(1-V**2))) )

# ---- fluid eq 2 (b, second), _,s = 0 ----
# (4V+N+3NV^2) u + (4(1+V^2+2NV)/(1-V^2)) v + (N(1-V^2)/A) A_,x + 4(1+V^2) N_,x + 2N(1+3V^2) = 0
eq2 = ( (4*V + N + 3*N*V**2)*u
        + (4*(1+V**2+2*N*V)/(1-V**2))*v
        + (N*(1-V**2)/A)*Ax
        + 4*(1+V**2)*Nx
        + 2*N*(1+3*V**2) )

# coefficient matrix M and RHS b:  M [u v]^T = b   (move constants to RHS)
a11 = sp.simplify(eq1.coeff(u)); a12 = sp.simplify(eq1.coeff(v))
a21 = sp.simplify(eq2.coeff(u)); a22 = sp.simplify(eq2.coeff(v))
b1  = sp.simplify(-(eq1.subs({u:0, v:0})))
b2  = sp.simplify(-(eq2.subs({u:0, v:0})))

M = sp.Matrix([[a11, a12],[a21, a22]])
print("M =")
sp.pprint(M)
detM = sp.simplify(M.det())
print("\ndet(M) =")
sp.pprint(detM)

# The det carries an overall 1/(1-V^2) from the v-column. Multiply through by (1-V^2)
# to get the polynomial locus:
det_poly = sp.simplify(detM*(1-V**2))
print("\ndet(M)*(1-V^2) =")
sp.pprint(sp.expand(det_poly))

# pre-derived candidate:
cand = 3 - V**2 - N**2 + 4*N*V + 3*N**2*V**2
print("\nPre-derived candidate  3 - V^2 - N^2 + 4NV + 3N^2V^2 :")
sp.pprint(cand)

# Is det_poly a scalar multiple of the candidate? Compute ratio.
ratio = sp.simplify(det_poly/cand)
print("\n det_poly / candidate  =")
sp.pprint(ratio)

# also print det_poly factored
print("\n factor(det_poly):")
sp.pprint(sp.factor(det_poly))
print("\n factor(candidate):")
sp.pprint(sp.factor(cand))
