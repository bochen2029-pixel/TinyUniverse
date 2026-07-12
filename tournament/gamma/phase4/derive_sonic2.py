# Part 2: (a) confirm det=0 <=> sound-speed condition; (b) L'Hopital limits at the
# removable singularity (numerators must vanish there -> gives V(0) locus + derivatives).

import sympy as sp

N, A, om, V = sp.symbols('N A omega V', real=True)
u, v = sp.symbols('u v')

Ax = A*(1 - A + (2*om/(1-V**2))*(1 + V**2/3))
Nx = N*(-2 + A - sp.Rational(2,3)*om)

eq1 = ( (1+N*V)*u + (4*(N+V)/(3*(1-V**2)))*v - (N*V/(3*A))*Ax
        + (sp.Rational(4,3)*V)*Nx + 2*N*(1 + 4*om/(9*(1-V**2))) )
eq2 = ( (4*V + N + 3*N*V**2)*u + (4*(1+V**2+2*N*V)/(1-V**2))*v + (N*(1-V**2)/A)*Ax
        + 4*(1+V**2)*Nx + 2*N*(1+3*V**2) )

a11 = eq1.coeff(u); a12 = eq1.coeff(v)
a21 = eq2.coeff(u); a22 = eq2.coeff(v)
b1  = -(eq1.subs({u:0, v:0}))
b2  = -(eq2.subs({u:0, v:0}))
M = sp.Matrix([[a11,a12],[a21,a22]])
detM = M.det()

# ------------------------------------------------------------------
# (a) sound speed condition. The sonic point: flow speed seen by the constant-x
# observer equals c_s = 1/sqrt(3). KHA's characteristic speeds for the fluid are
# (N + V)/(1 + NV) type combinations (the "N+V" and "NV" structure visible in M).
# The det factor is  D = 3 - V^2 - N^2 + 4NV + 3 N^2 V^2. Let's show D factors into
# the two characteristic conditions.
D = 3 - V**2 - N**2 + 4*N*V + 3*N**2*V**2
print("D = det locus:")
sp.pprint(sp.factor(D))

# The sound cone: characteristic speeds c_pm = (N + V +/- something). KHA / Ori-Piran:
# the sonic point is where the fluid 3-velocity relative to the similarity flow equals
# +- c_s. Solve D=0 for N in terms of V:
solN = sp.solve(D, N)
print("\nSonic locus, N(V) branches:")
for s in solN:
    sp.pprint(sp.simplify(s))

# Cross-check against the sound-speed statement. The local sound speed for p=rho/3 is
# c_s^2 = dp/drho = 1/3.  In KHA variables the relevant Lorentz-composition of the
# fluid velocity V with the "grid" velocity carried by N gives characteristics
#    lambda_pm : the sonic condition is  (V +/- c_s)/(1 +/- V c_s) = -N   (flow = -N carries
# the constant-x line). Let cs = 1/sqrt(3). Form both and see if their product = D (up to factor).
cs = 1/sp.sqrt(3)
# characteristic: the similarity-frame flow speed is N (the coefficient of omega_,x in the
# advective combination 1+NV etc). Sonic when the boosted sound speed matches:
lhs_p = (V + cs)/(1 + V*cs)   # forward sound ray velocity (special-rel velocity addition)
lhs_m = (V - cs)/(1 - V*cs)   # backward
cond_p = sp.simplify(lhs_p + N)   # = 0 ?
cond_m = sp.simplify(lhs_m + N)
prod = sp.simplify((cond_p*cond_m))
print("\n(V+cs)/(1+V cs) + N  times  (V-cs)/(1-V cs) + N  =")
sp.pprint(sp.simplify(prod))
# ratio to D:
print("\n that product / D :")
sp.pprint(sp.simplify(prod / D))
print("\n that product * (denominators) vs D — expand numerator of product*(1-3V^2)/... :")
num_prod = sp.simplify(prod * (1 + V*cs)*(1 - V*cs))
sp.pprint(sp.simplify(num_prod))
print(" times 3:")
sp.pprint(sp.expand(num_prod*3))

# ------------------------------------------------------------------
# (b) L'Hopital: at the sonic point det->0. For a regular (analytic) solution the
# Cramer numerators must ALSO ->0 there. Cramer:
#   u = det( [b1 a12; b2 a22] ) / detM ;  v = det( [a11 b1; a21 b2] ) / detM
num_u = sp.Matrix([[b1,a12],[b2,a22]]).det()
num_v = sp.Matrix([[a11,b1],[a21,b2]]).det()
print("\n\nnum_u (numerator for u=omega_,x/omega):")
sp.pprint(sp.simplify(num_u*(1-V**2)))   # clear the 1/(1-V^2)
print("\nnum_v (numerator for v=V_,x):")
sp.pprint(sp.simplify(num_v*(1-V**2)))

# The regularity condition at x=0: num_v = 0 (and num_u=0) simultaneously with D=0.
# Because A is fixed by the constraint algebra at the sonic point, this becomes a
# relation among (N, A, omega, V) at x=0. Print num_v=0 with det=0 substituted:
# Solve D=0 for A? A doesn't appear in D. A enters through num_u,num_v via Ax,Nx.
print("\n(regularity handled numerically in-tool via L'Hopital; symbolic forms above)")
