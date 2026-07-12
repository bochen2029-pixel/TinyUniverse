# Part 3: pin the sonic point completely.
# At the sonic point (x=0) we have FOUR unknowns (N0, A0, omega0, V0) and the conditions:
#   (i)  det D(N0,V0) = 0                         [sonic locus]
#   (ii) regularity: num_v(N0,A0,omega0,V0)=0     [removable singularity, L'Hopital]
#   (iii) num_u(...)=0                            [second removable-singularity condition]
# BUT num_u=0 and num_v=0 are NOT independent at D=0 (a regular singular point has ONE
# analytic branch through it after the L'Hopital derivative is chosen). Standard result
# (Ori-Piran / Brady / Gundlach): at the sonic point, given the locus D=0 plus a chosen
# root of the *quadratic* that determines the through-going slope, (A0, omega0) get fixed
# by the two constraint equations (a) once you also demand analyticity. The clean count:
#   free params at sonic point after imposing D=0 and regularity = ONE (call it V0).
#   Then A0, omega0 are DETERMINED, N0 = N0(V0) from D=0.
#
# We reduce the count concretely.

import sympy as sp

N, A, om, V = sp.symbols('N A omega V', real=True)
u, v = sp.symbols('u v')

D = 3*N**2*V**2 - N**2 + 4*N*V - V**2 + 3

# solve D=0 for N (two branches); the PHYSICAL sonic point for the growing/critical
# solution is the one with N>0, and (per KHA fig / Ori-Piran) V in (-1,0) region near center.
solN = sp.solve(D, N)
print("N branches from D=0:")
for i,s in enumerate(solN):
    print(f"  branch {i}:", sp.simplify(s))

# Numerically explore: for a grid of V, the two N-roots
import numpy as np
Vg = np.linspace(-0.95, 0.95, 39)
f0 = sp.lambdify(V, solN[0], 'numpy')
f1 = sp.lambdify(V, solN[1], 'numpy')
print("\n  V      N_branch0     N_branch1")
for vv in [-0.9,-0.7,-0.5,-0.3,-0.1,0.1,0.3,0.5,0.7,0.9]:
    print(f"  {vv:+.2f}   {float(f0(vv)):+.5f}     {float(f1(vv)):+.5f}")

# Ori-Piran / EC critical solution: the sonic point sits at a specific (N,V). For the
# radiation fluid the known sonic-point data (from the literature; Ori-Piran 1990,
# reproduced in Harada-Maeda and in Gundlach-Martin-Garcia reviews) has the fluid
# subsonic->supersonic transition. We'll let the shoot FIND V0; here we just confirm
# the locus is smooth and single-signed for N in a neighborhood.

# The other two conditions (num_u=0, num_v=0) at D=0: substitute a branch of N and see
# what relation between (A,omega,V) they force. Because the ODE has a regular singular
# point, ONE combination of {num_u,num_v} vanishes automatically given the L'Hopital
# slope; the OTHER fixes a relation. Let's compute num_u, num_v symbolically, substitute
# N=branch, and reduce.

Ax = A*(1 - A + (2*om/(1-V**2))*(1 + V**2/3))
Nx = N*(-2 + A - sp.Rational(2,3)*om)
eq1 = ( (1+N*V)*u + (4*(N+V)/(3*(1-V**2)))*v - (N*V/(3*A))*Ax
        + (sp.Rational(4,3)*V)*Nx + 2*N*(1 + 4*om/(9*(1-V**2))) )
eq2 = ( (4*V + N + 3*N*V**2)*u + (4*(1+V**2+2*N*V)/(1-V**2))*v + (N*(1-V**2)/A)*Ax
        + 4*(1+V**2)*Nx + 2*N*(1+3*V**2) )
a11=eq1.coeff(u); a12=eq1.coeff(v); a21=eq2.coeff(u); a22=eq2.coeff(v)
b1=-(eq1.subs({u:0,v:0})); b2=-(eq2.subs({u:0,v:0}))
num_u = sp.Matrix([[b1,a12],[b2,a22]]).det()
num_v = sp.Matrix([[a11,b1],[a21,b2]]).det()

# regularity: numerators vanish. clear denominators:
NUM_U = sp.simplify(num_u*9*(V+1)/(4*N))          # strip constant/pole factors
NUM_V = sp.simplify(num_v*9/( -N*(V**2-1) ))
print("\nNUM_U (=0 at sonic, cleared):")
print(sp.expand(NUM_U))
print("\nNUM_V (=0 at sonic, cleared):")
print(sp.expand(NUM_V))

# These are two equations in (A,omega) linear in omega and A? check degrees:
print("\ndeg NUM_U in omega:", sp.degree(sp.Poly(sp.expand(NUM_U), om)))
print("deg NUM_U in A    :", sp.degree(sp.Poly(sp.expand(NUM_U), A)))
print("deg NUM_V in omega:", sp.degree(sp.Poly(sp.expand(NUM_V), om)))
print("deg NUM_V in A    :", sp.degree(sp.Poly(sp.expand(NUM_V), A)))
