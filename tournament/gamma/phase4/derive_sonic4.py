# Part 4: solve (A0, omega0) from the 2x2 regularity system given (V0, N0=branch), and
# compute the L'Hopital slopes (omega', V') at the sonic point (the 0/0 limit).
# Emit a self-contained numerical function we can transcribe/verify against the C++.

import sympy as sp

N, A, om, V = sp.symbols('N A omega V', real=True)
u, v = sp.symbols('u v')

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

# regularity equations (numerators = 0). Multiply by (1-V^2)/... to clear poles -> polynomials
NU = sp.together(num_u)
NV = sp.together(num_v)
NU_n, _ = sp.fraction(NU)
NV_n, _ = sp.fraction(NV)
NU_n = sp.expand(NU_n)
NV_n = sp.expand(NV_n)

# These are linear in (A, omega). Solve the 2x2 for (A, omega) as functions of (N,V):
sol = sp.solve([NU_n, NV_n], [A, om], dict=True)
print("Number of (A,omega) solutions:", len(sol))
Asol = sp.simplify(sol[0][A])
omsol = sp.simplify(sol[0][om])
print("\nA0(N,V) =")
sp.pprint(Asol)
print("\nomega0(N,V) =")
sp.pprint(omsol)

# Now substitute the physical branch N = branch1 and evaluate numerically over V0.
Nbr = (-2*V + sp.sqrt(3)*(V**2-1))/(3*V**2-1)   # branch1 (N>0 region)
A_of_V = sp.lambdify(V, Asol.subs(N, Nbr), 'numpy')
om_of_V = sp.lambdify(V, omsol.subs(N, Nbr), 'numpy')
N_of_V = sp.lambdify(V, Nbr, 'numpy')

print("\n  V0      N0        A0         omega0")
for vv in [-0.6,-0.5,-0.4,-0.35,-0.3,-0.25,-0.2,-0.15,-0.1,-0.05]:
    try:
        print(f"  {vv:+.3f}   {float(N_of_V(vv)):+.5f}   {float(A_of_V(vv)):+.5f}   {float(om_of_V(vv)):+.5f}")
    except Exception as e:
        print(f"  {vv:+.3f}   err {e}")

# ---- L'Hopital slopes at the sonic point ----
# u = num_u/detM, v=num_v/detM both 0/0. L'Hopital in x:
#   at x=0, du/... -> we need d(num)/dx / d(det)/dx along the solution. Since everything
#   is a function of (N,A,omega,V) and we know N',A' exactly (eqs a) and want omega',V':
#   v = V' satisfies   V' = lim num_v/detM. Use L'Hopital: V' = (d num_v/dx)/(d detM/dx).
#   But d/dx = N' d/dN + A' d/dA + omega' d/domega + V' d/dV. This closes into a quadratic
#   for V' (Ori-Piran's two roots -> pick the analytic/through-going one).
# We'll DERIVE the quadratic for s := V'(0). Let w := omega'(0) = u*omega.
w = sp.symbols('w')  # omega'(0)
s = sp.symbols('s')  # V'(0)
# total x-derivative operator on a field F(N,A,omega,V):
def ddx(F):
    return sp.diff(F,N)*Nx + sp.diff(F,A)*Ax + sp.diff(F,om)*w + sp.diff(F,V)*s
# detM and num_v as functions (with u=w/omega, v=s):
detM = (a11*a22 - a12*a21)
# v = num_v/detM  => v*detM = num_v ; both ->0. differentiate: v'*detM + v*detM' = num_v'
# at sonic detM=0 so: v*detM' = num_v'  => s*ddx(detM) = ddx(num_v_with_uv_subst)
num_v_sub = num_v.subs({u: w/om, v: s})
num_u_sub = num_u.subs({u: w/om, v: s})
eqV = sp.Eq(s*ddx(detM), ddx(num_v_sub))
eqU = sp.Eq((w/om)*ddx(detM), ddx(num_u_sub))
print("\n(L'Hopital closure eqs eqV, eqU are quadratic in (s,w); solved in-tool numerically.)")
print("Symbolic degree check:")
print("  eqV poly in s:", sp.degree(sp.Poly(sp.expand((eqV.lhs-eqV.rhs)*om**0), s)))
