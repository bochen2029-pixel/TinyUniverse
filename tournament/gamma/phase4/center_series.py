# center_series.py — DERIVE (not guess) the regular-center local series, then trace the ingoing
# trajectory of the correct system to see exactly how it meets the singular locus Dson=0.
#
# Center at x->-inf. Let z=e^x->0. Regular center: A->1, V->0, and spacetime regular. Ansatz
#   N = nc/z * (1 + n1 z + n2 z^2 + ...)     (N ~ e^{-x} leading; N'/N = -2+A-2om/3 -> -1 as A->1,om->0? )
# Wait: N'/N = -2 + A - 2om/3. At center A->1,om->0 => N'/N -> -1 => N ~ e^{-x}=1/z. Good: N=nc/z(1+...).
#   A = 1 + a2 z^2 + ...     (even; A'->0)
#   om = w2 z^2 + ...        (om ~ 4pi r^2 a^2 rho, finite central density => om ~ r^2 ~ z^2)
#   V = v1 z + v3 z^3 + ...  (ingoing, odd)
# Plug into the 4 ODEs Y'=F(Y) with d/dx = z d/dz, match lowest orders to fix (n1,a2 vs w2, v1).
import sympy as sp
z = sp.symbols('z', positive=True)
nc, a2, w2, v1, n1, v3, a4, w4 = sp.symbols('nc a2 w2 v1 n1 v3 a4 w4', real=True)

# series
N = nc/z*(1 + n1*z)            # keep to needed order
A = 1 + a2*z**2
om = w2*z**2
V = v1*z + v3*z**3

# metric slopes as ODEs: A'/A, N'/N ; fluid slopes we take from the resolved expressions but near
# center it's cleaner to use the two fluid PDEs directly. Use metric eqns + the two constraints that
# the resolved om',V' satisfy. Simplest: impose the 4 slope eqns to leading order.
Ax_rhs = A*(1 - A + (2*om/(1-V**2))*(1 + V**2/3))
Nx_rhs = N*(-2 + A - sp.Rational(2,3)*om)
# d/dx = z d/dz
def ddx(f): return z*sp.diff(f, z)
eqN = sp.series((ddx(N) - Nx_rhs), z, 0, 3).removeO()
eqA = sp.series((ddx(A) - Ax_rhs), z, 0, 5).removeO()
print("N-eq (orders of z):", sp.expand(eqN))
print("A-eq (orders of z):", sp.expand(eqA))

# For the fluid: use the resolved slopes from the pickle (exact), expand in z.
import pickle
d = pickle.load(open("css_syms.pkl","rb"))
r = sp.Symbol('r', positive=True)
A0,N0,o0,V0 = sp.symbols('A0 N0 om0 V0')
omx = sp.sympify(d['omx']).subs(r,1)
Vxx = sp.sympify(d['Vxx']).subs(r,1)
sub = {A0:A, N0:N, o0:om, V0:V}
omx_s = omx.subs(sub); Vxx_s = Vxx.subs(sub)
eqom = sp.series(ddx(om) - omx_s, z, 0, 4).removeO()
eqV  = sp.series(ddx(V) - Vxx_s, z, 0, 4).removeO()
print("om-eq:", sp.expand(sp.simplify(eqom)))
print("V-eq :", sp.expand(sp.simplify(eqV)))

# Collect and solve order-by-order for (n1,a2,w2,v1,v3) leaving nc, a2 free? Let's see which are fixed.
print("\n--- solving leading balances ---")
alleqs = []
for e in [eqN, eqA, eqom, eqV]:
    p = sp.Poly(sp.expand(e), z)
    for c in p.all_coeffs():
        cc = sp.simplify(c)
        if cc != 0:
            alleqs.append(cc)
alleqs = list(sp.ordered(set(alleqs)))
for e in alleqs:
    print("  =0:", e)
sol = sp.solve(alleqs, [n1, w2, v1, v3], dict=True)
print("solution(s):", sol)
