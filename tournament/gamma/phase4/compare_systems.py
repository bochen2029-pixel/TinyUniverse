# compare_systems.py — is the re-derived css_symbols fluid system equal to the KHA-printed system
# under the velocity flip V -> -V (a candidate symmetry)? If yes, spectra match and the sign is a
# harmless relabeling. If no, the re-derivation changed the physics -> the source of the wrong spectrum.
import sympy as sp, pickle
N,A,w,V = sp.symbols('N A w V', real=True)
r=sp.Symbol('r',positive=True)

# --- KHA-printed resolved slopes (from kha_verbatim logic) ---
Ap = A*(1 - A + (2*w/(1-V**2))*(1+V**2/3))
Np = N*(-2 + A - sp.Rational(2,3)*w)
g,u = sp.symbols('g u')
f1 = (1+N*V)*g + 4*(N+V)*u/(3*(1-V**2)) - N*V*Ap/(3*A) + 4*V*Np/3 + 2*N*(1+4*w/(9*(1-V**2)))
f2 = (4*V+N+3*N*V**2)*g + 4*(1+V**2+2*N*V)*u/(1-V**2) + N*(1-V**2)*Ap/A + 4*(1+V**2)*Np + 2*N*(1+3*V**2)
sol=sp.solve([f1,f2],[g,u],dict=True)[0]
KHA_wp = sp.simplify(w*sol[g]); KHA_Vp = sp.simplify(sol[u])

# --- css_symbols re-derived resolved slopes ---
d=pickle.load(open("css_syms.pkl","rb"))
A0,N0,o0,V0=sp.symbols('A0 N0 om0 V0')
CSS_wp = sp.sympify(d['omx']).subs(r,1)
CSS_Vp = sp.sympify(d['Vxx']).subs(r,1)
# map css symbols (A0,N0,om0,V0) -> (A,N,w,V)
sub={A0:A,N0:N,o0:w,V0:V}
CSS_wp=CSS_wp.subs(sub); CSS_Vp=CSS_Vp.subs(sub)

print("Compare V' slopes:")
# under V->-V, dV/dx -> -dV/dx (since V flips and x same). So test: CSS_Vp(V) ?= -KHA_Vp(-V)
flip = {V:-V}
lhs = sp.simplify(CSS_Vp)
rhs = sp.simplify(-KHA_Vp.subs(flip))
diff = sp.simplify(lhs-rhs)
print("  CSS_Vp - (-KHA_Vp|_{V->-V}) =", diff, "   (0 => V-flip symmetry holds for V')")
# w' is even under V->-V (w scalar): CSS_wp(V) ?= KHA_wp(-V)
lhsw=sp.simplify(CSS_wp); rhsw=sp.simplify(KHA_wp.subs(flip))
diffw=sp.simplify(lhsw-rhsw)
print("  CSS_wp - (KHA_wp|_{V->-V})   =", diffw, "   (0 => V-flip symmetry holds for w')")

# also directly: are they equal WITHOUT flip? (i.e., did css just reproduce KHA?)
print("\nDirect (no flip):")
print("  CSS_Vp - KHA_Vp =", sp.simplify(CSS_Vp-KHA_Vp))
print("  CSS_wp - KHA_wp =", sp.simplify(CSS_wp-KHA_wp))
