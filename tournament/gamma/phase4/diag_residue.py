# diag_residue.py — understand the Fuchsian structure at the sonic point.
# Build the residue matrix R (L ~ R/(x_s - x) near sonic), get its eigenvalues (indicial exponents),
# and check whether the mu=1-2k branch is present with the right eigenvector.
import numpy as np, math, sympy as sp, pickle
import css_core as C

sp.init_printing()
S3 = math.sqrt(3.0)

# --- Build L(x) as a symbolic matrix along the EXACT background series about x=0 (sonic point) ---
_d = pickle.load(open("Lmat.pkl","rb"))
A0,N0,o0,V0 = sp.symbols('A0 N0 om0 V0')
Nx,Ax,ox,Vx = sp.symbols('N_x A_x om_x V_x')
k = sp.symbols('k')
L = sp.sympify(_d['Lmat'])
x = sp.symbols('x')

# EXACT background series (branch 1, verified by sonic_exact.py):
Ns = 2/sp.sqrt(3) + (-2/sp.sqrt(3))*x + sp.Rational(1,2)*(2/sp.sqrt(3))*x**2
As = sp.Rational(3,2) + 3*x + sp.Rational(1,2)*33*x**2
Os = sp.Rational(3,4) + sp.Rational(9,2)*x + sp.Rational(1,2)*sp.Rational(99,2)*x**2
Vs = 1/sp.sqrt(3) + (2/sp.sqrt(3))*x + sp.Rational(1,2)*(10*sp.sqrt(3)/3)*x**2
sub = {N0:Ns,A0:As,o0:Os,V0:Vs, Nx:sp.diff(Ns,x),Ax:sp.diff(As,x),ox:sp.diff(Os,x),Vx:sp.diff(Vs,x)}
Lx = L.subs(sub)

# The claim: L(x) has a simple pole at x=0 -> L ~ R/x + L_reg. But note the sonic point is at x=0 in
# THIS series parameterization. Dson(N(x),V(x)) ~ Dson'(0)*x for small x. Let's find the pole order.
print("=== Pole structure of L(x) at the sonic point x=0 ===")
Dson_series = sp.series(C.Dson(Ns,Vs), x, 0, 3).removeO() if False else None
# use symbolic Dson
Dsym = 3*Ns*Ns*Vs*Vs - Ns*Ns - 4*Ns*Vs - Vs*Vs + 3
Ds = sp.series(Dsym, x, 0, 3)
print("Dson(x) series about sonic:", sp.simplify(Ds))

# Build R = residue = lim x*L(x) as x->0 (Laurent leading), entrywise
R = sp.zeros(4,4)
poleorder = 0
for i in range(4):
    for j in range(4):
        lij = Lx[i,j]
        s = sp.series(lij, x, 0, 1)   # need to see if it has 1/x term
        # get the coefficient of 1/x
        ser = sp.series(lij*x, x, 0, 1).removeO()
        R[i,j] = sp.nsimplify(sp.simplify(ser.subs(x,0))) if ser != 0 else sp.Integer(0)
R = sp.simplify(R)
print("\n=== Residue matrix R (entrywise lim x*L) ===")
sp.pprint(R)

# Indicial equation: for hp ~ x^mu * v (with x measuring distance PAST... careful with sign),
# the ODE hp' = L hp = (R/x + ...) hp gives, to leading order, mu x^{mu-1} v = R x^{mu-1} v,
# so R v = mu v : indicial exponents = eigenvalues of R.
print("\n=== Eigenvalues of R (indicial exponents mu) ===")
ev = R.eigenvals()
for e,mult in ev.items():
    print("  mu =", sp.simplify(e), "  (mult", mult, ")")

print("\n=== Symbolic eigenvalues as function of k check ===")
# Expect exponents {0,0,0, 1-2k} OR {0, ...}. Print charpoly.
lam = sp.symbols('lam')
cp = sp.simplify(R.charpoly(lam).as_expr())
print("charpoly(lam) =", sp.factor(cp))

# The non-analytic exponent should be mu = 1-2k. Get its RIGHT eigenvector v_na (R v = mu v)
# and the LEFT eigenvector wL (wL R = mu wL) i.e. right-null of (R - mu I)^T.
mu_na = 1 - 2*k
print("\n=== Non-analytic eigenvector for mu = 1-2k ===")
Mr = R - mu_na*sp.eye(4)
right_ns = Mr.nullspace()
print("R - (1-2k)I  right-nullspace dim =", len(right_ns))
if right_ns:
    v_na = sp.simplify(right_ns[0])
    print("v_na (right eigenvector, RAW):")
    sp.pprint(v_na.T)
    # clear denominators to see pole structure
    print("v_na denominators per component:")
    for c in range(4):
        print("   comp",c,":", sp.denom(sp.together(v_na[c])))

left_ns = (Mr.T).nullspace()
print("\nwL (left eigenvector, RAW = right-null of (R-mu I)^T):")
if left_ns:
    wL = sp.simplify(left_ns[0].T)
    sp.pprint(wL)
    print("wL denominators per component:")
    for c in range(4):
        print("   comp",c,":", sp.denom(sp.together(wL[c])))
    # POLE-NORMALIZED left vector: multiply by (2k-3) to kill the k=3/2 pole seen in stageB_v2
    wL_norm = sp.simplify(sp.expand((2*k-3)*wL))
    print("\n(2k-3)*wL (pole-normalized):")
    sp.pprint(wL_norm)

# also: analytic left duals (mu=0) — we need wL for the NON-analytic direction to project it out.
# Verify the biorthogonality: wL . v_na should be nonzero; wL . (analytic right vecs) = 0.
print("\n=== Biorthogonality sanity (numeric at k=2.81) ===")
kv = 2.81055255
Rn = np.array(R.subs(k,kv)).astype(np.complex128)
evR, VR = np.linalg.eig(Rn)
print("numeric eig(R) at k=2.81:", np.round(evR,6))
print("expected non-analytic exponent 1-2k =", 1-2*kv)
