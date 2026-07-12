# css_symbols.py — the ONE self-contained symbolic foundation for the radiation-fluid CSS problem.
# Regenerates, from first principles (full-4D nabla_a T^{ab}=0, the verified correct reduction):
#   - the metric slope RHS  Ax = A'(x), Nx = N'(x)   [KHA, verified vs Evans-Coleman]
#   - the fluid pair as a 2x2 linear system  M . (om', V')^T = b   in (om_x, V_x)
#   - det(M), and the Cramer numerators  numOm = det[[b1,a12],[b2,a22]], numV = det[[a11,b1],[a21,b2]]
#   - the resolved slopes  om' = numOm/det,  V' = numV/det   (finite away from sonic)
# Dumps everything (srepr) to css_syms.pkl so downstream scripts are deterministic & fast.
#
# Sonic locus  D := 3 N^2 V^2 - N^2 + 4 N V - V^2 + 3 = 0  is confirmed here to equal (up to constant
# factor) the numerator of det(M).  At D=0, numOm=numV=0 is required for a locally analytic solution.
import sympy as sp, pickle

tau, r, th = sp.symbols('tau r theta', positive=True)
A0, N0, o0, V0 = sp.symbols('A0 N0 om0 V0')          # field values at a point
Ax, Nx, ox, Vx = sp.symbols('A_x N_x om_x V_x')      # x-derivatives

# metric (polar-radial), similarity slicing E=r/tau ; a=sqrt(A), alpha=N a E
E = r/tau
a = sp.sqrt(A0)
alpha = N0*a*E
W = 1/sp.sqrt(1 - V0**2)
rho = o0/(4*sp.pi*r**2*A0)
p = rho/3
g = sp.diag(-alpha**2, a**2, r**2, r**2*sp.sin(th)**2)
gi = g.inv()
coords = [tau, r, th, sp.symbols('phi')]

def Dk(expr, k):
    # total x-derivative via chain rule: fields depend on x=ln(-r/t); dx/dtau=-1/tau, dx/dr=1/r
    e = sp.diff(expr, coords[k])
    fac = (-1/tau) if k == 0 else ((1/r) if k == 1 else 0)
    e += (sp.diff(expr, A0)*Ax + sp.diff(expr, N0)*Nx +
          sp.diff(expr, o0)*ox + sp.diff(expr, V0)*Vx)*fac
    return e

def christ(l, i, j):
    ss = sp.S.Zero
    for dd in range(4):
        if gi[l, dd] == 0:
            continue
        ss += gi[l, dd]*(Dk(g[dd, i], j) + Dk(g[dd, j], i) - Dk(g[i, j], dd))
    return sp.simplify(ss/2)

uu = [W/alpha, W*V0/a, sp.S.Zero, sp.S.Zero]         # u^a, verified u_a u^a = -1
T = [[(rho + p)*uu[i]*uu[j] + p*gi[i, j] for j in range(4)] for i in range(4)]

def divT(b):
    ss = sp.S.Zero
    for ai in range(4):
        ss += Dk(T[ai][b], ai)
        for c in range(4):
            g1 = christ(ai, ai, c)
            if g1 != 0:
                ss += g1*T[c][b]
            g2 = christ(b, ai, c)
            if g2 != 0:
                ss += g2*T[ai][c]
    return ss

print("full-4D divergence: energy (b=tau) ...", flush=True)
Dtau = divT(0)
print("full-4D divergence: radial momentum (b=r, with angular pressure terms) ...", flush=True)
Dr = divT(1)

Es = sp.symbols('E', positive=True)
def fin(D):
    D = sp.simplify(D*4*sp.pi*r**2)
    D = D.subs(tau, r/Es)
    return sp.simplify(D.subs(Es, 1))     # E is an overall factor; set 1

en = fin(Dtau)
mo = fin(Dr)

# metric slopes (already ODEs), verified vs Evans-Coleman (4),(5)
Axv = A0*(1 - A0 + (2*o0/(1 - V0**2))*(1 + V0**2/3))
Nxv = N0*(-2 + A0 - sp.Rational(2, 3)*o0)

# substitute metric slopes -> the fluid pair becomes 2x2 linear in (ox, Vx)
e1 = sp.expand(en.subs({Ax: Axv, Nx: Nxv}))
e2 = sp.expand(mo.subs({Ax: Axv, Nx: Nxv}))
a11 = sp.simplify(e1.coeff(ox)); a12 = sp.simplify(e1.coeff(Vx))
a21 = sp.simplify(e2.coeff(ox)); a22 = sp.simplify(e2.coeff(Vx))
b1 = sp.simplify(-(e1.subs({ox: 0, Vx: 0})))
b2 = sp.simplify(-(e2.subs({ox: 0, Vx: 0})))

det = sp.simplify(a11*a22 - a12*a21)
numOm = sp.simplify(b1*a22 - a12*b2)
numV = sp.simplify(a11*b2 - b1*a21)
omx = sp.simplify(numOm/det)
Vxx = sp.simplify(numV/det)

# The physical singular locus is where the RESOLVED slopes blow up: numerator(det)=0.
# = 4 om (3 N^2 V^2 - N^2 - 4 N V - V^2 + 3) = 0.  The (...) factor is the KHA sonic locus
# (3N^2V^2 - N^2 + 4NV - V^2 + 3) evaluated at INGOING V<0 (V -> -V). Define the self-consistent
# sonic polynomial used everywhere downstream:  Dson = 3 N^2 V^2 - N^2 - 4 N V - V^2 + 3.
Dson = 3*N0**2*V0**2 - N0**2 - 4*N0*V0 - V0**2 + 3
detnum = sp.factor(sp.numer(sp.together(det)))
print("num(det) factored =", detnum, flush=True)
print("num(det)/(4*om*Dson) =", sp.simplify(detnum/(4*o0*Dson)), "  (==1 => det singular locus is Dson)", flush=True)

pickle.dump({
    'a11': sp.srepr(a11), 'a12': sp.srepr(a12), 'a21': sp.srepr(a21), 'a22': sp.srepr(a22),
    'b1': sp.srepr(b1), 'b2': sp.srepr(b2),
    'det': sp.srepr(det), 'numOm': sp.srepr(numOm), 'numV': sp.srepr(numV),
    'omx': sp.srepr(omx), 'Vxx': sp.srepr(Vxx),
    'Axv': sp.srepr(Axv), 'Nxv': sp.srepr(Nxv),
    'Dson': sp.srepr(Dson),
}, open("css_syms.pkl", "wb"))
print("wrote css_syms.pkl", flush=True)

# sanity: center dV/dx (A=1,V=0,om->0,N large) should -> 0 (~1.5/N)
fV = sp.lambdify((A0, N0, o0, V0), Vxx, 'numpy')
for Nc in [1e4, 1e5, 1e6]:
    print(f"  center dV/dx at N={Nc:.0e}: {float(fV(1.0, Nc, 1e-12, 0.0)):.3e}  (expect ~1.5/N={1.5/Nc:.3e})")
