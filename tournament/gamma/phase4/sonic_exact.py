# sonic_exact.py — the EXACT sonic point and its separatrix (L'Hopital) slopes, 2nd-order series.
# Sonic point (exact): V0=1/sqrt3, N0=2/sqrt3, A0=3/2, om0=3/4.  Dson=numV=numOm=0 there.
# Metric slopes: A'=3, N'=-2/sqrt3.
#
# Fluid slopes (om',V') = L'Hopital limit of (numOm/det, numV/det). Both -> 0/0. Along the
# trajectory to 1st order: numV ~ (grad numV . Y')*x, det ~ (grad det . Y')*x, so
#   V' = (grad numV . Y')/(grad det . Y') , om' = (grad numOm . Y')/(grad det . Y'),
# with Y'=(N',A',om',V'). N',A' known -> two coupled eqns in (om',V') that reduce to a QUADRATIC.
# Two roots = the two separatrix slopes. Pick the branch matching the center-side approach (V'<0
# since V decreases outward toward sonic from the high-V small-a2 side? — determined numerically).
import sympy as sp, pickle
d = pickle.load(open("css_syms.pkl","rb"))
r = sp.Symbol('r', positive=True)
A0,N0,o0,V0 = sp.symbols('A0 N0 om0 V0')
numV = sp.sympify(d['numV']).subs(r,1)
numOm= sp.sympify(d['numOm']).subs(r,1)
det  = sp.sympify(d['det']).subs(r,1)
Axv  = sp.sympify(d['Axv']); Nxv = sp.sympify(d['Nxv'])

pt = {V0: 1/sp.sqrt(3), N0: 2/sp.sqrt(3), A0: sp.Rational(3,2), o0: sp.Rational(3,4)}
Npr = sp.simplify(Nxv.subs(pt)); Apr = sp.simplify(Axv.subs(pt))
print("exact sonic point: V0=1/sqrt3, N0=2/sqrt3, A0=3/2, om0=3/4")
print(f"  A'={Apr}, N'={Npr}")

# gradients wrt (N,A,om,V)
def gradvec(expr):
    return [sp.diff(expr, v) for v in (N0, A0, o0, V0)]
gnumV = [sp.simplify(g.subs(pt)) for g in gradvec(numV)]
gnumO = [sp.simplify(g.subs(pt)) for g in gradvec(numOm)]
gdet  = [sp.simplify(g.subs(pt)) for g in gradvec(det)]
print("  grad numV  =", gnumV)
print("  grad numOm =", gnumO)
print("  grad det   =", gdet)

op, vp = sp.symbols('op vp')                       # om', V'
Yp = [Npr, Apr, op, vp]
DnV = sum(gnumV[i]*Yp[i] for i in range(4))
DnO = sum(gnumO[i]*Yp[i] for i in range(4))
Dd  = sum(gdet[i]*Yp[i] for i in range(4))
# V' * Dd = DnV ; om' * Dd = DnO
eq1 = sp.expand(vp*Dd - DnV)
eq2 = sp.expand(op*Dd - DnO)
sol = sp.solve([eq1, eq2], [op, vp], dict=True)
print("\nseparatrix slopes (om', V'):")
roots = []
for s in sol:
    ov = sp.nsimplify(sp.simplify(s[op])); vv = sp.nsimplify(sp.simplify(s[vp]))
    ovn = complex(s[op]); vvn = complex(s[vp])
    print(f"  om'={ov} = {ovn.real:+.6f}{'' if abs(ovn.imag)<1e-12 else f'{ovn.imag:+.6f}i'} ,  "
          f"V'={vv} = {vvn.real:+.6f}{'' if abs(vvn.imag)<1e-12 else f'{vvn.imag:+.6f}i'}")
    if abs(ovn.imag) < 1e-9 and abs(vvn.imag) < 1e-9:
        roots.append((float(ovn.real), float(vvn.real)))

# ---- 2nd order series. Fields Y_i(x)=Y0 + Y1 x + Y2 x^2/2. For N,A (regular): Y'' = d/dx of slope.
# For om,V we get Y2 by matching O(x^2) of  numV - det*V' = 0  and numOm - det*om' = 0.
# Build a small-x series of numV,det,numOm along Y(x) with unknown 2nd-order coeffs and solve.
xs = sp.symbols('x')
def branch_series(om1, V1):
    N2,A2,om2,V2 = sp.symbols('N2 A2 om2 V2')
    N1 = Npr; A1 = Apr
    Nser = pt[N0] + N1*xs + N2*xs**2/2
    Aser = pt[A0] + A1*xs + A2*xs**2/2
    Oser = pt[o0] + om1*xs + om2*xs**2/2
    Vser = pt[V0] + V1*xs + V2*xs**2/2
    sub = {N0:Nser, A0:Aser, o0:Oser, V0:Vser}
    # N''=d/dx(Nxv), A''=d/dx(Axv): differentiate the RHS along the trajectory
    Nx_along = Nxv.subs(sub); Ax_along = Axv.subs(sub)
    eN = sp.series(sp.diff(Nser,xs) - Nx_along, xs, 0, 2).removeO()
    eA = sp.series(sp.diff(Aser,xs) - Ax_along, xs, 0, 2).removeO()
    # fluid: numV - det*V' = 0 and numOm - det*om' along traj, expand to O(x^2)
    nV_along = numV.subs(sub); dd_along = det.subs(sub); nO_along = numOm.subs(sub)
    eV = sp.series(nV_along - dd_along*sp.diff(Vser,xs), xs, 0, 3).removeO()
    eO = sp.series(nO_along - dd_along*sp.diff(Oser,xs), xs, 0, 3).removeO()
    eqs = []
    for e in [eN, eA, eV, eO]:
        p = sp.Poly(sp.expand(e), xs)
        for c in p.all_coeffs():
            cc = sp.simplify(c)
            if cc != 0: eqs.append(cc)
    s2 = sp.solve(eqs, [N2,A2,om2,V2], dict=True)
    return s2

print("\n2nd-order coeffs for each real separatrix branch:")
for (om1, V1) in roots:
    s2 = branch_series(sp.nsimplify(om1, rational=False), sp.nsimplify(V1, rational=False))
    print(f"  branch (om'={om1:+.5f}, V'={V1:+.5f}):")
    if s2:
        ss = s2[0]
        for k,v in ss.items():
            print(f"     {k} = {sp.simplify(v)}  = {float(v):+.6f}")
    else:
        print("     (2nd-order solve empty — check)")
