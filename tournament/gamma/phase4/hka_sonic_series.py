# hka_sonic_series.py — local Taylor series of the HKA (4.1) background about the sonic point x=0.
#
# HKA give sonic DATA (4.7-4.9) but OMIT the higher coefficients ("we do not give explicit
# expressions"). We generate them by expanding (4.1a,b,d,e) about x=0.
#
# Method: write A(x)=A0+A1 x+A2 x^2+..., same for N,om,V. The metric eqs (4.1a,b) are explicit
# first-order ODEs -> give A1,N1,... directly once om,V known. The fluid pair (4.1d,e) is the
# 2x2 system M(x).(om_x,V_x) = rhs(x); at x=0 M is SINGULAR (sonic). Expanding to next order,
# the O(x^0) fluid equations give a QUADRATIC for (om1,V1) [the L'Hopital slopes]; higher orders
# are linear. We pick the branch whose trajectory flows to the regular center (verified downstream).
#
# Everything symbolic in V0 first, then substitute the numeric V0 from the Stage-A shoot.
# Eq numbers per HKA_beta_equations.md.

import sympy as sp

g = sp.Rational(4, 3)
cs = sp.sqrt(g - 1)
x = sp.symbols('x', real=True)

def sonic_data(V0):
    """HKA (4.7)-(4.9): sonic point data as exact functions of V0 (rational/radical)."""
    N0 = (1 - cs*V0)/(cs - V0)
    A0 = (g**2 + 4*g - 4 + 8*(g-1)**sp.Rational(3,2)*V0 - (3*g-2)*(2-g)*V0**2)/(g**2*(1-V0**2))
    om0 = 2*cs*(cs - V0)*(1 + cs*V0)/(g**2*(1-V0**2))
    return A0, N0, om0

def build_series(V0val, order=4):
    """Return dict of Taylor coeffs {A:[A0..], N:[...], om:[...], V:[...]} at x=0 for given numeric V0.
    Uses exact rational V0 if V0val is a sympy Rational; else works with high-precision floats symbolically.
    Branch: pick the (om1,V1) root that is real and gives the center-directed flow (V decreasing toward
    the center means, on the ingoing branch V0<0, the trajectory has |V| growing away from sonic going to
    +x and V->0 going to -x; concretely we pick the root matching the sign the two-sided shoot needs and
    verify by integration)."""
    A0v, N0v, om0v = sonic_data(V0val)
    A0v, N0v, om0v = sp.nsimplify(A0v), sp.nsimplify(N0v), sp.nsimplify(om0v)

    # unknown coeffs
    ncA = [A0v] + list(sp.symbols(f'A1:{order+1}', real=True))
    ncN = [N0v] + list(sp.symbols(f'N1:{order+1}', real=True))
    ncO = [om0v] + list(sp.symbols(f'O1:{order+1}', real=True))
    ncV = [sp.nsimplify(V0val)] + list(sp.symbols(f'V1:{order+1}', real=True))

    Aser = sum(ncA[i]*x**i for i in range(order+1))
    Nser = sum(ncN[i]*x**i for i in range(order+1))
    Oser = sum(ncO[i]*x**i for i in range(order+1))
    Vser = sum(ncV[i]*x**i for i in range(order+1))

    Axs = sp.diff(Aser, x); Nxs = sp.diff(Nser, x); Oxs = sp.diff(Oser, x); Vxs = sp.diff(Vser, x)

    # Residuals of (4.1a,b,d,e). Multiply through by (1-V^2) and om to clear denominators.
    one_m_V2 = (1 - Vser**2)
    # (4.1a): A_x - A[1 - A + 2 om(1+(g-1)V^2)/(1-V^2)] = 0  -> * (1-V^2)
    R_a = sp.expand((Axs - Aser*(1 - Aser))*one_m_V2 - Aser*2*Oser*(1+(g-1)*Vser**2))
    # (4.1b): N_x - N[-2 + A - (2-g)om] = 0
    R_b = sp.expand(Nxs - Nser*(-2 + Aser - (2-g)*Oser))
    # (4.1d): (1+NV) om_x/om + g(N+V)V_x/(1-V^2) - RHS_d = 0  -> * om*(1-V^2)
    RHS_d = (3*(2-g)/2)*Nser*Vser - ((2+g)/2)*Aser*Nser*Vser + (2-g)*Nser*Vser*Oser
    R_d = sp.expand((1+Nser*Vser)*Oxs*one_m_V2 + g*(Nser+Vser)*Vxs*Oser - RHS_d*Oser*one_m_V2)
    # (4.1e): (g-1)(N+V) om_x/om + g(1+NV)V_x/(1-V^2) - RHS_e = 0 -> * om*(1-V^2)
    RHS_e = (2-g)*(g-1)*Nser*Oser + ((7*g-6)/2)*Nser + ((2-3*g)/2)*Aser*Nser
    R_e = sp.expand((g-1)*(Nser+Vser)*Oxs*one_m_V2 + g*(1+Nser*Vser)*Vxs*Oser - RHS_e*Oser*one_m_V2)

    unknowns = ncA[1:] + ncN[1:] + ncO[1:] + ncV[1:]
    sol = {}
    # Solve order by order in x. At each power p, collect coeff of x^p from each residual.
    # Metric eqs R_a,R_b give A_{p+1},N_{p+1} linearly. Fluid eqs at lowest order give the quadratic
    # for (O1,V1); at higher order linear for (O_{p+1},V_{p+1}).
    for p in range(0, order):
        eqs = []
        for R in (R_a, R_b, R_d, R_e):
            Rp = R.subs(sol)
            c = Rp.coeff(x, p)
            if c != 0:
                eqs.append(sp.expand(c))
        # unknowns appearing at this order: A_{p+1},N_{p+1},O_{p+1},V_{p+1}
        newun = [ncA[p+1], ncN[p+1], ncO[p+1], ncV[p+1]]
        newun = [u for u in newun if u in unknowns]
        s = sp.solve(eqs, newun, dict=True)
        if not s:
            raise RuntimeError(f"no solution at order {p}")
        if len(s) > 1:
            # quadratic branch (order 0): pick the real root. Report both, choose by V1 sign.
            reals = []
            for cand in s:
                vals = {k: sp.nsimplify(v) for k, v in cand.items()}
                # numeric check real
                ok = all(abs(sp.im(sp.N(v))) < 1e-9 for v in vals.values())
                if ok:
                    reals.append(vals)
            if not reals:
                raise RuntimeError(f"no real branch at order {p}: {s}")
            # choose branch: the center-directed one. Heuristic: on ingoing branch (V0<0),
            # the EC solution has V passing through the sonic point with V1 of a definite sign.
            # Try both; the downstream shoot verifies. Here pick the one with the LARGER |slope|
            # is wrong in general — instead we keep BOTH and let build_series_both expose them.
            chosen = reals  # return list; caller picks
            return {'branches': chosen, 'A0': A0v, 'N0': N0v, 'om0': om0v, 'V0': sp.nsimplify(V0val),
                    'ncA': ncA, 'ncN': ncN, 'ncO': ncO, 'ncV': ncV,
                    'R': (R_a, R_b, R_d, R_e), 'order': order, 'x': x, 'p_quad': p, 'sol_before': dict(sol)}
        sol.update({k: sp.nsimplify(v) for k, v in s[0].items()})
    # no quadratic encountered (shouldn't happen); return single
    return {'branches': [sol], 'A0': A0v, 'N0': N0v, 'om0': om0v, 'V0': sp.nsimplify(V0val)}


def finish_branch(info, branch_sol):
    """Given the quadratic-branch choice at order p_quad, continue solving higher orders linearly."""
    R_a, R_b, R_d, R_e = info['R']
    ncA, ncN, ncO, ncV = info['ncA'], info['ncN'], info['ncO'], info['ncV']
    order = info['order']; p0 = info['p_quad']
    sol = dict(info['sol_before']); sol.update(branch_sol)
    for p in range(p0+1, order):
        eqs = []
        for R in (R_a, R_b, R_d, R_e):
            Rp = R.subs(sol); c = Rp.coeff(x, p)
            if c != 0:
                eqs.append(sp.expand(c))
        newun = [ncA[p+1], ncN[p+1], ncO[p+1], ncV[p+1]]
        newun = [u for u in newun if u not in sol]
        s = sp.solve(eqs, newun, dict=True)
        if not s:
            raise RuntimeError(f"no solution at order {p}")
        sol.update({k: sp.nsimplify(v) for k, v in s[0].items()})
    # assemble coefficient lists
    def coefs(nc):
        return [nc[0]] + [sol.get(nc[i], sp.Integer(0)) for i in range(1, order+1)]
    return {'A': coefs(ncA), 'N': coefs(ncN), 'om': coefs(ncO), 'V': coefs(ncV)}


if __name__ == "__main__":
    # Explore branches at a trial V0 to understand the quadratic.
    for v0 in [sp.Rational(-1,4)]:
        info = build_series(v0, order=4)
        print(f"V0={float(v0):+.4f}  A0={float(info['A0']):.5f} N0={float(info['N0']):.5f} om0={float(info['om0']):.5f}")
        print(f"quadratic branch at order p={info['p_quad']}, {len(info['branches'])} real root(s):")
        for bi, br in enumerate(info['branches']):
            v1 = br.get(info['ncV'][1]); o1 = br.get(info['ncO'][1])
            a1 = br.get(info['ncA'][1]); n1 = br.get(info['ncN'][1])
            print(f"  branch {bi}: A1={float(a1):+.5f} N1={float(n1):+.5f} om1={float(o1):+.5f} V1={float(v1):+.5f}")
            full = finish_branch(info, br)
            print("     A:", [f"{float(c):+.4g}" for c in full['A']])
            print("     N:", [f"{float(c):+.4g}" for c in full['N']])
            print("     om:", [f"{float(c):+.4g}" for c in full['om']])
            print("     V:", [f"{float(c):+.4g}" for c in full['V']])
