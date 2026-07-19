# nr_ec2.py — the TRUE Evans-Coleman background, by HKA's OWN construction (rflanl.tex
# sec-ss.practice): the sonic point is a ONE-PARAMETER family in V0 (eqs 4.7-4.9 give
# (N0, A0, om0) as functions of V0); shoot sonic -> center; the desired solution sits at
# the transition between the two failure modes:
#   case (1): A < 1 somewhere            (V0 + 0 side)
#   case (2): another sonic point, det=0 (V0 - 0 side)
# EC = the transition root with exactly ONE zero of V (Bogoyavlensky index 1).
#
# WHY THIS FILE EXISTS (D-032): the prior "EC" background (hka_ec.py, Stage-A golden) is
# actually the collapsing flat FRIEDMANN solution — V0 = -1/sqrt(3) = -c_s is precisely the
# Friedmann point of the sonic line (4.7-4.9 at V0=-c_s give (3/2, 2/sqrt3, 3/4)); measured
# fingerprint: V = -sqrt(1-1/A) (free-fall identity) holds to 1.9e-10 along it. The true EC
# crossing has V0 > 0 (V has one zero) and N̄'(sonic) ~= -0.35699 (the paper's footnote value).
import numpy as np, math, sys
from scipy.integrate import solve_ivp
import hka_ec as E
import nr_sonic as NS

G = 4.0/3.0; SG = math.sqrt(G - 1.0)

def sonic_values(V0, g=G):
    """(A0, N0, om0) at the sonic point, parametrized by V0 (rflanl.tex 4.7-4.9)."""
    sg = math.sqrt(g - 1.0)
    N0 = (1.0 - sg*V0)/(sg - V0)
    A0 = (g*g + 4*g - 4 + 8*sg**3*V0 - (3*g - 2)*(2 - g)*V0*V0)/(g*g*(1 - V0*V0))
    om0 = 2*sg*(sg - V0)*(1 + sg*V0)/(g*g*(1 - V0*V0))
    return A0, N0, om0

def bg_series_at(V0, order=6):
    """Sonic Taylor series at the GENERAL sonic point V0 (nr_sonic recursion, un-hardcoded).
    The order-1 solvability is QUADRATIC -> returns a list of branch solutions
    [(A,N,om,V) coeff arrays], one per real root (sorted by root value)."""
    A0, N0, om0 = sonic_values(V0)
    K = order + 2
    out = []
    N = NS.sconst(K, N0); om = NS.sconst(K, om0); V = NS.sconst(K, V0)
    A = NS._A_series(N, om, V)
    (M11, M12, M21, M22), _ = NS._fluid_MB(A, N, om, V)
    M0 = np.array([[M11[0], M12[0]], [M21[0], M22[0]]], complex)
    U, s, Vh = np.linalg.svd(M0)
    n_R = np.conj(Vh[-1]); n_L = U[:, -1]; M0pinv = np.linalg.pinv(M0)

    def build(branch_root):
        N = NS.sconst(K, N0); om = NS.sconst(K, om0); V = NS.sconst(K, V0)

        def r_of():
            A_ = NS._A_series(N, om, V)
            (m11, m12, m21, m22), (b1, b2) = NS._fluid_MB(A_, N, om, V)
            wpo = NS.sderiv_shift(om); wpv = NS.sderiv_shift(V)
            P1 = NS.smul(m11, wpo) + NS.smul(m12, wpv) - b1
            P2 = NS.smul(m21, wpo) + NS.smul(m22, wpv) - b2
            return P1, P2

        def nLr(m):
            P1, P2 = r_of(); return n_L.conj().dot(-np.array([P1[m], P2[m]]))

        p_next = None
        for m in range(0, order+1):
            if m == 1:
                om[1], V[1] = p_next[0] + branch_root*n_R[0], p_next[1] + branch_root*n_R[1]
            elif m >= 2:
                om[m], V[m] = p_next[0], p_next[1]; r0 = nLr(m)
                om[m], V[m] = p_next[0] + n_R[0], p_next[1] + n_R[1]; r1 = nLr(m)
                denom = r1 - r0
                alpha = -r0/denom if abs(denom) > 1e-11 else 0.0
                om[m], V[m] = p_next[0] + alpha*n_R[0], p_next[1] + alpha*n_R[1]
            if m < order:
                A_ = NS._A_series(N, om, V)
                rhsN = NS.smul(N, NS.sconst(K, -2.0) + A_ - (2.0 - G)*om)
                N[m+1] = rhsN[m]/(m+1)
                P1, P2 = r_of(); r_m = -np.array([P1[m], P2[m]])
                w_next = (M0pinv.dot(r_m))/(m+1)
                p_next = (w_next[0], w_next[1])
                om[m+1], V[m+1] = 0.0, 0.0
        A_ = NS._A_series(N, om, V)
        return A_[:order+1], N[:order+1], om[:order+1], V[:order+1]

    # order-1 quadratic solvability: sample nLr(1) at 3 alphas, solve the parabola
    N = NS.sconst(K, N0); om = NS.sconst(K, om0); V = NS.sconst(K, V0)
    A_ = NS._A_series(N, om, V)
    (m11, m12, m21, m22), (b1, b2) = NS._fluid_MB(A_, N, om, V)
    r1v = -np.array([(NS.smul(m11, NS.sderiv_shift(om)) + NS.smul(m12, NS.sderiv_shift(V)) - b1)[0],
                     (NS.smul(m21, NS.sderiv_shift(om)) + NS.smul(m22, NS.sderiv_shift(V)) - b2)[0]])
    # p_next at order 0:
    w1 = (M0pinv.dot(r1v))/1.0
    samp = []
    for a in (-2.0, 0.0, 2.0):
        Nt = NS.sconst(K, N0); omt = NS.sconst(K, om0); Vt = NS.sconst(K, V0)
        # N1 from the metric row:
        At = NS._A_series(Nt, omt, Vt)
        rhsN = NS.smul(Nt, NS.sconst(K, -2.0) + At - (2.0 - G)*omt)
        Nt[1] = rhsN[0]
        omt[1], Vt[1] = w1[0] + a*n_R[0], w1[1] + a*n_R[1]
        At = NS._A_series(Nt, omt, Vt)
        (m11t, m12t, m21t, m22t), (b1t, b2t) = NS._fluid_MB(At, Nt, omt, Vt)
        P1 = NS.smul(m11t, NS.sderiv_shift(omt)) + NS.smul(m12t, NS.sderiv_shift(Vt)) - b1t
        P2 = NS.smul(m21t, NS.sderiv_shift(omt)) + NS.smul(m22t, NS.sderiv_shift(Vt)) - b2t
        samp.append(n_L.conj().dot(-np.array([P1[1], P2[1]])))
    C = samp[1]; Aq = (samp[0] + samp[2] - 2*C)/8.0; Bq = (samp[2] - samp[0])/4.0
    if abs(Aq) > 1e-12:
        roots = np.roots([Aq, Bq, C])
    else:
        roots = np.array([-C/Bq])
    roots = np.sort(np.real_if_close(roots, tol=1e6))
    branches = []
    for rt in roots:
        if abs(np.imag(rt)) > 1e-8*max(1.0, abs(np.real(rt))):
            continue                        # complex root: branch not allowed (paper req. 3)
        branches.append(build(float(np.real(rt))))
    return branches

def seed_state(coeffs, t):
    A, N, om, V = coeffs
    def ev(a): return float(np.real(sum(a[k]*t**k for k in range(len(a)))))
    return np.array([ev(N), ev(om), ev(V)])

def _ev_A1(x, Y):    # case (1): A < 1
    return E.A_of(Y[0], Y[1], Y[2]) - 1.0
_ev_A1.terminal = True; _ev_A1.direction = -1
def _ev_det(x, Y):   # case (2): another sonic point
    return E.Dson(Y[0], Y[2])
_ev_det.terminal = True; _ev_det.direction = 0

def shoot_center(V0, branch=0, delta=1e-3, x_end=-14.0, order=6, rtol=1e-11, atol=1e-13):
    """Integrate sonic(V0, branch) -> center. Returns dict with:
    outcome in {'case1' (A<1), 'case2' (det=0 / Lipschitz stall), 'blowup', 'deep', 'nobranch'}
    F = NV + 1/2 at the deepest healthy point (the center fixed point's unstable direction —
    the SIGNED shooting function; EC root = its zero)."""
    brs = bg_series_at(V0, order)
    if branch >= len(brs):
        return dict(outcome='nobranch', F=None, sol=None)
    Y0 = seed_state(brs[branch], -delta)
    with np.errstate(all='ignore'):
        sol = solve_ivp(E.rhs3, [-delta, x_end], Y0, method='DOP853', rtol=rtol, atol=atol,
                        events=[_ev_A1, _ev_det], dense_output=True, max_step=0.05)
    N, om, V = sol.y[:, -1]
    A = E.A_of(N, om, V); det = E.Dson(N, V)
    F = N*V + 0.5
    if sol.status == 1 and len(sol.t_events[0]) > 0:
        out = 'case1'
    elif (sol.status == 1 and len(sol.t_events[1]) > 0) or (sol.status < 0 and abs(det) < 0.2):
        out = 'case2'
    elif sol.status < 0:
        out = 'blowup'
    else:
        out = 'deep'                      # reached x_end without a verdict
    return dict(outcome=out, F=F, sol=sol, x_stop=sol.t[-1], A=A, det=det, om=om)

def scan(v_lo=-0.45, v_hi=0.50, n=39, branch=0, verbose=True):
    outs = []
    if verbose: print(f"scan V0 in [{v_lo},{v_hi}] branch={branch}:")
    for V0 in np.linspace(v_lo, v_hi, n):
        if abs(V0 - (-1/math.sqrt(3))) < 1e-9: continue
        r = shoot_center(V0, branch)
        outs.append((V0, r))
        if verbose and r['outcome'] != 'nobranch':
            print(f"  V0={V0:+.4f}: {r['outcome']:8s} x_stop={r['x_stop']:+7.3f}  "
                  f"F=NV+1/2={r['F']:+9.4f}  A={r['A']:8.3f} det={r['det']:+8.3f} om={r['om']:.2e}")
    return outs

def bisect_ec(v_lo, v_hi, branch=0, tol=5e-14):
    """brentq-style bisection on sign(F) (the center unstable direction)."""
    from scipy.optimize import brentq
    def f(V0):
        r = shoot_center(V0, branch)
        return r['F'] if r['F'] is not None else float('nan')
    return brentq(f, v_lo, v_hi, xtol=tol, rtol=1e-14)

def validate(V0, branch=0, x_end=-14.0):
    r = shoot_center(V0, branch, x_end=x_end)
    o, sol = r['outcome'], r['sol']
    xs = np.linspace(-1e-3, max(x_end, sol.t[-1]), 4000)
    Y = sol.sol(xs)
    N, om, V = Y
    A = E.A_of(N, om, V)
    nzero = int(np.sum(np.diff(np.sign(V)) != 0))
    A0, N0, om0 = sonic_values(V0)
    Nbp = -2.0 + A0 - (2.0 - G)*om0
    fried = np.max(np.abs(V + np.sqrt(np.maximum(1 - 1/A, 0))))
    print(f"V0 = {V0:.12f}  (branch {branch}, outcome {o}, reached x={sol.t[-1]:.2f})")
    print(f"  sonic: A0={A0:.8f} N0={N0:.8f} om0={om0:.8f}")
    print(f"  Nbar'(sonic) = -2+A0-(2-g)om0 = {Nbp:+.6f}   (paper's sonic-gauge fingerprint: -0.35699)")
    print(f"  zeros of V in (center, sonic): {nzero}   (Evans-Coleman: exactly 1)")
    print(f"  max|V+sqrt(1-1/A)| (Friedmann identity; should be LARGE now): {fried:.3e}")
    print(f"  center check: A->{A[-1]:.6f} (1)  NV->{N[-1]*V[-1]:.6f} (-1/2)  om->{om[-1]:.2e} (0)")
    return sol

def nzeros(sol):
    xs = np.linspace(sol.t[0], sol.t[-1], 3000)
    V = sol.sol(xs)[2]
    return int(np.sum(np.diff(np.sign(V)) != 0))

def shoot_phys(V0, **kw):
    """Evaluate BOTH order-1 branches, return the PHYSICAL one: outcome case1/case2 with the
    deepest reach (the sorted-root 'branch' index flickers along V0 — select by behavior)."""
    best = None
    for b in (0, 1):
        r = shoot_center(V0, b, **kw)
        if r['sol'] is None: continue
        if r['outcome'] not in ('case1', 'case2'): continue
        if r['x_stop'] > -0.5: continue          # immediate-death spurious sibling
        if best is None or r['x_stop'] < best['x_stop']:
            best = r; best['branch'] = b
    return best

def ec_root(v_lo=0.110, v_hi=0.116, tol=1e-13):
    """The EC root: bisection on the physical branch's OUTCOME (the paper's criterion):
    case2 (another sonic point) = below the root, case1 (A<1) = above. Nudge-tolerant
    (at rare V0 both order-1 branches land spurious; probe tiny offsets)."""
    def side(V0):
        for dv in (0.0, 1e-9, -1e-9, 3e-9, -3e-9, 1e-8, -1e-8, 1e-7, -1e-7):
            r = shoot_phys(V0 + dv)
            if r is not None:
                return +1 if r['outcome'] == 'case2' else -1
            # near the root the physical run goes deep ON the EC orbit (A(x_end)~1) and then
            # either reaches x_end ('deep') or stalls in the peel ('blowup') — classify by the
            # peel direction sign(F): + = case2 side (om->0, V->+), - = case1 side.
            # Exclude the spurious N->0 attractor (A(x_end) >> 1, F = +0.5 exactly).
            for b in (0, 1):
                rr = shoot_center(V0 + dv, b)
                if rr['sol'] is None or rr['outcome'] not in ('deep', 'blowup'): continue
                if rr['A'] < 2.0 and abs(rr['F'] - 0.5) > 1e-3:
                    return +1 if rr['F'] > 0 else -1
        return 0
    s_lo, s_hi = side(v_lo), side(v_hi)
    assert s_lo == +1 and s_hi == -1, f"not bracketed: side({v_lo})={s_lo}, side({v_hi})={s_hi}"
    while v_hi - v_lo > tol:
        vm = 0.5*(v_lo + v_hi)
        s = side(vm)
        if s == 0:
            vm = v_lo + 0.37*(v_hi - v_lo)      # deterministic off-center retry
            s = side(vm)
            if s == 0:
                raise RuntimeError(f"unclassifiable region near {vm}")
        if s == +1: v_lo = vm
        else: v_hi = vm
    return 0.5*(v_lo + v_hi)

def tail_coeffs(V0, x_fit=None):
    """Extract the center-asymptotic coefficients (eq:ss-asymp1) from the deepest clean tail:
    om ~ om_inf e^{2x}, N ~ N_inf e^{-x}, V ~ V_inf e^{x}; constraints A_inf = 2 om_inf/3,
    N_inf V_inf = -1/2 validate the fit."""
    r = shoot_phys(V0)
    sol = r['sol']
    if x_fit is None: x_fit = sol.t[-1] + 0.7      # a bit above the stall
    N, om, V = sol.sol(x_fit)
    A = E.A_of(N, om, V)
    om_inf = om*math.exp(-2*x_fit); N_inf = N*math.exp(x_fit); V_inf = V*math.exp(-x_fit)
    A_inf = (A - 1.0)*math.exp(-2*x_fit)
    return dict(x_fit=x_fit, om_inf=om_inf, N_inf=N_inf, V_inf=V_inf, A_inf=A_inf,
                c1=A_inf/(2*om_inf/3.0), c2=N_inf*V_inf/(-0.5), r=r)

_EC = None
def build_ec(V0=None, x0=-16.0, verbose=True):
    """THE TRUE EC BACKGROUND: root-find V0 (if not given), extract center data from the tail,
    re-launch center -> sonic (STABLE direction) for a dense global profile. Returns dict with
    dense sol, xs (the Dson=0 stop), and bg_state(x) in hka_beta4 format."""
    global _EC
    if _EC is not None and V0 is None: return _EC
    if V0 is None: V0 = ec_root()
    tc = tail_coeffs(V0)
    z0 = math.exp(x0)
    Y0 = E.center_state3(tc['N_inf'], tc['om_inf'], z0)
    with np.errstate(all='ignore'):
        sol = solve_ivp(E.rhs3, [x0, 1.0], Y0, method='DOP853', rtol=1e-12, atol=1e-14,
                        events=E._event_Dson, dense_output=True, max_step=0.02)
    xs = sol.t_events[0][0] if len(sol.t_events[0]) else sol.t[-1]
    Ns, oms, Vs = sol.sol(xs - 1e-9)
    As = E.A_of(Ns, oms, Vs)
    A0, N0, om0 = sonic_values(V0)
    if verbose:
        print(f"[build_ec] V0* = {V0:.12f}")
        print(f"  center data from tail: om_inf(oi_EC) = {tc['om_inf']:.8f}  N_inf = {tc['N_inf']:.8f}  "
              f"(constraint checks c1={tc['c1']:.6f} c2={tc['c2']:.6f}; both should be ~1)")
        print(f"  re-launched center->sonic: hits Dson=0 at x_s = {xs:.8f}")
        print(f"  sonic state reached: A={As:.6f} N={Ns:.6f} om={oms:.6f} V={Vs:.6f}")
        print(f"  sonic target (4.7-4.9): A0={A0:.6f} N0={N0:.6f} om0={om0:.6f} V0={V0:.6f}")
        print(f"  Nbar'(sonic) = {-2.0 + A0 - (2.0-G)*om0:+.6f}   (paper: -0.35699)")
    import hka_pert_core as PC
    def bg_state(x):
        N_, om_, V_ = sol.sol(x)
        return PC.bg_fields(N_, om_, V_)
    _EC = dict(V0=V0, xs=xs, sol=sol, bg_state=bg_state, tc=tc)
    return _EC

def fine(v_lo=0.06, v_hi=0.16, step=0.002):
    """Fine scan, BOTH branches, with V-zero counts and reach depth — locate the EC transition
    (reach diverges, zeros=1, fingerprint -> -0.35699)."""
    print(f"{'V0':>8} | {'b0 outcome':>10} {'x_stop':>7} {'F':>8} {'nz':>3} | {'b1 outcome':>10} {'x_stop':>7} {'F':>8} {'nz':>3}")
    for V0 in np.arange(v_lo, v_hi + step/2, step):
        row = f"{V0:8.4f} |"
        for b in (0, 1):
            r = shoot_center(V0, b)
            if r['sol'] is None:
                row += f" {'nobranch':>10} {'':>7} {'':>8} {'':>3} |"
            else:
                nz = nzeros(r['sol'])
                row += f" {r['outcome']:>10} {r['x_stop']:7.2f} {r['F']:8.3f} {nz:3d} |"
        print(row, flush=True)

if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "find"
    if mode == "fine":
        a = float(sys.argv[2]) if len(sys.argv) > 2 else 0.06
        b = float(sys.argv[3]) if len(sys.argv) > 3 else 0.16
        st = float(sys.argv[4]) if len(sys.argv) > 4 else 0.002
        fine(a, b, st)
    elif mode == "scan":
        b = int(sys.argv[2]) if len(sys.argv) > 2 else 0
        scan(branch=b)
    elif mode == "friedcheck":
        # control: V0 = -1/sqrt3 must reproduce the Friedmann values (3/2, 2/sqrt3, 3/4)
        A0, N0, om0 = sonic_values(-1/math.sqrt(3))
        print(f"V0=-1/sqrt3: A0={A0:.9f} (1.5)  N0={N0:.9f} ({2/math.sqrt(3):.9f})  om0={om0:.9f} (0.75)")
    else:
        # find: coarse scan branch 0+1, bisect every sign change of F=NV+1/2, validate
        for b in (0, 1):
            outs = scan(branch=b)
            for i in range(1, len(outs)):
                (v1, r1), (v2, r2) = outs[i-1], outs[i]
                if r1['F'] is None or r2['F'] is None: continue
                if np.isfinite(r1['F']) and np.isfinite(r2['F']) and r1['F']*r2['F'] < 0:
                    try:
                        v = bisect_ec(v1, v2, branch=b)
                    except Exception as e:
                        print(f"  [bisect {v1:+.3f}..{v2:+.3f} b{b} failed: {e}]"); continue
                    print(f"\n=== F-root (branch {b}) at V0 = {v:.12f} ===")
                    validate(v, branch=b)
