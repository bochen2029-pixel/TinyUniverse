# nr_bgcheck.py — VALIDATE the exact background sonic Taylor series (hka_sonic_series) against the
# accurate integrated Evans-Coleman background (hka_ec). This is step-1's foundation: if the exact
# series reproduces the true background near the sonic point to ~1e-10, it replaces the inaccurate
# one-sided polyfit (hka_pert_sonic.bg_series_near_sonic) that poisoned the Laurent/Frobenius.
import numpy as np, sympy as sp, math
import hka_ec as E
import hka_sonic_series as SS
from scipy.integrate import solve_ivp

def main(order=6):
    V0f = -1.0/math.sqrt(3.0)
    r = E.shoot_to_sonic(1.0, 0.375, rtol=1e-12, atol=1e-14)
    xs = r['x']
    print(f"x_s = {xs:.8f}   sonic (A,N,om,V) = ({r['A']:.6f},{r['N']:.6f},{r['om']:.6f},{r['V']:.6f})")
    print(f"  exact sonic want (1.5, 1.154701, 0.75, -0.577350)")

    # accurate reference background: integrate center -> just before sonic, dense
    z0 = math.exp(-16.0); Y0 = E.center_state3(1.0, 0.375, z0)
    sol = solve_ivp(E.rhs3, [-16.0, xs-1e-8], Y0, method='DOP853',
                    rtol=1e-12, atol=1e-14, dense_output=True)
    def refY(x):
        N, om, V = sol.sol(x); A = E.A_of(N, om, V)
        return np.array([A, N, om, V])

    # exact series about the sonic point in t = x - x_s
    V0s = -1/sp.sqrt(3)
    print(f"\nbuilding exact series at V0=-1/sqrt(3), order={order} ...", flush=True)
    info = SS.build_series(V0s, order=order)
    print(f"quadratic branch at order p={info['p_quad']}, {len(info['branches'])} real root(s)")
    ts = [-0.02, -0.05, -0.10, -0.15, -0.20]
    for bi, br in enumerate(info['branches']):
        o1 = br.get(info['ncO'][1]); v1 = br.get(info['ncV'][1])
        full = SS.finish_branch(info, br)
        cA = [complex(sp.N(c, 30)) for c in full['A']]
        cN = [complex(sp.N(c, 30)) for c in full['N']]
        cO = [complex(sp.N(c, 30)) for c in full['om']]
        cV = [complex(sp.N(c, 30)) for c in full['V']]
        def serY(t):
            return np.array([sum(cA[k]*t**k for k in range(len(cA))),
                             sum(cN[k]*t**k for k in range(len(cN))),
                             sum(cO[k]*t**k for k in range(len(cO))),
                             sum(cV[k]*t**k for k in range(len(cV)))]).real
        print(f"\n branch{bi}: om1={float(o1):+.6f} V1={float(v1):+.6f}")
        print(f"   A coeffs: {[f'{c.real:+.5g}' for c in cA]}")
        print(f"   V coeffs: {[f'{c.real:+.5g}' for c in cV]}")
        worst = 0.0
        for t in ts:
            err = np.abs(serY(t) - refY(xs+t)); worst = max(worst, err.max())
            print(f"   t={t:+.3f}: max|series-ref| = {err.max():.2e}")
        print(f"   >>> branch{bi} worst over t in [-0.2,-0.02]: {worst:.2e}")

if __name__ == "__main__":
    import sys
    order = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    main(order)
