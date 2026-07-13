# nr_frob.py — accurate sonic Laurent residue R + analytic Frobenius modes, built on the EXACT
# background series (nr_sonic) and the VERIFIED operator (hka_pert_hka99). This replaces the poisoned
# hka_pert_sonic polyfit path. Reconciles the residue indicial exponent mu (was 3 disagreeing values).
import numpy as np, math
import nr_sonic as NS
import hka_pert_core as PC, hka_pert_hka99 as H99
PC.Lnum = H99.Lnum                      # route the shared machinery through the verified operator
import hka_pert_sonic as PS
import hka_frobenius as FR
import hka_ec as E

_XS = None
def xs():
    global _XS
    if _XS is None:
        _XS = E.shoot_to_sonic(1.0, 0.375, rtol=1e-12, atol=1e-14)['x']
    return _XS

_CACHE = {}
def bgser(order=14):
    if order not in _CACHE:
        A, N, om, V = NS.bg_series(order)
        _CACHE[order] = dict(A=np.array(A), N=np.array(N), om=np.array(om), V=np.array(V), xs=xs())
    return _CACHE[order]

def bg_fld(coeffs_dict, t):
    """(A,N,om,V,obar_x,V_x) at t=x-x_s from the exact series (slopes from series derivatives)."""
    A, N, om, V = coeffs_dict['A'], coeffs_dict['N'], coeffs_dict['om'], coeffs_dict['V']
    def ev(a): return sum(a[k]*t**k for k in range(len(a)))
    def dev(a): return sum(k*a[k]*t**(k-1) for k in range(1, len(a)))
    omv = ev(om)
    return (ev(A), ev(N), omv, ev(V), dev(om)/omv, dev(V))

def laurent(kappa, order=10, eps=1e-3):
    """R (t^-1 coeff) and the regular Laurent coeffs of L(t;kappa), from the ACCURATE background."""
    return PS.L_laurent(kappa, bgser(order), order=order, eps=eps)

def residue_report(kappas=(2.81055255, 1.0, 0.35699)):
    print("Residue R indicial exponents (nonzero eig) on the EXACT background:")
    for k in kappas:
        lau = laurent(k, order=8, eps=1e-3); R = lau[0]
        ev = np.linalg.eigvals(R)
        nz = ev[np.argmax(np.abs(ev))]
        print(f"  kappa={k:<10}: eig(R)={np.sort_complex(ev).round(4)}")
        print(f"      nonzero mu={nz.real:+.5f}   vs 1-2k={1-2*k:+.5f}  1+2k={1+2*k:+.5f}  -(1+2k)? {-(1+2*k):+.5f}")

def gate(kappa=2.81055255, order=14, tm=-0.03, eps=1e-3):
    """Frobenius GATE: each analytic mode must solve Psi'=L Psi at the match point to <1e-8."""
    bs = bgser(order)
    modes, R = FR.analytic_modes(kappa, bs, order=order, eps=eps)
    fld = bg_fld(bs, tm); L = H99.Lnum(fld, complex(kappa))
    print(f"kappa={kappa}: {len(modes)} analytic modes; residue eig={np.sort_complex(np.linalg.eigvals(R)).round(3)}")
    worst = 0.0
    for i, a in enumerate(modes):
        Psi = sum(a[n]*tm**n for n in range(len(a)))
        dPsi = sum(n*a[n]*tm**(n-1) for n in range(1, len(a)))
        res = np.linalg.norm(dPsi - L.dot(Psi))/max(np.linalg.norm(Psi), 1e-30)
        worst = max(worst, res)
        print(f"  mode{i}: |Psi'-L Psi|/|Psi| at t={tm} = {res:.2e}")
    print(f"  >>> worst analytic-mode residual = {worst:.2e}  ({'PASS' if worst<1e-8 else 'FAIL'} <1e-8)")
    return worst

if __name__ == "__main__":
    import sys
    residue_report()
    print()
    for tm in (-0.02, -0.03, -0.05):
        gate(2.81055255, order=16, tm=tm)
        print()
