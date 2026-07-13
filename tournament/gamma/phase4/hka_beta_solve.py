# hka_beta_solve.py — fluid-beta warm-up eigenvalue via the AUTHORITATIVE HKA99 operator.
#
# Method (robust, physically grounded): build the UNIQUE physical analytic solution at the sonic
# point (3 Frobenius modes fixed by gauge Nbar_p(x_s)=0 + the algebraic identity eq:alg-PP), then
# integrate it TOWARD the center (x decreasing). If it contains the e^{-2x} growing-at-center mode
# it blows up by ~e^{2|x_c-x_s|}; at an eigenvalue that mode's coefficient is 0 and the solution stays
# bounded. So Delta(kappa) = log10 |y(x_c)| has a SHARP MINIMUM at each eigenvalue. Gauge mode shows
# only at kappa = -Nbar_x(x_s) ~ 0.357 (the sonic gauge); the physical fluid-beta mode is kappa~2.81.
import numpy as np, math, sys
from scipy.integrate import solve_ivp
from scipy.linalg import null_space
import hka_pert_core as PC, hka_pert_hka99 as H99
PC.Lnum = H99.Lnum
import hka_beta4 as B, hka_frobenius as FR
B.bg(); B.bg_path(); XS = B.bg()['xs']
G = 4.0/3.0

def Lx(x, k): return H99.Lnum(B.bg_state(x), complex(k))

def identity_coeffs(k, x):
    """eq:alg-PP:  (k-A)Abar_p + cN Nbar_p + com obar_p + cV V_p = 0."""
    A, N, om, V, obx, Vx = [float(np.real(z)) for z in B.bg_state(x)]
    oV2 = 1.0 - V*V
    cA = k - A
    cN = 2.0*G*N*V*om/oV2
    com = 2.0*om*(1.0 + (G-1.0)*V*V + G*N*V)/oV2
    cV = 2.0*G*om*(N*(1.0 + V*V) + 2.0*V)/oV2**2
    return np.array([cA, cN, com, cV], complex)

def v_sonic(k, order=16, dt=-0.02):
    """Physical analytic solution just inside the sonic point (x = x_s + dt, dt<0)."""
    modes, R = FR.analytic_modes(k, B.bg(), order=order)
    if len(modes) < 3: return None
    S0 = [m[0] for m in modes]                       # a_0 = ker(R) basis (t=0 values)
    idc = identity_coeffs(k, XS - 1e-7)
    Cmat = np.array([[S0[i][1] for i in range(3)],                    # gauge: Nbar_p(0)=0
                     [complex(np.dot(idc, S0[i])) for i in range(3)]], complex)  # identity at x_s
    ns = null_space(Cmat, rcond=1e-10)
    if ns.shape[1] < 1: return None
    c = ns[:, 0]
    vs = sum(c[i]*FR.eval_mode(modes[i], dt) for i in range(3))
    return vs/np.linalg.norm(vs)

def growth(k, xc=-14.0, order=16):
    vs = v_sonic(k, order=order)
    if vs is None: return None
    x0 = XS - 0.02
    sol = solve_ivp(lambda x, y: Lx(x, k).dot(y), [x0, xc], vs.astype(complex),
                    method='DOP853', rtol=1e-9, atol=1e-12)
    if not sol.success: return None
    return math.log10(np.linalg.norm(sol.y[:, -1]))

def refine(k0, dk=0.05, xc=-14.0, order=16, iters=40):
    """golden-section minimize growth around k0 (real kappa)."""
    a, b = k0-dk, k0+dk
    gr = (math.sqrt(5)-1)/2
    c, d = b-gr*(b-a), a+gr*(b-a)
    fc, fd = growth(c, xc, order), growth(d, xc, order)
    for _ in range(iters):
        if fc is None or fd is None: return None
        if fc < fd: b, d, fd = d, c, fc; c = b-gr*(b-a); fc = growth(c, xc, order)
        else:       a, c, fc = c, d, fd; d = a+gr*(b-a); fd = growth(d, xc, order)
        if b-a < 1e-9: break
    return 0.5*(a+b)

if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "refine":
        k0 = float(sys.argv[2])
        kstar = refine(k0)
        print(f"refined kappa = {kstar:.8f}   beta = 1/kappa = {1.0/kstar:.8f}   growth={growth(kstar):.2f}")
    else:
        lo, hi, step = (2.0, 3.6, 0.05)
        print(f"x_s={XS:.5f}. Delta=log10|y(x_c)| toward center; SHARP MIN at an eigenvalue (gauge ~0.357).")
        rows = []
        k = lo
        while k <= hi + 1e-9:
            g = growth(k); rows.append((k, g)); print(f"  kappa={k:.3f}  log|y|={g:.3f}" if g is not None else f"  kappa={k:.3f} fail")
            k += step
        good = [(k, g) for k, g in rows if g is not None]
        print("\nlocal minima:")
        for i in range(1, len(good)-1):
            if good[i][1] < good[i-1][1] and good[i][1] < good[i+1][1]:
                print(f"   kappa~{good[i][0]:.3f}  log|y|={good[i][1]:.3f}")
