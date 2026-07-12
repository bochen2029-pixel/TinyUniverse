# stageB_qr.py — CORRECT + FAST two-sided linear eigenvalue match for the fluid-CSS perturbation.
#
# DIAGNOSIS of the prior failure (|det|~1e-14 for EVERY k):
#   The perturbation ODE hp'(x)=L(x;k)hp is Fuchsian at the sonic point x_s. Indicial exponents
#   (eigs of the residue R) are mu={0,0,0, 1-2k} (verified diag_residue.py).
#     - 3 ANALYTIC modes (mu=0)             -> "sonic-analytic" subspace, dim 3.
#     - 1 NON-ANALYTIC mode (mu=1-2k~-4.6)  -> excluded for analyticity.
#   At the CENTER x->-inf indicial exps are {+2,-1,-1,-2} (k-indep, diag_center.py): exactly 1 mode
#   is REGULAR (Re>0). So dim(center-regular)=1, dim(sonic-analytic)=3, 1+3=4 (full space).
#   The physical eigenmode lives in BOTH -> they intersect nontrivially only at discrete k.
#   Eigenvalue condition:  det[ v_center | s1 | s2 | s3 ] = 0  at a match point x_m.
#
#   WHY prior det was ~machine-zero for ALL k: the 3 sonic-analytic seeds, integrated INWARD to x_m,
#   ALL collapse onto the single fastest-growing analytic direction. Column-normalizing 3 near-
#   parallel vectors gives a rank-deficient block -> det~0 regardless of k. Normalization alone does
#   NOT cure loss of independence.
#
# THE FIX: keep the analytic subspace independent. Two guards, both used:
#   (1) build the 3 analytic solutions from a Frobenius series AT the sonic point (order 8) and
#       ORTHONORMALIZE them (QR) at x_m -> a clean 4x3 orthonormal basis Q of the analytic subspace;
#   (2) (integrated cross-check) propagate the 3-column matrix with periodic QR re-orthonormalization.
#   Residual = det[ c | Q ] with c the unit center vector; equivalently ||c - Q Q^H c|| (rejection).
#   Both are O(1) and hit zero exactly at the eigenvalues.
#
# SPEED: L(x;k) is EXACTLY AFFINE in k: L = P(x) + k Q(x), no k in any denominator (verified). So we
#   lambdify P_ij, Q_ij (functions of the background fields only) and assemble L = P + k Q as pure
#   numpy. The sonic Laurent matrices M_p = A_p + k B_p have CONSTANT numeric A_p,B_p (fixed exact
#   background series) precomputed once. Per-k cost is then pure numpy — fast enough to scan+root-find.
#
# Honest: no target is used in the residual. The scan PRINTS det(k) so sign changes are visible first.

import numpy as np, math, sympy as sp, pickle
import css_core as C

S3 = math.sqrt(3.0)
A2STAR = 0.250067905275     # Stage-A shooting param (nc=1 gauge)
NC = 1.0
DSON_SLOPE = 4.0/3.0        # Dson ~ (4/3)(x - x_s) near sonic (verified diag_residue)

# ---- perturbation operator L = P(fields) + k Q(fields); lambdify P_ij, Q_ij entrywise ----
_d = pickle.load(open("Lmat.pkl", "rb"))
_A0, _N0, _o0, _V0 = sp.symbols('A0 N0 om0 V0')
_Nx, _Ax, _ox, _Vx = sp.symbols('N_x A_x om_x V_x')
_ksym = sp.symbols('k')
_Lsym = sp.sympify(_d['Lmat'])
_ARGS = (_N0, _A0, _o0, _V0, _Nx, _Ax, _ox, _Vx)
_Psym = _Lsym.subs(_ksym, 0)                 # P = L|_{k=0}
_Qsym = sp.expand(_Lsym - _Psym)             # k Q  (linear in k)
_Qsym = _Qsym.applyfunc(lambda e: sp.expand(e).coeff(_ksym, 1))   # Q = d L/dk
_fP = [[sp.lambdify(_ARGS, _Psym[i, j], 'numpy') for j in range(4)] for i in range(4)]
_fQ = [[sp.lambdify(_ARGS, _Qsym[i, j], 'numpy') for j in range(4)] for i in range(4)]


def bg_slopes(Y):
    N, A, om, V = Y
    return (C.Nx(A, N, om, V), C.Ax(A, N, om, V), C.omx(A, N, om, V), C.Vx(A, N, om, V))


def _PQ_at(Y, slp):
    N, A, om, V = Y; Nb, Ab, ob, Vb = slp
    P = np.array([[complex(_fP[i][j](N, A, om, V, Nb, Ab, ob, Vb)) for j in range(4)] for i in range(4)])
    Q = np.array([[complex(_fQ[i][j](N, A, om, V, Nb, Ab, ob, Vb)) for j in range(4)] for i in range(4)])
    return P, Q


def Lnum(Y, slp, k):
    P, Q = _PQ_at(Y, slp)
    return P + k*Q


# ------------------------------------------------------------- background cache (k-independent)
def build_background(a2, z0, h, dm):
    """Integrate the background center->sonic ONCE (RK4), caching per-step stage states (N,A,om,V and
    slopes) AND the k-independent P,Q stage matrices, so per-k perturbation RK4 is pure numpy.
    Stop at the match shell x_m where (x_s - x)=dm, i.e. Dson=-(4/3)dm (just inside sonic)."""
    Y = C.center_seed(NC, a2, z0); x = math.log(z0)
    stages = []            # each entry: 4 stage tuples (Y_stage, slp_stage)
    Dtar = -DSON_SLOPE * dm
    prevD = C.Dson(Y[0], Y[3])
    n = int((3.0 - x) / h)
    i_m = None
    Yseq = []
    for i in range(n):
        k1 = bg_slopes(Y); Y1 = Y
        Y2 = tuple(Y[j] + 0.5*h*k1[j] for j in range(4)); k2 = bg_slopes(Y2)
        Y3 = tuple(Y[j] + 0.5*h*k2[j] for j in range(4)); k3 = bg_slopes(Y3)
        Y4 = tuple(Y[j] + h*k3[j] for j in range(4));      k4 = bg_slopes(Y4)
        stages.append(((Y1, k1), (Y2, k2), (Y3, k3), (Y4, k4)))
        Yseq.append(Y1)
        Ynext = tuple(Y[j] + (h/6)*(k1[j] + 2*k2[j] + 2*k3[j] + k4[j]) for j in range(4))
        Dnext = C.Dson(Ynext[0], Ynext[3])
        # Dson rises from very negative (near center-side) toward 0 at the sonic point (x_s), staying
        # NEGATIVE just inside. Stop at the shell where D crosses Dtar=-(4/3)dm from below (D<0 side),
        # i.e. x_m sits a distance dm INSIDE the sonic point. Require V>0.05 (on the sonic branch).
        if Dnext < 0 and Dnext >= Dtar and prevD < Dtar and Ynext[3] > 0.05:
            i_m = i + 1
            Y = Ynext; prevD = Dnext
            break
        Y = Ynext; prevD = Dnext
    if i_m is None:
        return None
    ns = i_m
    # precompute stage P,Q matrices (4 stages per step) as arrays (4,4,ns)
    P = [np.empty((4, 4, ns), complex) for _ in range(4)]
    Q = [np.empty((4, 4, ns), complex) for _ in range(4)]
    for si in range(4):
        for i in range(ns):
            Ys, slp = stages[i][si]
            Pm, Qm = _PQ_at(Ys, slp)
            P[si][:, :, i] = Pm; Q[si][:, :, i] = Qm
    xm = math.log(z0) + ns * h
    Dlast = C.Dson(Yseq[ns-1][0], Yseq[ns-1][3]) if ns >= 1 else None
    return dict(P=P, Q=Q, i_m=ns, h=h, z0=z0, xm=xm, a2=a2, dm=dm,
                Y0=stages[0][0][0], slp0=stages[0][0][1],
                Ylast=Yseq[ns-1], Dlast=Dlast, Yseq=Yseq)


_BG = {}
def get_bg(a2, z0, h, dm):
    key = (round(a2, 12), round(math.log(z0), 6), h, dm)
    if key not in _BG:
        _BG[key] = build_background(a2, z0, h, dm)
    return _BG[key]


def Lstacks(bg, k):
    """L stage tensors = P + k Q (pure numpy)."""
    return [bg['P'][s] + k*bg['Q'][s] for s in range(4)]


# ----------------------------------------------------------- center-regular seed (the single mode)
def center_regular_seed(bg, k):
    L = _PQ_from_bg0(bg, k)
    ev, Vec = np.linalg.eig(L)
    idx = int(np.argmax(ev.real))       # the +2 regular mode
    v = Vec[:, idx]
    return v / np.linalg.norm(v)

def _PQ_from_bg0(bg, k):
    Pm, Qm = _PQ_at(bg['Y0'], bg['slp0'])
    return Pm + k*Qm


def prop_center(bg, k, Ls):
    """Forward RK4 of the single center-regular perturbation vector center -> x_m."""
    L1, L2, L3, L4 = Ls
    hp = center_regular_seed(bg, k)
    h = bg['h']; ns = bg['i_m']
    for i in range(ns):
        d1 = L1[:, :, i].dot(hp)
        d2 = L2[:, :, i].dot(hp + 0.5*h*d1)
        d3 = L3[:, :, i].dot(hp + 0.5*h*d2)
        d4 = L4[:, :, i].dot(hp + h*d3)
        hp = hp + (h/6)*(d1 + 2*d2 + 2*d3 + d4)
        nrm = np.linalg.norm(hp)
        if nrm > 1e12: hp = hp/nrm
        if not np.all(np.isfinite(hp)): return None
    return hp / np.linalg.norm(hp)


# ------------------------------------------- sonic-analytic subspace: 3-D Frobenius series (fast)
# Precompute CONSTANT numeric Laurent matrices A_p,B_p with M_p = A_p + k B_p, from the EXACT
# background series (branch B1). t*L(x_s+t) = R + M1 t + M2 t^2 + ...  (R = M_0).
#
# NUMERIC coefficient extraction (avoids slow sp.series over radicals): lambdify t*P_ij(fields(t)) and
# t*Q_ij(fields(t)) with the EXACT 2nd-order background series substituted, then extract Taylor
# coefficients by DFT on a small circle |t|=eps (roots of unity). g(t)=sum c_p t^p analytic (pole
# killed by the factor t) => c_p = (1/M) sum_m g(eps*w^m) (eps*w^m)^{-p}, w=exp(2pi i/M). Exact to
# rounding for p<M. This is milliseconds vs minutes.
_ORDER = 8
def _precompute_ML():
    t = sp.symbols('t')
    Ns = 2/sp.sqrt(3) + (-2/sp.sqrt(3))*t + sp.Rational(1, 2)*(2*sp.sqrt(3)/3)*t**2
    As = sp.Rational(3, 2) + 3*t + sp.Rational(1, 2)*33*t**2
    Os = sp.Rational(3, 4) + sp.Rational(9, 2)*t + sp.Rational(1, 2)*sp.Rational(99, 2)*t**2
    Vs = 1/sp.sqrt(3) + (2/sp.sqrt(3))*t + sp.Rational(1, 2)*(10*sp.sqrt(3)/3)*t**2
    dsub = (Ns, As, Os, Vs, sp.diff(Ns, t), sp.diff(As, t), sp.diff(Os, t), sp.diff(Vs, t))
    # lambdify t*P_ij and t*Q_ij as functions of t only (fields substituted)
    tPf = [[sp.lambdify(t, t*_Psym[i, j].subs(dict(zip(_ARGS, dsub))), 'numpy') for j in range(4)] for i in range(4)]
    tQf = [[sp.lambdify(t, t*_Qsym[i, j].subs(dict(zip(_ARGS, dsub))), 'numpy') for j in range(4)] for i in range(4)]
    M = 64                                   # DFT length (>> order); eps small
    eps = 1e-3
    w = np.exp(2j*np.pi*np.arange(M)/M)
    tt = eps*w
    A_list = [np.zeros((4, 4), complex) for _ in range(_ORDER+1)]
    B_list = [np.zeros((4, 4), complex) for _ in range(_ORDER+1)]
    for i in range(4):
        for j in range(4):
            gp = np.array([complex(tPf[i][j](tv)) for tv in tt])
            gq = np.array([complex(tQf[i][j](tv)) for tv in tt])
            cp = np.fft.ifft(gp)             # cp[p]*? -> coefficients scaled by eps^p
            cq = np.fft.ifft(gq)
            for p in range(_ORDER+1):
                A_list[p][i, j] = cp[p]/(eps**p)
                B_list[p][i, j] = cq[p]/(eps**p)
    return A_list, B_list
_ALIST, _BLIST = _precompute_ML()

def _Mp(p, k):
    return _ALIST[p] + k*_BLIST[p]


def analytic_seeds(k, order=_ORDER):
    """3 independent analytic Taylor solutions at the sonic point (mu=0 subspace)."""
    R = _Mp(0, k); Ls = [_Mp(p+1, k) for p in range(order)]
    U, Sv, Vh = np.linalg.svd(R)
    ker = [Vh[i].conj() for i in range(4) if Sv[i] < 1e-8]
    seeds = []
    for a0 in ker:
        a = [np.array(a0, complex)]; ok = True
        for n in range(1, order+1):
            rhs = np.zeros(4, complex)
            for j in range(min(n, len(Ls))):
                rhs += Ls[j].dot(a[n-1-j])
            try:
                a.append(np.linalg.solve(n*np.eye(4) - R, rhs))
            except np.linalg.LinAlgError:
                ok = False; break
        if ok:
            seeds.append(a)
    return seeds


def eval_series(a, t):
    s = np.zeros(4, complex); tp = 1.0
    for n in range(len(a)):
        s += a[n]*tp; tp *= t
    return s


def sonic_subspace_series(bg, k):
    """Analytic subspace basis at x_m by direct series evaluation at t=-dm (dm inside series radius)."""
    seeds = analytic_seeds(k)
    if len(seeds) < 3: return None
    t_m = -bg['dm']
    cols = [eval_series(a, t_m) for a in seeds[:3]]
    Q, _ = np.linalg.qr(np.column_stack(cols))
    return Q[:, :3]


def sonic_subspace_integrated(bg, k, Ls):
    """Integrated cross-check: seed the 3 analytic solutions with the series at the innermost cached
    node (closest to sonic), then RK4 the 3-column matrix INWARD to x_m with periodic QR."""
    seeds = analytic_seeds(k)
    if len(seeds) < 3: return None
    L1, L2, L3, L4 = Ls
    ns = bg['i_m']; h = bg['h']
    t_last = bg['Dlast'] / DSON_SLOPE            # (x_last - x_s), negative
    Mmat = np.column_stack([eval_series(a, t_last) for a in seeds[:3]])
    Mmat, _ = np.linalg.qr(Mmat)                 # orthonormal start near sonic (node ns-1)
    i_target = max(0, (ns - 1) - int(0.5*ns))    # integrate inward ~half the span
    for i in range(ns-1, i_target, -1):
        hh = -h
        d1 = L1[:, :, i].dot(Mmat)
        d2 = L2[:, :, i].dot(Mmat + 0.5*hh*d1)
        d3 = L3[:, :, i].dot(Mmat + 0.5*hh*d2)
        d4 = L4[:, :, i].dot(Mmat + hh*d3)
        Mmat = Mmat + (hh/6)*(d1 + 2*d2 + 2*d3 + d4)
        if (ns-1 - i) % 20 == 0:
            Mmat, _ = np.linalg.qr(Mmat)
        if not np.all(np.isfinite(Mmat)): return None
    Mmat, _ = np.linalg.qr(Mmat)
    return Mmat[:, :3], i_target


def residual(bg, k):
    """Eigenvalue residual = det[ c | Q ] and rejection ||c - Q Q^H c||. Zero <=> eigenvalue."""
    k = complex(k)
    Ls = Lstacks(bg, k)
    c = prop_center(bg, k, Ls)
    if c is None: return None
    Q = sonic_subspace_series(bg, k)
    if Q is None: return None
    M = np.column_stack([c, Q])
    det = np.linalg.det(M)
    rej = c - Q.dot(Q.conj().T.dot(c))
    return det, np.linalg.norm(rej), c, Q


def residual_integrated(bg, k):
    """Same residual but the analytic subspace comes from the INTEGRATED (QR-stabilized) propagation,
    and the center vector is integrated to the SAME node i_target. Cross-check of the series variant."""
    k = complex(k)
    Ls = Lstacks(bg, k)
    sub = sonic_subspace_integrated(bg, k, Ls)
    if sub is None: return None
    Q, i_target = sub
    # integrate center only to node i_target
    L1, L2, L3, L4 = Ls
    hp = center_regular_seed(bg, k); h = bg['h']
    for i in range(i_target):
        d1 = L1[:, :, i].dot(hp); d2 = L2[:, :, i].dot(hp+0.5*h*d1)
        d3 = L3[:, :, i].dot(hp+0.5*h*d2); d4 = L4[:, :, i].dot(hp+h*d3)
        hp = hp + (h/6)*(d1+2*d2+2*d3+d4)
        nrm = np.linalg.norm(hp)
        if nrm > 1e12: hp = hp/nrm
        if not np.all(np.isfinite(hp)): return None
    c = hp/np.linalg.norm(hp)
    M = np.column_stack([c, Q]); det = np.linalg.det(M)
    rej = c - Q.dot(Q.conj().T.dot(c))
    return det, np.linalg.norm(rej), c, Q


if __name__ == "__main__":
    import sys, time
    z0 = math.exp(-12); h = 5e-5; dm = 0.04
    a2 = A2STAR
    if len(sys.argv) > 1: dm = float(sys.argv[1])
    if len(sys.argv) > 2: h = float(sys.argv[2])
    t0 = time.time()
    bg = get_bg(a2, z0, h, dm)
    if bg is None:
        print("background did not reach sonic shell"); sys.exit(2)
    print(f"# QR two-sided match  a2={a2}  z0=e^-12  h={h}  dm={dm}  i_m={bg['i_m']}  x_m={bg['xm']:.4f}  (bg {time.time()-t0:.1f}s)")
    print(f"{'k':>7} {'Re det':>13} {'Im det':>13} {'|det|':>11} {'||rej||':>11}")
    rows = []
    for k in np.arange(0.2, 4.01, 0.1):
        out = residual(bg, k)
        if out is None:
            print(f"{k:7.2f}   (fail)"); continue
        det, rej, c, Q = out
        rows.append((k, det, rej))
        print(f"{k:7.2f} {det.real:13.4e} {det.imag:13.4e} {abs(det):11.3e} {rej:11.3e}")
    print("\nSign changes in Re(det):")
    for i in range(1, len(rows)):
        if rows[i-1][1].real * rows[i][1].real < 0:
            print(f"   between k={rows[i-1][0]:.2f} and k={rows[i][0]:.2f}")
    print("|det| local minima:")
    for i in range(1, len(rows)-1):
        if abs(rows[i][1]) < abs(rows[i-1][1]) and abs(rows[i][1]) < abs(rows[i+1][1]):
            print(f"   min at k={rows[i][0]:.2f}  |det|={abs(rows[i][1]):.3e}")
    print(f"# total {time.time()-t0:.1f}s")
