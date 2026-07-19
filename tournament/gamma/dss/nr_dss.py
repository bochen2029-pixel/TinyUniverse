# nr_dss.py — RESEARCH SCAFFOLD (pre-contract validation) for the scalar-DSS critical
# solution (Gundlach gr-qc/9604019 = tournament/gamma/GUN96_src/chop2.tex): construct the
# discretely-self-similar Choptuik background and read the echoing period Delta ~ 3.4453.
#
# Status: research code in tournament/ (the sanctioned pre-C++ pattern, as the fluid-beta
# thread). The module contract (contracts/similarity_nexus.contract.md v0.1.0 DRAFT) is
# awaiting operator review; NO substrate/*.cpp is written from this until approval.
#
# FORMULATION (paper sec 2.3 + app A/C/D):
#   fields Z = (a, g, X+, X-) on the cylinder (zeta, tau), tau-period Delta;
#   coordinates tau = ln t, zeta = ln(r/t) - xi0(tau); SSH at zeta = 0 (D0=1 gauge).
#   evolution in zeta:  X±,z = [C± ± E g X±,t] / [1 ± P E g],   E = e^{zeta+xi0}, P = 1+xi0'
#                       g,z  = (1-a^2) g
#   a NOT evolved: reconstructed from the constraint (a,t eq) each zeta — with f = a^{-2}:
#                       f,t + P f = P - P(X+^2+X-^2) - (E g)^{-1}... see a_from_constraint.
#   C± = [(1-a^2)/2 - a^2 X∓^2] X± - X∓.
#   Frequency parity (kappa=0 boundedness): a, g, xi0 carry EVEN harmonics (0,2,4,..);
#   X± carry ODD harmonics (1,3,..). tau-translation gauge: sin part of xi0's k=2 harmonic = 0.
#   BCs: center (zeta -> -inf, app A): a = 1 + a2 e^{2z}, g = 1 - a2 e^{2z}, a2 = Y1^2/3,
#        X = X2 e^{2z}, Y = Y1 e^z + Y3 e^{3z} (X± = X ± Y), free periodic Y1(tau);
#        SSH (zeta = 0, app A): free X+0(tau); g0 = e^{-xi0}/P (D0=1); a0, X-0, X-1 from
#        closed-form periodic linear ODEs; first-order coefficients algebraic.
#   Two-sided implicit-midpoint march (explicit is UNSTABLE, app D) to zeta_m; match
#   (g, X+, X-) harmonics (a matches via the constraint); unknowns {Delta, xi0, Y1, X+0}.
#
# CONTROLS: (1) vacuum: Y1 = X+0 = xi0 = 0 -> a = g = 1, X± = 0 must satisfy everything
# exactly (machine-precision residual); (2) the constraint-vs-evolution consistency of a
# (a,z = C1 monitored along the march); (3) Delta vs 3.4453 +- 0.0005 (Gundlach).
import numpy as np, sys, time
from numpy.fft import rfft, irfft

# ---------------------------------------------------------------- Fourier layer (period D)
class F:
    """Collocation on M uniform tau-points over [0, D). Real fields. All transforms are
    PRECOMPUTED dense matmuls (numpy rfft has ~10x call overhead at M~32; the residual is
    evaluated thousands of times in the search — measured 10x speedup)."""
    def __init__(self, M, D):
        self.M = M; self.D = D
        self.tau = np.arange(M)*D/M
        K = M//2 + 1
        k = np.arange(K)
        self.omega = 2*np.pi*k/D            # rfft frequencies
        n = np.arange(M)
        self.Wf = np.exp(-2j*np.pi*np.outer(k, n)/M)          # forward (== np.fft.rfft)
        c = np.full(K, 2.0); c[0] = 1.0
        if M % 2 == 0: c[-1] = 1.0
        self.Wi = (np.exp(2j*np.pi*np.outer(n, k)/M)*c)/M     # inverse (apply, take .real)
        # masks
        kk = np.arange(K)
        self.mask_even = ((kk % 2 == 0) & (kk <= M//3)).astype(float)
        self.mask_odd  = ((kk % 2 == 1) & (kk <= M//3)).astype(float)
        # one-time self-check vs numpy rfft/irfft
        t = np.cos(2*np.pi*3*n/M) + 0.3*np.sin(2*np.pi*5*n/M)
        assert np.abs(self.Wf@t - rfft(t)).max() < 1e-10
        assert np.abs((self.Wi@(self.Wf@t)).real - t).max() < 1e-10
    def fwd(self, f):  return self.Wf@f
    def inv(self, F_): return (self.Wi@F_).real
    def dtau(self, f):
        return self.inv(1j*self.omega*self.fwd(f))
    def proj(self, f, parity):
        """Keep only even (parity=0) or odd (parity=1) harmonics; also cut the top 1/3."""
        m = self.mask_even if parity == 0 else self.mask_odd
        return self.inv(self.fwd(f)*m)
    def solve_periodic(self, P, h):
        """Unique periodic solution of f' + P f + h = 0 (app C closed form).
        Requires mean(P) != 0 — except the fully degenerate case (P ~ 0 and h ~ 0,
        e.g. the vacuum X-0 equation), where f = 0 is returned."""
        P0 = np.mean(P); Pt = P - P0
        if abs(P0) < 1e-13:
            if np.abs(h).max() < 1e-12:
                return np.zeros_like(h)
            raise FloatingPointError("solve_periodic: mean(P)=0 with h!=0 (no periodic solution)")
        FP = self.fwd(Pt)
        I0 = np.zeros_like(FP); I0[1:] = FP[1:]/(1j*self.omega[1:])
        w = np.exp(self.inv(I0))            # integrating factor (periodic part)
        R = self.fwd(w*h)
        out = R/(P0 + 1j*self.omega)
        # NOTE: when |mean P| is small (the SSH X-0 equation with a0~1), the DC division is
        # near-0/0 — but for parity-constrained fields the DC is a SYMMETRY ZERO; callers
        # must parity-project the result (an epsilon regularization here broke the SSH
        # regularity condition — measured).
        return -self.inv(out)/w

# ---------------------------------------------------------------- the field equations
def a_from_constraint(F_, ab, g, Xp, Xm, xi0, dxi0, zeta):
    """Reconstruct a from the constraint (paper eq 'constraint'), f = a^{-2} linearizes:
       a,t = (Eg)^{-1} C2 + P C1  with  C1 = a/2[(1-a^2)+a^2 S+],  C2 = a^3/2 S-,
       S± = Xp^2 ± Xm^2  =>  f' + P f + [P(S+ - 1) + Q S-] = 0,  Q = e^{-(zeta+xi0)}/g."""
    P = 1.0 + dxi0
    Q = np.exp(-(zeta + xi0))/g
    Sp = Xp*Xp + Xm*Xm
    Sm = Xp*Xp - Xm*Xm
    h = P*(Sp - 1.0) + Q*Sm
    f = F_.solve_periodic(P, h)
    f = F_.proj(f, 0)                       # a^{-2} is EVEN-harmonic (mask includes DC)
    if not np.all(np.isfinite(f)):
        raise FloatingPointError("a_from_constraint: nan")
    # off-branch excursions (f -> 0/negative) are CLAMPED, not raised: a hard cliff gives the
    # outer LM zero gradient; a clamp turns near-branch inconsistency into a large-but-smooth
    # matching mismatch (the true solution has f >= ~0.2 and never engages the clamp)
    return 1.0/np.sqrt(np.maximum(f, 0.01)), (ab is None)

def rhs_z(F_, g, Xp, Xm, xi0, dxi0, zeta):
    """zeta-derivatives of (g, X+, X-) with a from the constraint."""
    a, _ = a_from_constraint(F_, None, g, Xp, Xm, xi0, dxi0, zeta)
    a2 = a*a
    E = np.exp(zeta + xi0)
    P = 1.0 + dxi0
    Cp = (0.5*(1.0 - a2) - a2*Xm*Xm)*Xp - Xm
    Cm = (0.5*(1.0 - a2) - a2*Xp*Xp)*Xm - Xp
    dXp = (Cp + E*g*F_.dtau(Xp))/(1.0 + P*E*g)
    dXm = (Cm - E*g*F_.dtau(Xm))/(1.0 - P*E*g)
    dg = (1.0 - a2)*g
    return dg, dXp, dXm, a

class MarchFail(FloatingPointError):
    """Raised when a march leaves the physical branch; carries the survived fraction s in
    [0,1) so the residual can return a SLOPED penalty (a flat cliff gives LM zero gradient)."""
    def __init__(self, s):
        super().__init__(f"march failed at fraction {s:.3f}")
        self.s = s

def march(F_, g, Xp, Xm, xi0, dxi0, z0, z1, nstep, picard=4):
    """Implicit-midpoint march zeta z0 -> z1 (app D; explicit is unstable). Picard iteration
    on the midpoint state. Parity-projects each step."""
    h = (z1 - z0)/nstep
    z = z0
    for istep in range(nstep):
        try:
            zm = z + 0.5*h
            gm, Xpm, Xmm = g, Xp, Xm
            for _ in range(picard):
                dg, dXp, dXm, _ = rhs_z(F_, gm, Xpm, Xmm, xi0, dxi0, zm)
                gm = g + 0.5*h*dg
                Xpm = Xp + 0.5*h*dXp
                Xmm = Xm + 0.5*h*dXm
            dg, dXp, dXm, _ = rhs_z(F_, gm, Xpm, Xmm, xi0, dxi0, zm)
            g2 = F_.proj(g + h*dg, 0)
            Xp2 = F_.proj(Xp + h*dXp, 1)
            Xm2 = F_.proj(Xm + h*dXm, 1)
            if not (np.all(np.isfinite(g2)) and np.all(np.isfinite(Xp2)) and np.all(np.isfinite(Xm2))):
                raise FloatingPointError("nan")
            g, Xp, Xm = g2, Xp2, Xm2
        except FloatingPointError as e:
            if isinstance(e, MarchFail): raise
            raise MarchFail(istep/nstep)
        z += h
    return g, Xp, Xm

# ---------------------------------------------------------------- boundary seeds (app A)
def seed_center(F_, Y1, xi0, dxi0, zL):
    """Center expansion at zeta = zL: a=1+a2 e^{2z}, g=1-a2 e^{2z}, X=X2 e^{2z},
    Y=Y1 e^z + Y3 e^{3z}; X± = X ± Y."""
    e1 = np.exp(zL); e2 = e1*e1; e3 = e2*e1
    exi = np.exp(xi0)
    a2c = Y1*Y1/3.0
    X2 = (exi/3.0)*(F_.dtau(Y1) - (1.0 + dxi0)*Y1)
    Y3 = -(2.0/3.0)*Y1**3 + exi*(0.5*F_.dtau(X2) - (1.0 + dxi0)*X2)
    X = X2*e2
    Y = Y1*e1 + Y3*e3
    g = 1.0 - a2c*e2
    return g, X + Y, X - Y

def seed_ssh(F_, Xp0, xi0, dxi0, delta):
    """SSH expansion at zeta = -delta (app A): D0=1 => g0 = e^{-xi0}/P; a0 from
    f=a0^{-2}: f' + P f + P(2 X+0^2 - 1) = 0; X-0 from ODE1; first order: D1, a1, X+1
    algebraic, X-1 from ODE2; g(z) ~ g0 e^{-z}(D0 + D1 z)."""
    P = 1.0 + dxi0
    g0 = np.exp(-xi0)/P
    f = F_.proj(F_.solve_periodic(P, P*(2.0*Xp0*Xp0 - 1.0)), 0)     # a0^{-2}: EVEN
    a0 = 1.0/np.sqrt(np.maximum(f, 1e-12))
    a02 = a0*a0
    # ODE1: X-0' + P[a0^2(1/2+X+0^2)-1/2] X-0 + P X+0 = 0.  X-0 is ODD-harmonic: its DC is a
    # SYMMETRY ZERO — the solve's near-0/0 DC must be projected away (measured blowup source).
    Xm0 = F_.proj(F_.solve_periodic(P*(a02*(0.5 + Xp0*Xp0) - 0.5), P*Xp0), 1)
    D1 = 2.0 - a02
    a1 = 0.5*a0*((1.0 - a02) + a02*(Xp0*Xp0 + Xm0*Xm0))
    D2 = 0.5*(D1*(2.0 - a02) - 2.0*a0*a1)          # D'' /2 from D' = D(2-a^2)
    Xp1 = 0.5*(-a02*Xp0*(0.5 + Xm0*Xm0) + 0.5*Xp0 - Xm0 + F_.dtau(Xp0)/P)
    # ODE2: X-1' + P[-5/2 + a0^2(3/2+X+0^2)] X-1
    #       + P[2 a0 a1 X-0 (1/2+X+0^2) + 2 a0^2 X-0 X+0 X+1 + X+1] + (2-a0^2) X-0' = 0
    hh = P*(2.0*a0*a1*Xm0*(0.5 + Xp0*Xp0) + 2.0*a02*Xm0*Xp0*Xp1 + Xp1) \
         + (2.0 - a02)*F_.dtau(Xm0)
    Xm1 = F_.proj(F_.solve_periodic(P*(-2.5 + a02*(1.5 + Xp0*Xp0)), hh), 1)
    z = -delta
    g = g0*np.exp(-z)*(1.0 + D1*z + D2*z*z)
    Xp = Xp0 + Xp1*z
    Xm = Xm0 + Xm1*z
    return g, Xp, Xm

# ---------------------------------------------------------------- packing + residual
class Prob:
    """Unknowns: Delta; xi0 harmonics (even: c0,c2,s2*,c4,s4; s2 = 0 fixed = tau gauge);
    Y1 and X+0 harmonics (odd: c1,s1,c3,s3). Total 1 + 4 + 4 + 4 = 13."""
    def __init__(self, M=32, zL=-6.0, delta=0.02, zm=-2.0, nstepL=90, nstepR=60):
        self.M = M; self.zL = zL; self.delta = delta; self.zm = zm
        self.nstepL = nstepL; self.nstepR = nstepR
    def fields(self, u):
        D = u[0]
        F_ = F(self.M, D)
        t = F_.tau; w = 2*np.pi/D
        xi0 = u[1] + u[2]*np.cos(2*w*t) + u[3]*np.cos(4*w*t) + u[4]*np.sin(4*w*t)
        Y1  = u[5]*np.cos(w*t) + u[6]*np.sin(w*t) + u[7]*np.cos(3*w*t) + u[8]*np.sin(3*w*t)
        Xp0 = u[9]*np.cos(w*t) + u[10]*np.sin(w*t) + u[11]*np.cos(3*w*t) + u[12]*np.sin(3*w*t)
        return F_, xi0, F_.dtau(xi0), Y1, Xp0
    def march_right(self, F_, gR, XpR, XmR, xi0, dxi0):
        """SSH -> zm with GEOMETRICALLY GRADED steps near zeta=0: the X- equation's
        denominator (1-D) ~ |D1 zeta| makes the stiffness ~ 1/|zeta| — uniform steps
        overshoot at the second node (measured)."""
        z = -self.delta
        while z > -0.3:
            z2 = max(1.25*z, -0.3)          # z negative: 1.25*z steps outward geometrically
            gR, XpR, XmR = march(F_, gR, XpR, XmR, xi0, dxi0, z, z2, 2)
            z = z2
        n = max(6, int((z - self.zm)/0.03))
        return march(F_, gR, XpR, XmR, xi0, dxi0, z, self.zm, n)

    def residual(self, u):
        sL = sR = 0.0
        try:
            F_, xi0, dxi0, Y1, Xp0 = self.fields(u)
            gL, XpL, XmL = seed_center(F_, Y1, xi0, dxi0, self.zL)
            gL, XpL, XmL = march(F_, gL, XpL, XmL, xi0, dxi0, self.zL, self.zm, self.nstepL)
            sL = 1.0
            gR, XpR, XmR = seed_ssh(F_, Xp0, xi0, dxi0, self.delta)
            gR, XpR, XmR = self.march_right(F_, gR, XpR, XmR, xi0, dxi0)
            sR = 1.0
            if not (np.all(np.isfinite(gL)) and np.all(np.isfinite(gR))):
                raise FloatingPointError("march nan")
        except MarchFail as e:
            if sL == 0.0: sL = e.s
            else: sR = e.s
            return np.full(13, 8.0*(3.0 - sL - sR))    # SLOPED cliff: deeper survival = lower
        except FloatingPointError:
            return np.full(13, 24.0)
        # (march_right's per-segment fractions make the cliff piecewise; LM needs
        #  diff_step large enough to see it — set in attempt())
        def harm(f, parity, kmax):
            R = F_.fwd(f)/self.M
            out = []
            for k in range(0, kmax+1):
                if k % 2 != parity: continue
                out.append(R[k].real)
                if k > 0: out.append(R[k].imag)
            return out
        r = harm(gL - gR, 0, 4) + harm(XpL - XpR, 1, 3) + harm(XmL - XmR, 1, 3)
        return np.array(r)

def vacuum_control(M=32):
    """Y1 = X+0 = xi0 = 0: flat space. Residual must be ~machine zero."""
    p = Prob(M=M)
    u = np.zeros(13); u[0] = 3.4453
    r = p.residual(u)
    print(f"[vacuum] max|residual| = {np.abs(r).max():.2e}   (must be ~1e-14)")
    return np.abs(r).max()

def attempt(seed=None, M=32, verbose=True, max_nfev=600, fix_delta=False):
    from scipy.optimize import least_squares
    p = Prob(M=M)
    if seed is None:
        u0 = np.zeros(13)
        u0[0] = 3.4453           # Delta (lit)
        u0[9], u0[10] = 0.35, 0.0   # X+0 ~ 0.35 cos(w tau)
        u0[5], u0[6] = 0.0, 0.35    # Y1  ~ 0.35 sin(w tau)
    else:
        u0 = np.array(seed)
    lo = np.full(13, -4.0); hi = np.full(13, 4.0)
    lo[0], hi[0] = (u0[0]-1e-9, u0[0]+1e-9) if fix_delta else (3.0, 3.9)
    t0 = time.time()
    sol = least_squares(p.residual, u0, method='trf', bounds=(lo, hi),
                        x_scale=np.concatenate([[0.1], np.full(12, 0.5)]),
                        diff_step=0.02,     # cliff-visible finite differences (survival slope)
                        xtol=1e-14, ftol=1e-12, max_nfev=max_nfev, verbose=2 if verbose else 0)
    if verbose:
        print(f"[attempt] {time.time()-t0:.0f}s  status={sol.status}  |r|={np.linalg.norm(sol.fun):.2e}")
        print(f"  Delta = {sol.x[0]:.6f}   (Gundlach 3.4453 +- 0.0005)")
        print(f"  xi0   = {sol.x[1:5].round(4)}")
        print(f"  Y1    = {sol.x[5:9].round(4)}")
        print(f"  X+0   = {sol.x[9:13].round(4)}")
    return sol

def basin_search(M=32, Dfix=3.4453, topn=8,
                 XI0=(-1.2, -0.8, -0.4, 0.0), AYs=(0.4, 0.6, 1.0, 1.5),
                 AXs=(0.1, 0.2, 0.3, 0.45, 0.7),
                 PHR=(-2.75, -2.36, -1.96, -1.57, -1.18, -0.79, 0.0, 0.79, 1.57, 2.36)):
    """Coarse global search over seed amplitudes at FIXED Delta (lit value). The residual
    is invariant under a global tau-shift (measured: 4-fold degeneracy in round 1), so the
    GLOBAL phase is fixed (Y1 pure sin) and only the RELATIVE phase of X+0 is scanned."""
    p = Prob(M=M)
    t0 = time.time(); results = []
    n = 0
    for x0c in XI0:
        for AY in AYs:
            for AX in AXs:
                for ph in PHR:
                    u = np.zeros(13); u[0] = Dfix
                    u[1] = x0c
                    u[5], u[6] = 0.0, AY                       # Y1 = AY sin(w tau)
                    u[9], u[10] = AX*np.sin(ph), AX*np.cos(ph) # X+0 phase-shifted
                    r = p.residual(u)
                    results.append((np.linalg.norm(r), x0c, AY, AX, ph, u.copy()))
                    n += 1
    results.sort(key=lambda t: t[0])
    print(f"[basin] {n} evals in {time.time()-t0:.0f}s (fixed Delta={Dfix}); top {topn}:")
    for (rn, x0c, AY, AX, ph, _) in results[:topn]:
        print(f"  |r|={rn:9.4f}  xi0c0={x0c:+.2f}  A_Y1={AY:.2f}  A_X+0={AX:.2f}  rel-ph={ph:+.2f}")
    return results

def continuation(M=32, Dfix=3.4453, cs=(0.4, 0.7, 1.0, 1.4, 1.9, 2.5), max_nfev=450):
    """Amplitude-constrained continuation: the UNCONSTRAINED problem's global minimum is the
    VACUUM (measured: the basin search slides to zero amplitude — flat space matches at the
    1.35e-6 seed floor). Pin RMS(Y1) = c by an appended penalty residual and solve the
    constrained LM at each c: the isolated DSS solution appears as a dip of min|r| at its
    true amplitude. Delta held at the literature value during the scan."""
    from scipy.optimize import least_squares
    p = Prob(M=M)
    out = []
    for c in cs:
        def rc(u):
            base = p.residual(u)
            F_, xi0, dxi0, Y1, Xp0 = p.fields(u)
            rms = np.sqrt(np.mean(Y1*Y1))
            return np.concatenate([base, [10.0*(rms - c)]])
        u0 = np.zeros(13); u0[0] = Dfix
        u0[6] = c*np.sqrt(2.0)          # Y1 ~ c*sqrt2 sin -> RMS ~ c
        u0[10] = 0.25*c                 # X+0 seeded at the basin's preferred rel-phase
        lo = np.full(13, -6.0); hi = np.full(13, 6.0)
        lo[0], hi[0] = Dfix - 1e-9, Dfix + 1e-9
        t0 = time.time()
        sol = least_squares(rc, u0, method='trf', bounds=(lo, hi),
                            x_scale=0.5, xtol=1e-13, ftol=1e-11, max_nfev=max_nfev)
        rphys = np.linalg.norm(sol.fun[:-1])
        print(f"  c={c:4.2f}: min|r_phys|={rphys:9.5f}  (pen={sol.fun[-1]:+.1e})  "
              f"{time.time()-t0:4.0f}s  nfev={sol.nfev}", flush=True)
        out.append((c, rphys, sol.x.copy()))
    return out

if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "vac"
    if mode == "vac":
        vacuum_control()
    elif mode == "cont":
        vacuum_control()
        continuation()
    elif mode == "basin":
        vacuum_control()
        basin_search()
    elif mode == "polish":
        # LM polish from the basin's best candidates (Delta freed)
        res = basin_search()
        for i in range(3):
            print(f"\n=== polish candidate {i} (|r0|={res[i][0]:.3f}) ===")
            attempt(seed=res[i][6])
    elif mode == "attempt":
        vacuum_control()
        attempt()
