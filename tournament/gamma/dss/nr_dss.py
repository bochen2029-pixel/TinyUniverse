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
    """Collocation on M uniform tau-points over [0, D). Real fields; operations via rfft."""
    def __init__(self, M, D):
        self.M = M; self.D = D
        self.tau = np.arange(M)*D/M
        k = np.arange(M//2 + 1)
        self.omega = 2*np.pi*k/D            # rfft frequencies
    def dtau(self, f):
        return irfft(1j*self.omega*rfft(f), self.M)
    def proj(self, f, parity):
        """Keep only even (parity=0) or odd (parity=1) harmonics; also cut the top 1/3."""
        F_ = rfft(f)
        k = np.arange(len(F_))
        mask = (k % 2 == parity)
        mask &= (k <= self.M//3)            # dealias-lite
        return irfft(F_*mask, self.M)
    def solve_periodic(self, P, h):
        """Unique periodic solution of f' + P f + h = 0 (app C closed form).
        Requires mean(P) != 0 — except the fully degenerate case (P ~ 0 and h ~ 0,
        e.g. the vacuum X-0 equation), where f = 0 is returned."""
        P0 = np.mean(P); Pt = P - P0
        if abs(P0) < 1e-13:
            if np.abs(h).max() < 1e-12:
                return np.zeros_like(h)
            raise FloatingPointError("solve_periodic: mean(P)=0 with h!=0 (no periodic solution)")
        I0P = irfft(np.concatenate([[0.0], (rfft(Pt)[1:]/(1j*self.omega[1:]))]), self.M)
        w = np.exp(I0P)                     # integrating factor (periodic part)
        rhs = w*h
        R = rfft(rhs)
        out = np.zeros_like(R)
        out[0] = R[0]/P0
        out[1:] = R[1:]/(P0 + 1j*self.omega[1:])
        return -irfft(out, self.M)/w

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
    return 1.0/np.sqrt(np.maximum(f, 1e-12)), (ab is None)

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

def march(F_, g, Xp, Xm, xi0, dxi0, z0, z1, nstep, picard=4):
    """Implicit-midpoint march zeta z0 -> z1 (app D; explicit is unstable). Picard iteration
    on the midpoint state. Parity-projects each step."""
    h = (z1 - z0)/nstep
    z = z0
    for _ in range(nstep):
        zm = z + 0.5*h
        gm, Xpm, Xmm = g, Xp, Xm
        for _ in range(picard):
            dg, dXp, dXm, _ = rhs_z(F_, gm, Xpm, Xmm, xi0, dxi0, zm)
            gm = g + 0.5*h*dg
            Xpm = Xp + 0.5*h*dXp
            Xmm = Xm + 0.5*h*dXm
        dg, dXp, dXm, _ = rhs_z(F_, gm, Xpm, Xmm, xi0, dxi0, zm)
        g = F_.proj(g + h*dg, 0)
        Xp = F_.proj(Xp + h*dXp, 1)
        Xm = F_.proj(Xm + h*dXm, 1)
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
    f = F_.solve_periodic(P, P*(2.0*Xp0*Xp0 - 1.0))
    a0 = 1.0/np.sqrt(np.maximum(f, 1e-12))
    a02 = a0*a0
    # ODE1: X-0' + P[a0^2(1/2+X+0^2)-1/2] X-0 + P X+0 = 0
    Xm0 = F_.solve_periodic(P*(a02*(0.5 + Xp0*Xp0) - 0.5), P*Xp0)
    D1 = 2.0 - a02
    a1 = 0.5*a0*((1.0 - a02) + a02*(Xp0*Xp0 + Xm0*Xm0))
    D2 = 0.5*(D1*(2.0 - a02) - 2.0*a0*a1)          # D'' /2 from D' = D(2-a^2)
    Xp1 = 0.5*(-a02*Xp0*(0.5 + Xm0*Xm0) + 0.5*Xp0 - Xm0 + F_.dtau(Xp0)/P)
    # ODE2: X-1' + P[-5/2 + a0^2(3/2+X+0^2)] X-1
    #       + P[2 a0 a1 X-0 (1/2+X+0^2) + 2 a0^2 X-0 X+0 X+1 + X+1] + (2-a0^2) X-0' = 0
    hh = P*(2.0*a0*a1*Xm0*(0.5 + Xp0*Xp0) + 2.0*a02*Xm0*Xp0*Xp1 + Xp1) \
         + (2.0 - a02)*F_.dtau(Xm0)
    Xm1 = F_.solve_periodic(P*(-2.5 + a02*(1.5 + Xp0*Xp0)), hh)
    z = -delta
    g = g0*np.exp(-z)*(1.0 + D1*z + D2*z*z)
    Xp = Xp0 + Xp1*z
    Xm = Xm0 + Xm1*z
    return g, Xp, Xm

# ---------------------------------------------------------------- packing + residual
class Prob:
    """Unknowns: Delta; xi0 harmonics (even: c0,c2,s2*,c4,s4; s2 = 0 fixed = tau gauge);
    Y1 and X+0 harmonics (odd: c1,s1,c3,s3). Total 1 + 4 + 4 + 4 = 13."""
    def __init__(self, M=32, zL=-6.0, delta=0.02, zm=-1.0, nstepL=250, nstepR=60):
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
    def residual(self, u):
        F_, xi0, dxi0, Y1, Xp0 = self.fields(u)
        gL, XpL, XmL = seed_center(F_, Y1, xi0, dxi0, self.zL)
        gL, XpL, XmL = march(F_, gL, XpL, XmL, xi0, dxi0, self.zL, self.zm, self.nstepL)
        gR, XpR, XmR = seed_ssh(F_, Xp0, xi0, dxi0, self.delta)
        gR, XpR, XmR = march(F_, gR, XpR, XmR, xi0, dxi0, -self.delta, self.zm, self.nstepR)
        def harm(f, parity, kmax):
            R = rfft(f)/self.M
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

def attempt(seed=None, M=32, verbose=True):
    from scipy.optimize import least_squares
    p = Prob(M=M)
    if seed is None:
        u0 = np.zeros(13)
        u0[0] = 3.4453           # Delta (lit)
        u0[9], u0[10] = 0.35, 0.0   # X+0 ~ 0.35 cos(w tau)
        u0[5], u0[6] = 0.0, 0.35    # Y1  ~ 0.35 sin(w tau)
    else:
        u0 = np.array(seed)
    t0 = time.time()
    sol = least_squares(p.residual, u0, method='lm', xtol=1e-14, ftol=1e-14, max_nfev=4000)
    if verbose:
        print(f"[attempt] {time.time()-t0:.0f}s  status={sol.status}  |r|={np.linalg.norm(sol.fun):.2e}")
        print(f"  Delta = {sol.x[0]:.6f}   (Gundlach 3.4453 +- 0.0005)")
        print(f"  xi0   = {sol.x[1:5].round(4)}")
        print(f"  Y1    = {sol.x[5:9].round(4)}")
        print(f"  X+0   = {sol.x[9:13].round(4)}")
    return sol

if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "vac"
    if mode == "vac":
        vacuum_control()
    elif mode == "attempt":
        vacuum_control()
        attempt()
