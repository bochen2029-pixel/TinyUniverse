# nr_relax.py — GLOBAL RELAXATION for the scalar-field DSS critical solution (Gundlach
# gr-qc/9604019 app D: "To solve this system of algebraic equations to machine precision
# together with algebraic boundary conditions on two sides, I use a standard relaxation
# algorithm"). The two-sided SHOOTING formulation is structurally hostile (measured, see
# nr_dss.py + git log): off the solution manifold the a-constraint has no real solution,
# so march-based residuals are cliff-ridden. Relaxation makes every grid value an unknown
# and every equation algebraic — and the evolver cylinder seeds ALL of them directly.
#
# Unknowns:  u = [Delta, xi0-even-coeffs(NE-1; k=2 sine excluded = tau gauge),
#                 per node k=0..Nz-1: g-even(NE), X+-odd(NO), X--odd(NO)]
# Residuals: interior implicit-midpoint equations for (g, X+, X-)   NPF*(Nz-1)
#            center BCs at z_L (app A pointwise relations, no extra Y1 unknowns):
#              g = 1 - Y^2/3 (even NE),  X = e^z (e^xi0/3)(Y' - (1+xi0')Y) (odd NO)
#            SSH at z=0 (the LAST node — the singular point is never evaluated):
#              D0 = (1+xi0') e^{xi0} g - 1 = 0 (even NE)
#              regularity: C- - e^{xi0} g X-,tau = 0 (odd NO)
#   count: NPF*Nz + NE unknowns = NPF*(Nz-1) + (NE+NO) + (NE+NO) residuals -> SQUARE.
# a is reconstructed at every midpoint/node from the constraint (closed-form periodic
# solve, f = a^-2 linear; even-projected + clamped for smooth off-branch excursions).
#
# Truncation is CONFIGURABLE: configure(M, KE, KO). Defaults (32, 10, 9) reproduce the
# session-1 numbers; Gundlach's own content reaches k ~ 21 (2N = 128) — raise for the
# strong-field solve.
import numpy as np, sys, time
from numpy.fft import rfft, irfft

M = 32; KE = 10; KO = 9
NE = 1 + KE                    # c0 + (c,s) per even k>=2:  1 + 2*(KE/2)
NO = KO + 1                    # (c,s) per odd k:           2*((KO+1)/2)
NPF = NE + 2*NO
_k = np.arange(M//2 + 1)

def configure(M_=32, KE_=10, KO_=9):
    """Set the tau-truncation. KE even, KO odd; KE, KO <= M/3 (dealias headroom)."""
    global M, KE, KO, NE, NO, NPF, _k
    assert KE_ % 2 == 0 and KO_ % 2 == 1 and KE_ <= M_//3 + 1 and KO_ <= M_//3 + 1
    M, KE, KO = M_, KE_, KO_
    NE = 1 + KE
    NO = KO + 1
    NPF = NE + 2*NO
    _k = np.arange(M//2 + 1)

def bases(Delta):
    """Evaluation + derivative matrices for even/odd harmonic vectors at M tau-points."""
    tau = np.arange(M)*Delta/M
    w = 2*np.pi/Delta
    Be = np.zeros((M, NE)); Bde = np.zeros((M, NE))
    Be[:, 0] = 1.0
    j = 1
    for k in range(2, KE+1, 2):
        Be[:, j] = np.cos(k*w*tau);  Bde[:, j] = -k*w*np.sin(k*w*tau)
        Be[:, j+1] = np.sin(k*w*tau); Bde[:, j+1] = k*w*np.cos(k*w*tau)
        j += 2
    Bo = np.zeros((M, NO)); Bdo = np.zeros((M, NO))
    j = 0
    for k in range(1, KO+1, 2):
        Bo[:, j] = np.cos(k*w*tau);  Bdo[:, j] = -k*w*np.sin(k*w*tau)
        Bo[:, j+1] = np.sin(k*w*tau); Bdo[:, j+1] = k*w*np.cos(k*w*tau)
        j += 2
    Pe = np.linalg.pinv(Be); Po = np.linalg.pinv(Bo)
    return Be, Bde, Bo, Bdo, Pe, Po

# ---- batched periodic solver on VALUE rows: f' + P f + h = 0 per row -------------------
def solve_periodic_rows(P, h, Delta):
    om = 2*np.pi*_k/Delta
    P0 = P.mean(axis=1, keepdims=True)
    FP = rfft(P - P0, axis=1)
    I0 = np.zeros_like(FP); I0[:, 1:] = FP[:, 1:]/(1j*om[1:])
    w = np.exp(irfft(I0, M, axis=1))
    R = rfft(w*h, axis=1)
    out = R/(P0 + 1j*om[None, :])
    return -irfft(out, M, axis=1)/w

class Relax:
    def __init__(self, Nz=40, zL=-6.0):
        self.Nz = Nz; self.zL = zL
        self.z = np.linspace(zL, 0.0, Nz)
        self.h = self.z[1] - self.z[0]
        self.NE, self.NO, self.NPF, self.M = NE, NO, NPF, M
        self.nu = 1 + (NE - 1) + NPF*Nz

    # ---- packing ----
    def unpack(self, u):
        Delta = u[0]
        x = u[1:NE]                              # NE-1 xi0 coeffs (s2 excluded)
        xi0c = np.zeros(NE)
        xi0c[0] = x[0]; xi0c[1] = x[1]           # c0, c2  (index 2 = s2 = 0 gauge)
        xi0c[3:] = x[2:]
        F = u[NE:].reshape(self.Nz, NPF)
        gc = F[:, :NE]; xpc = F[:, NE:NE+NO]; xmc = F[:, NE+NO:]
        return Delta, xi0c, gc, xpc, xmc

    def fields(self, u):
        Delta, xi0c, gc, xpc, xmc = self.unpack(u)
        B = bases(Delta)
        Be, Bde, Bo, Bdo, Pe, Po = B
        xi0 = Be@xi0c; dxi0 = Bde@xi0c
        G = gc@Be.T
        Xp = xpc@Bo.T; Xm = xmc@Bo.T
        dXp = xpc@Bdo.T; dXm = xmc@Bdo.T
        return Delta, xi0, dxi0, G, Xp, Xm, dXp, dXm, B

    def a_rows(self, G, Xp, Xm, xi0, dxi0, zrow, Delta):
        P = np.broadcast_to(1.0 + dxi0, G.shape).copy()
        Q = np.exp(-(zrow[:, None] + xi0[None, :]))/G
        Sp = Xp*Xp + Xm*Xm; Sm = Xp*Xp - Xm*Xm
        hh = P*(Sp - 1.0) + Q*Sm
        f = solve_periodic_rows(P, hh, Delta)
        Ff = rfft(f, axis=1)
        mask = ((_k % 2 == 0) & (_k <= M//3)).astype(float)
        f = irfft(Ff*mask, M, axis=1)
        return 1.0/np.sqrt(np.maximum(f, 5e-3))

    def rhs_rows(self, G, Xp, Xm, dXp, dXm, xi0, dxi0, zrow, Delta):
        a = self.a_rows(G, Xp, Xm, xi0, dxi0, zrow, Delta)
        a2 = a*a
        E = np.exp(zrow[:, None] + xi0[None, :])
        P = 1.0 + dxi0
        Cp = (0.5*(1.0 - a2) - a2*Xm*Xm)*Xp - Xm
        Cm = (0.5*(1.0 - a2) - a2*Xp*Xp)*Xm - Xp
        dG = (1.0 - a2)*G
        dXpz = (Cp + E*G*dXp)/(1.0 + P[None, :]*E*G)
        dXmz = (Cm - E*G*dXm)/(1.0 - P[None, :]*E*G)
        return dG, dXpz, dXmz

    def residual(self, u, pin=None):
        """pin=(c, w): append w*(RMS(X+ at the SSH node) - c) — the Gundlach normalization
        device, used to EXCLUDE THE VACUUM (measured: the unpinned Newton slid to flat
        space; the 1e-28 residual was the tell — a real solution floors at truncation)."""
        try:
            Delta, xi0, dxi0, G, Xp, Xm, dXp, dXm, B = self.fields(u)
            Be, Bde, Bo, Bdo, Pe, Po = B
            Nz, h = self.Nz, self.h
            Gm = 0.5*(G[:-1] + G[1:]); Xpm = 0.5*(Xp[:-1] + Xp[1:]); Xmm = 0.5*(Xm[:-1] + Xm[1:])
            dXpm = 0.5*(dXp[:-1] + dXp[1:]); dXmm = 0.5*(dXm[:-1] + dXm[1:])
            zm = 0.5*(self.z[:-1] + self.z[1:])
            dG, dXpz, dXmz = self.rhs_rows(Gm, Xpm, Xmm, dXpm, dXmm, xi0, dxi0, zm, Delta)
            Rg = (G[1:] - G[:-1]) - h*dG
            Rp = (Xp[1:] - Xp[:-1]) - h*dXpz
            Rm = (Xm[1:] - Xm[:-1]) - h*dXmz
            rg = Rg@Pe.T; rp = Rp@Po.T; rm = Rm@Po.T
            Y0 = 0.5*(Xp[0] - Xm[0]); X0 = 0.5*(Xp[0] + Xm[0])
            dY0 = 0.5*(dXp[0] - dXm[0])
            bc_g = G[0] - (1.0 - (Y0*Y0)/3.0)
            bc_X = X0 - np.exp(self.zL)*(np.exp(xi0)/3.0)*(dY0 - (1.0 + dxi0)*Y0)
            rbc_g = Pe@bc_g
            rbc_X = Po@bc_X
            aN = self.a_rows(G[-1:], Xp[-1:], Xm[-1:], xi0, dxi0, self.z[-1:], Delta)[0]
            a2N = aN*aN
            D0 = (1.0 + dxi0)*np.exp(xi0)*G[-1] - 1.0
            CmN = (0.5*(1.0 - a2N) - a2N*Xp[-1]*Xp[-1])*Xm[-1] - Xp[-1]
            REG = CmN - np.exp(xi0)*G[-1]*dXm[-1]
            r_ssh_e = Pe@D0
            r_ssh_o = Po@REG
            out = [rg.ravel(), rp.ravel(), rm.ravel(), rbc_g, rbc_X, r_ssh_e, r_ssh_o]
            if pin is not None:
                c, w = pin
                out.append(np.array([w*(np.sqrt(np.mean(Xp[-1]**2)) - c)]))
            return np.concatenate(out)
        except FloatingPointError:
            n = NPF*(self.Nz - 1) + 2*(NE + NO) + (1 if pin is not None else 0)
            return np.full(n, 30.0)

    # ---- sparsity for grouped finite differences ----
    def sparsity(self, pin=False):
        from scipy.sparse import lil_matrix
        nr = NPF*(self.Nz - 1) + 2*(NE + NO) + (1 if pin else 0)
        Sp = lil_matrix((nr, self.nu), dtype=np.int8)
        Sp[:, :NE] = 1                                   # Delta + xi0: dense columns
        for k in range(self.Nz - 1):
            cols = slice(NE + NPF*k, NE + NPF*(k + 2))
            Sp[NE*k: NE*(k+1), cols] = 1
            r1 = NE*(self.Nz - 1)
            Sp[r1 + NO*k: r1 + NO*(k+1), cols] = 1
            r2 = r1 + NO*(self.Nz - 1)
            Sp[r2 + NO*k: r2 + NO*(k+1), cols] = 1
        rb = NPF*(self.Nz - 1)
        Sp[rb:rb + NE + NO, NE:NE + NPF] = 1                     # center BCs: node 0
        Sp[rb + NE + NO:, NE + NPF*(self.Nz - 1):] = 1          # SSH (+pin): last node
        return Sp

def seed_from_cylinder(rx, cyl, xi0_tau, Delta):
    """Interpolate the evolver cylinder onto the relaxation grid IN THE GAUGED FRAME
    (zeta_G = zeta_raw - xi0(tau)), tau-rotated so xi0's k=2 sine vanishes (the gauge)."""
    Be, Bde, Bo, Bdo, Pe, Po = bases(Delta)
    zc = cyl['z']
    Mc = len(xi0_tau)
    th = 2*np.pi*np.arange(Mc)/Mc
    c2 = 2*np.mean(xi0_tau*np.cos(2*th)); s2 = 2*np.mean(xi0_tau*np.sin(2*th))
    ph2 = 0.5*np.arctan2(s2, c2)
    def rot(f):
        F_ = np.fft.rfft(f)
        kk = np.arange(len(F_))
        return np.fft.irfft(F_*np.exp(1j*kk*ph2), len(f))
    xi0_tau = rot(xi0_tau)
    cyl = dict(cyl)
    for key in ('X', 'Y', 'g', 'a'):
        cyl[key] = np.apply_along_axis(rot, 0, np.array(cyl[key]))
    # resample the cylinder tau-grid onto M points if needed
    if Mc != M:
        tsrc = np.arange(Mc)/Mc; tdst = np.arange(M)/M
        xi0_tau = np.interp(tdst, np.concatenate([tsrc, [1.0]]), np.concatenate([xi0_tau, xi0_tau[:1]]))
        for key in ('X', 'Y', 'g', 'a'):
            A = cyl[key]
            A2 = np.empty((M, A.shape[1]))
            for jz in range(A.shape[1]):
                col = A[:, jz]
                A2[:, jz] = np.interp(tdst, np.concatenate([tsrc, [1.0]]), np.concatenate([col, col[:1]]))
            cyl[key] = A2
    u = np.zeros(rx.nu)
    u[0] = Delta
    xc = Pe@xi0_tau
    u[1] = xc[0]; u[2] = xc[1]; u[3:NE] = xc[3:]
    F = np.zeros((rx.Nz, NPF))
    for k, zG in enumerate(rx.z):
        gv = np.empty(M); xv = np.empty(M); yv = np.empty(M)
        for i in range(M):
            zr = zG + xi0_tau[i]
            gv[i] = np.interp(zr, zc, cyl['g'][i])
            xv[i] = np.interp(zr, zc, cyl['X'][i])
            yv[i] = np.interp(zr, zc, cyl['Y'][i])
        F[k, :NE] = Pe@gv
        F[k, NE:NE+NO] = Po@(xv + yv)
        F[k, NE+NO:] = Po@(xv - yv)
    u[NE:] = F.ravel()
    return u

def upsample_u(rx1, u1, rx2):
    """Interpolate a solution between Nz grids (same truncation)."""
    u2 = np.zeros(rx2.nu)
    u2[:NE] = u1[:NE]
    F1 = u1[NE:].reshape(rx1.Nz, NPF)
    F2 = np.zeros((rx2.Nz, NPF))
    for c in range(NPF):
        F2[:, c] = np.interp(rx2.z, rx1.z, F1[:, c])
    u2[NE:] = F2.ravel()
    return u2

def vacuum_control(Nz=40):
    rx = Relax(Nz=Nz)
    u = np.zeros(rx.nu); u[0] = 3.4453
    for k in range(rx.Nz):
        u[NE + NPF*k] = 1.0                    # g = 1 (c0)
    r = rx.residual(u)
    print(f"[vac-relax] M={M} KE={KE} KO={KO}: |r|max = {np.abs(r).max():.2e}")
    return np.abs(r).max()

if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "vac"
    if mode == "vac":
        vacuum_control()
    elif mode == "vac48":
        configure(48, 14, 13)
        vacuum_control()
