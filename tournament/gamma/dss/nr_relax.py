# nr_relax.py — GLOBAL RELAXATION for the scalar-field DSS critical solution (Gundlach
# gr-qc/9604019 app D: "To solve this system of algebraic equations to machine precision
# together with algebraic boundary conditions on two sides, I use a standard relaxation
# algorithm"). The two-sided SHOOTING formulation is structurally hostile (measured, see
# nr_dss.py + git log): off the solution manifold the a-constraint has no real solution,
# so march-based residuals are cliff-ridden. Relaxation makes every grid value an unknown
# and every equation algebraic — and the evolver cylinder seeds ALL of them directly.
#
# Unknowns:  u = [Delta, xi0-even-coeffs(10; k=2 sine excluded = tau gauge),
#                 per node k=0..Nz-1: g-even(11), X+-odd(10), X--odd(10)]
# Residuals: interior implicit-midpoint equations for (g, X+, X-)   31*(Nz-1)
#            center BCs at z_L (app A pointwise relations, no extra Y1 unknowns):
#              g = 1 - Y^2/3 (even 11),  X = (e^xi0/3)(Y' - (1+xi0')Y) (odd 10)
#              with X=(X+ + X-)/2, Y=(X+ - X-)/2 evaluated from node-0 fields
#            SSH at z=0 (the LAST node — the singular point is never evaluated):
#              D0 = (1+xi0') e^{xi0} g - 1 = 0 (even 11)
#              regularity: C- - e^{xi0} g X-,tau = 0 (odd 10)
#   count: 31 Nz + 11 unknowns = 31(Nz-1) + 21 + 21 residuals  ->  SQUARE.
# a is reconstructed at every midpoint/node from the constraint (closed-form periodic
# solve, f = a^-2 linear; even-projected).
import numpy as np, sys, time
from numpy.fft import rfft, irfft

M = 32                 # tau collocation points
KE = 10                # top even harmonic (c0 + (c2,s2)..(c10,s10) -> 11 coeffs)
KO = 9                 # top odd harmonic  ((c1,s1)..(c9,s9)       -> 10 coeffs)
NE, NO = 11, 10

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
    # projection (values -> coeffs): least-squares pinv is exact here (orthogonal columns)
    Pe = np.linalg.pinv(Be); Po = np.linalg.pinv(Bo)
    return Be, Bde, Bo, Bdo, Pe, Po

# ---- batched periodic solver on VALUE rows: f' + P f + h = 0 per row -------------------
_k = np.arange(M//2 + 1)
def solve_periodic_rows(P, h, Delta):
    """P, h: (R, M) value rows; returns f rows. Closed-form integrating factor per row."""
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
        self.nu = 1 + 10 + 31*Nz

    # ---- packing ----
    def unpack(self, u):
        Delta = u[0]
        x = u[1:11]
        xi0c = np.zeros(NE)
        xi0c[0] = x[0]; xi0c[1] = x[1]          # c0, c2   (s2 = 0 gauge)
        xi0c[3:] = x[2:]                        # c4..s10
        F = u[11:].reshape(self.Nz, 31)
        gc = F[:, :11]; xpc = F[:, 11:21]; xmc = F[:, 21:31]
        return Delta, xi0c, gc, xpc, xmc

    def fields(self, u):
        Delta, xi0c, gc, xpc, xmc = self.unpack(u)
        Be, Bde, Bo, Bdo, Pe, Po = bases(Delta)
        xi0 = Be@xi0c; dxi0 = Bde@xi0c
        G = gc@Be.T                              # (Nz, M) values
        Xp = xpc@Bo.T; Xm = xmc@Bo.T
        dXp = xpc@Bdo.T; dXm = xmc@Bdo.T         # tau-derivative VALUES
        return Delta, xi0, dxi0, G, Xp, Xm, dXp, dXm, (Be, Bde, Bo, Bdo, Pe, Po)

    def a_rows(self, G, Xp, Xm, xi0, dxi0, zrow, Delta):
        """a at value-rows (R, M) with per-row z: constraint f = a^-2 (linear, even-proj)."""
        P = np.broadcast_to(1.0 + dxi0, G.shape).copy()
        Q = np.exp(-(zrow[:, None] + xi0[None, :]))/G
        Sp = Xp*Xp + Xm*Xm; Sm = Xp*Xp - Xm*Xm
        hh = P*(Sp - 1.0) + Q*Sm
        f = solve_periodic_rows(P, hh, Delta)
        # even-projection + clamp (off-branch excursions stay smooth for the outer solver)
        Ff = rfft(f, axis=1)
        mask = ((_k % 2 == 0) & (_k <= M//3)).astype(float)
        f = irfft(Ff*mask, M, axis=1)
        return 1.0/np.sqrt(np.maximum(f, 5e-3))

    def rhs_rows(self, G, Xp, Xm, dXp, dXm, xi0, dxi0, zrow, Delta):
        """(dG, dXp, dXm)/dz at value-rows (R, M)."""
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

    def residual(self, u):
        try:
            Delta, xi0, dxi0, G, Xp, Xm, dXp, dXm, B = self.fields(u)
            Be, Bde, Bo, Bdo, Pe, Po = B
            Nz, h = self.Nz, self.h
            # interior: implicit midpoint between consecutive nodes
            Gm = 0.5*(G[:-1] + G[1:]); Xpm = 0.5*(Xp[:-1] + Xp[1:]); Xmm = 0.5*(Xm[:-1] + Xm[1:])
            dXpm = 0.5*(dXp[:-1] + dXp[1:]); dXmm = 0.5*(dXm[:-1] + dXm[1:])
            zm = 0.5*(self.z[:-1] + self.z[1:])
            dG, dXpz, dXmz = self.rhs_rows(Gm, Xpm, Xmm, dXpm, dXmm, xi0, dxi0, zm, Delta)
            Rg = (G[1:] - G[:-1]) - h*dG
            Rp = (Xp[1:] - Xp[:-1]) - h*dXpz
            Rm = (Xm[1:] - Xm[:-1]) - h*dXmz
            # project to harmonic spaces
            rg = Rg@Pe.T; rp = Rp@Po.T; rm = Rm@Po.T          # (Nz-1, 11/10/10)
            # center BCs at node 0 (app A): Y = (X+-X-)/2, X = (X++X-)/2
            Y0 = 0.5*(Xp[0] - Xm[0]); X0 = 0.5*(Xp[0] + Xm[0])
            dY0 = 0.5*(dXp[0] - dXm[0])
            # In node-0 VALUES (Y = Y1 e^z): g = 1 - Y^2/3  and  X = e^{z}(e^{xi0}/3)(Y' - (1+xi0')Y)
            bc_g = G[0] - (1.0 - (Y0*Y0)/3.0)
            bc_X = X0 - np.exp(self.zL)*(np.exp(xi0)/3.0)*(dY0 - (1.0 + dxi0)*Y0)
            rbc_g = Pe@bc_g
            rbc_X = Po@bc_X
            # SSH at the last node z=0: D0 = (1+xi0') e^{xi0} g - 1 = 0 (even);
            # regularity: C- - e^{xi0} g X-,tau = 0 (odd)
            aN = self.a_rows(G[-1:], Xp[-1:], Xm[-1:], xi0, dxi0, self.z[-1:], Delta)[0]
            a2N = aN*aN
            D0 = (1.0 + dxi0)*np.exp(xi0)*G[-1] - 1.0
            CmN = (0.5*(1.0 - a2N) - a2N*Xp[-1]*Xp[-1])*Xm[-1] - Xp[-1]
            REG = CmN - np.exp(xi0)*G[-1]*dXm[-1]
            r_ssh_e = Pe@D0
            r_ssh_o = Po@REG
            return np.concatenate([rg.ravel(), rp.ravel(), rm.ravel(),
                                   rbc_g, rbc_X, r_ssh_e, r_ssh_o])
        except FloatingPointError:
            return np.full(31*(self.Nz - 1) + 42, 30.0)

    # ---- sparsity for grouped finite differences ----
    def sparsity(self):
        from scipy.sparse import lil_matrix
        nr = 31*(self.Nz - 1) + 42
        Sp = lil_matrix((nr, self.nu), dtype=np.int8)
        Sp[:, :11] = 1                                   # Delta + xi0: dense columns
        for k in range(self.Nz - 1):
            r0 = 0
            cols = slice(11 + 31*k, 11 + 31*(k + 2))
            Sp[r0 + 11*k: r0 + 11*(k+1), cols] = 1                       # rg block k
            r1 = 11*(self.Nz - 1)
            Sp[r1 + 10*k: r1 + 10*(k+1), cols] = 1                       # rp block k
            r2 = r1 + 10*(self.Nz - 1)
            Sp[r2 + 10*k: r2 + 10*(k+1), cols] = 1                       # rm block k
        rb = 31*(self.Nz - 1)
        Sp[rb:rb+21, 11:11+31] = 1                                       # center BCs: node 0
        Sp[rb+21:, 11 + 31*(self.Nz-1):] = 1                             # SSH: last node
        return Sp

def seed_from_cylinder(rx, cyl, xi0_tau, Delta):
    """Interpolate the evolver cylinder onto the relaxation grid IN THE GAUGED FRAME
    (zeta_G = zeta_raw - xi0(tau)) and project to harmonics. The whole cylinder is first
    tau-ROTATED by the common shift that kills xi0's k=2 sine (the solver's gauge)."""
    Be, Bde, Bo, Bdo, Pe, Po = bases(Delta)
    zc = cyl['z']
    # gauge rotation angle from xi0's k=2 harmonic
    th = 2*np.pi*np.arange(len(xi0_tau))/len(xi0_tau)
    c2 = 2*np.mean(xi0_tau*np.cos(2*th)); s2 = 2*np.mean(xi0_tau*np.sin(2*th))
    ph2 = 0.5*np.arctan2(s2, c2)
    def rot(f):
        F_ = np.fft.rfft(f)
        kk = np.arange(len(F_))
        return np.fft.irfft(F_*np.exp(1j*kk*ph2), len(f))
    xi0_tau = rot(xi0_tau)
    cyl = dict(cyl)
    for key in ('X', 'Y', 'g', 'a'):
        cyl[key] = np.array([_ for _ in cyl[key]])
        cyl[key] = np.apply_along_axis(rot, 0, cyl[key])   # rotate along the tau axis
    u = np.zeros(rx.nu)
    u[0] = Delta
    xc = Pe@xi0_tau
    u[1] = xc[0]; u[2] = xc[1]; u[3:11] = xc[3:]
    F = np.zeros((rx.Nz, 31))
    for k, zG in enumerate(rx.z):
        gv = np.empty(M); xv = np.empty(M); yv = np.empty(M)
        for i in range(M):
            zr = zG + xi0_tau[i]
            gv[i] = np.interp(zr, zc, cyl['g'][i])
            xv[i] = np.interp(zr, zc, cyl['X'][i])
            yv[i] = np.interp(zr, zc, cyl['Y'][i])
        F[k, :11] = Pe@gv
        F[k, 11:21] = Po@(xv + yv)
        F[k, 21:31] = Po@(xv - yv)
    u[11:] = F.ravel()
    return u

def vacuum_control(Nz=40):
    rx = Relax(Nz=Nz)
    u = np.zeros(rx.nu); u[0] = 3.4453
    for k in range(rx.Nz):
        u[11 + 31*k] = 1.0                    # g = 1 (c0)
    r = rx.residual(u)
    print(f"[vac-relax] |r|max = {np.abs(r).max():.2e}  (flat space must satisfy everything exactly)")
    return np.abs(r).max()

if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "vac"
    if mode == "vac":
        vacuum_control()
