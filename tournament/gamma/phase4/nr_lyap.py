# nr_lyap.py — fluid-beta eigenvalue by HKA's SECOND method: the LYAPUNOV TIME-EVOLUTION
# (rflanl.tex sec V.G + app-lyapunov), on the PROVEN operator (D-031, nr_rederive diff=0).
#
# WHY THIS WORKS WHERE THE SHOOT WALLED: the s-evolution has NO sonic singularity — M_s,fl
# ([[As,Bs],[Cs,Ds]]) is invertible everywhere (det = g[1-(g-1)V^2]/(1-V^2) > 0); M_x,fl is
# only MULTIPLIED, never inverted. Eigenmodes emerge in DESCENDING Re kappa: the relevant
# mode (kappa=2.81055255) dominates at late s; the gauge mode sits at EXACTLY kappa=1 in the
# center gauge Nbar_p(-inf)=0 (sec V.G) — a built-in analytic control.
#
# STRUCTURE (eq:EOM-var, verified coefficients hka_pert_hka99 = nr_rederive):
#   metric rows (NO d_s -> spatial constraints, slaved each evaluation):
#       d_x Abar = G1*Abar + G3*obar + G4*V        BC (center-regular): Abar -> 0
#       d_x Nbar = H1*Abar + H3*obar               BC (gauge, sec V.G): Nbar -> 0 at center
#   fluid rows (the dynamics), w = (obar, V):
#       M_s d_s w = (E.Psi, F.Psi) - M_x d_x w
#       => d_s w = S.(E.Psi, F.Psi) - Axi * xi * d_xi w,   S = M_s^{-1}, Axi = S.M_x
#
# STABILIZING COORDINATE (HKA's "appropriate coordinate transformations"): evolve on a
# UNIFORM grid in xi = e^x (so d_x = xi d_xi). In x the advection speed blows up ~ e^{-x}
# toward the center (N ~ 1/xi) making CFL catastrophic; in xi the speeds are BOUNDED
# (xi*speed_x -> O(1) at the center, -> 0 incoming at the sonic point = the free-BC fact).
#
# Domain xi in (0, xi_s] (xi_s = e^{x_s} ~ 0.866). Right edge = SONIC POINT, free BC
# (one-sided interior differences — no information enters from outside, sec V.G).
# Left edge = center, regularity BC (bounded modes -> const; zero-gradient ghosts).
# The two irregular center modes (1/xi, 1/xi^2 = e^{-x}, e^{-2x}) are excluded by the BC.
#
# Scheme: method of lines, RK4, 4th-order central differences + Kreiss-Oliger dissipation
# (house standard; the paper used Lax-Wendroff — both are stable dissipative 2nd+-order
# schemes; MOL-RK4+KO lets the same code do high-accuracy convergence runs).
# Benettin/Gram-Schmidt n-frame (app-lyapunov): kappa_i = mean growth of the i-th
# orthogonalized direction. kappa_1 = relevant -> beta = 1/kappa_1. kappa_2 = gauge = 1.
#
# Controls:
#   1. `python nr_lyap.py gauge`  — initialize with the EXACT gauge mode (obar_x, V_x);
#      growth must be exactly 1 and the shape invariant (validates advection+source+slaving).
#   2. default run — 3-frame Benettin from generic data: expect (2.81055, 1.0000, <1).
#   3. refinement — kappa stationary under N doubling (G-CONVERGE), unique relevant (G-UNIQUE).
import numpy as np, sys, time
import hka_pert_core as PC, hka_pert_hka99 as H99
PC.Lnum = H99.Lnum
import hka_beta4 as B

B.bg(); B.bg_path()
XS = B.bg()['xs']

# ---------------------------------------------------------------- setup on the xi grid
class Grid:
    def __init__(self, N=800, xi_edge_frac=1.0-1e-6, bg_state=None, xs=None):
        """bg_state(x) -> (A,N,om,V,obar_x,V_x) and xs (sonic x) override the default
        (hka_beta4 = the FRIEDMANN background, kept for the record); pass the true EC
        background via nr_ec2.build_ec()."""
        self.N = N
        if bg_state is None: bg_state = B.bg_state
        if xs is None: xs = XS
        xi_s = np.exp(xs)
        self.xi = np.linspace(xi_s/N*0.5, xi_s*xi_edge_frac, N)   # first pt ~ half-cell off 0
        self.dxi = self.xi[1] - self.xi[0]
        x = np.log(self.xi)
        # background + operator coefficients per grid point
        cof = {k: np.zeros(N) for k in
               ['G1','G3','G4','H1','H3','E1','E2','E3','E4','F1','F2','F3','F4']}
        S = np.zeros((N,2,2)); AXI = np.zeros((N,2,2))
        obx = np.zeros(N); Vx = np.zeros(N)
        for j in range(N):
            fld = bg_state(x[j])
            c = H99.coeffs(fld)
            for k in cof: cof[k][j] = c[k]
            Ms = np.array([[c['As'], c['Bs']], [c['Cs'], c['Ds']]])
            Mx = np.array([[c['Ax'], c['Bx']], [c['Cx'], c['Dx']]])
            Si = np.linalg.inv(Ms)
            S[j] = Si; AXI[j] = self.xi[j]*Si.dot(Mx)             # xi-space advection matrix
            obx[j] = float(np.real(fld[4])); Vx[j] = float(np.real(fld[5]))
        self.cof = cof; self.S = S; self.AXI = AXI
        self.obx = obx; self.Vx = Vx
        # the FIFTH equation (rflanl.tex l.5090, "We also have" — the linearized momentum
        # constraint, an EVOLUTION equation for Abar):
        #     d_s Abar = - d_x Abar - c1*(Nbar + obar) - c2*V
        # (its eigen-form + row 0 of eq:EOM-var is exactly how the paper derives eq:alg-PP:
        #  cN = c1, com = G3 + c1, cV = G4 + c2 — verified to match identity_coeffs.)
        g43 = 4.0/3.0
        Abg = np.zeros(N); Nbg = np.zeros(N); ombg = np.zeros(N); Vbg = np.zeros(N)
        for j in range(N):
            A_, N_, om_, V_, _, _ = [float(np.real(z)) for z in bg_state(x[j])]
            Abg[j], Nbg[j], ombg[j], Vbg[j] = A_, N_, om_, V_
        oV2 = 1.0 - Vbg*Vbg
        self.c1 = 2.0*g43*Nbg*Vbg*ombg/oV2
        self.c2 = 2.0*g43*Nbg*ombg*(1.0 + Vbg*Vbg)/oV2**2
        # background log-slopes (for the exact gauge mode in B): Abar_x (Eq1), Nbar_x (Eq2)
        self.Abx = 1.0 - Abg + 2.0*ombg/oV2*(1.0 + Vbg*Vbg/3.0)
        self.Nbx = -2.0 + Abg - 2.0*ombg/3.0
        # integrating factor for the Abar slaving: mu = exp(int G1/xi dxi), mu[0]=1
        g1oxi = cof['G1']/self.xi
        self.mu = np.exp(np.concatenate([[0.0], np.cumsum(0.5*(g1oxi[1:]+g1oxi[:-1])*self.dxi)]))
        self.spd = max(abs(np.linalg.eigvals(AXI[j]).real).max() for j in range(0, N, max(1, N//64)))
        self.spd = max(self.spd, self.xi[-1])     # B: Abar advects at speed xi (unit speed in x)

    def cumtrapz0(self, f):
        """cumulative trapezoid with F[0]=0 on the uniform xi grid."""
        return np.concatenate([[0.0], np.cumsum(0.5*(f[1:]+f[:-1])*self.dxi)])

    def slave(self, ob, V):
        """(Abar, Nbar) from the spatial constraints. Abar(center)=0 (regular), Nbar(center)=0 (gauge)."""
        q = (self.cof['G3']*ob + self.cof['G4']*V)/self.xi
        Ab = self.mu*self.cumtrapz0(q/self.mu)
        Nb = self.cumtrapz0((self.cof['H1']*Ab + self.cof['H3']*ob)/self.xi)
        return Ab, Nb

    def dxi4(self, f):
        """4th-order first derivative (interior); one-sided at both edges (stable base scheme).
        Sonic edge = free BC (sec V.G). Center edge: the physical regularity datum for the one
        incoming characteristic is enforced separately (V-penalty in rhs), not via ghosts —
        parity ghosts at the staggered corner proved violently unstable (measured: kappa~+13 pair)."""
        N = self.N; d = np.empty(N); h = self.dxi
        d[2:-2] = (f[:-4] - 8*f[1:-3] + 8*f[3:-1] - f[4:])/(12*h)
        d[0] = (-3*f[0] + 4*f[1] - f[2])/(2*h)
        d[1] = (f[2] - f[0])/(2*h)
        d[-2] = (f[-1] - f[-3])/(2*h)
        d[-1] = (3*f[-1] - 4*f[-2] + f[-3])/(2*h)
        return d

    def rhs(self, ob, V, pen=0.0):
        Ab, Nb = self.slave(ob, V)
        EP = self.cof['E1']*Ab + self.cof['E2']*Nb + self.cof['E3']*ob + self.cof['E4']*V
        FP = self.cof['F1']*Ab + self.cof['F2']*Nb + self.cof['F3']*ob + self.cof['F4']*V
        dob = self.dxi4(ob); dV = self.dxi4(V)
        S = self.S; A = self.AXI
        dsob = S[:,0,0]*EP + S[:,0,1]*FP - (A[:,0,0]*dob + A[:,0,1]*dV)
        dsV  = S[:,1,0]*EP + S[:,1,1]*FP - (A[:,1,0]*dob + A[:,1,1]*dV)
        if pen > 0.0:
            # center regularity (eq:PPasmp1): bounded modes have V_p -> 0; relax the corner value
            # toward it (SAT-style penalty for the incoming characteristic; gentle, rate pen/s).
            dsV[0] -= pen*V[0]
        return dsob, dsV

    def ko(self, f, eps=0.5):
        """Kreiss-Oliger dissipation: 4th-order interior, 2nd-order (1,-2,1) at the edge points
        (uniform positive sign — a sign-flipped edge term is ANTI-dissipative, caught in review)."""
        g = f.copy()
        g[2:-2] -= (eps/16.0)*(f[:-4] - 4*f[1:-3] + 6*f[2:-2] - 4*f[3:-1] + f[4:])
        g[1]  -= (eps/4.0)*(f[0] - 2*f[1] + f[2])
        g[-2] -= (eps/4.0)*(f[-1] - 2*f[-2] + f[-3])
        return g

    def step(self, ob, V, ds, pen=0.0):
        k1o, k1v = self.rhs(ob, V, pen)
        k2o, k2v = self.rhs(ob + 0.5*ds*k1o, V + 0.5*ds*k1v, pen)
        k3o, k3v = self.rhs(ob + 0.5*ds*k2o, V + 0.5*ds*k2v, pen)
        k4o, k4v = self.rhs(ob + ds*k3o, V + ds*k3v, pen)
        ob = ob + (ds/6)*(k1o + 2*k2o + 2*k3o + k4o)
        V = V + (ds/6)*(k1v + 2*k2v + 2*k3v + k4v)
        return self.ko(ob), self.ko(V)

    # ---------------- formulation B: the PHYSICALLY COMPLETE system ----------------
    # Dynamical (Abar, obar, V); Abar evolves by the FIFTH equation (momentum constraint,
    # rflanl.tex l.5090); Nbar slaved by row 1 (gauge Nbar(center)=0, sec V.G); row 0 of
    # eq:EOM-var becomes the MONITORED constraint (preserved on physical solutions).
    def slaveN(self, Ab, ob):
        return self.cumtrapz0((self.cof['H1']*Ab + self.cof['H3']*ob)/self.xi)

    def rhsB(self, Ab, ob, V, pen=0.0):
        Nb = self.slaveN(Ab, ob)
        EP = self.cof['E1']*Ab + self.cof['E2']*Nb + self.cof['E3']*ob + self.cof['E4']*V
        FP = self.cof['F1']*Ab + self.cof['F2']*Nb + self.cof['F3']*ob + self.cof['F4']*V
        dob = self.dxi4(ob); dV = self.dxi4(V); dAb = self.dxi4(Ab)
        S = self.S; A = self.AXI
        dsob = S[:,0,0]*EP + S[:,0,1]*FP - (A[:,0,0]*dob + A[:,0,1]*dV)
        dsV  = S[:,1,0]*EP + S[:,1,1]*FP - (A[:,1,0]*dob + A[:,1,1]*dV)
        dsAb = -self.xi*dAb - self.c1*(Nb + ob) - self.c2*V
        if pen > 0.0:
            dsV[0]  -= pen*V[0]     # center regularity: V_p -> 0 (eq:PPasmp1)
            dsAb[0] -= pen*Ab[0]    # center regularity: Abar_p -> 0 (eq:PPasmp1)
        return dsAb, dsob, dsV

    def stepB(self, Ab, ob, V, ds, pen=0.0):
        k1a, k1o, k1v = self.rhsB(Ab, ob, V, pen)
        k2a, k2o, k2v = self.rhsB(Ab + 0.5*ds*k1a, ob + 0.5*ds*k1o, V + 0.5*ds*k1v, pen)
        k3a, k3o, k3v = self.rhsB(Ab + 0.5*ds*k2a, ob + 0.5*ds*k2o, V + 0.5*ds*k2v, pen)
        k4a, k4o, k4v = self.rhsB(Ab + ds*k3a, ob + ds*k3o, V + ds*k3v, pen)
        Ab = Ab + (ds/6)*(k1a + 2*k2a + 2*k3a + k4a)
        ob = ob + (ds/6)*(k1o + 2*k2o + 2*k3o + k4o)
        V = V + (ds/6)*(k1v + 2*k2v + 2*k3v + k4v)
        return self.ko(Ab), self.ko(ob), self.ko(V)

    def row0_res(self, Ab, ob, V):
        """Row-0 constraint residual (physicality monitor): xi*d_xi Abar - (G1 Ab + G3 ob + G4 V)."""
        R0 = self.xi*self.dxi4(Ab) - (self.cof['G1']*Ab + self.cof['G3']*ob + self.cof['G4']*V)
        den = max(np.abs(Ab).max(), np.abs(ob).max(), np.abs(V).max(), 1e-30)
        return np.abs(R0[2:-2]).max()/den

def nrm(g, ob, V):
    return np.sqrt(np.sum(ob*ob + V*V)*g.dxi)

# ---------------------------------------------------------------- controls & runs
def run_gauge(N=800, s_end=3.0, cfl=0.4):
    """CONTROL: init with the exact gauge mode (obar_x, V_x) -> growth EXACTLY 1, shape frozen."""
    g = Grid(N)
    ds = cfl*g.dxi/g.spd
    ob, V = g.obx.copy(), g.Vx.copy()
    n0 = nrm(g, ob, V); ob /= n0; V /= n0
    shape0 = np.concatenate([ob, V])
    # bonus static check: slaved (Abar,Nbar) must reproduce the kbar=1 gauge member
    Ab, Nb = g.slave(g.obx, g.Vx)
    x = np.log(g.xi)
    Abx = np.zeros(g.N); Nbx = np.zeros(g.N)
    for j in range(g.N):
        A_, N_, om_, V_, _, _ = [float(np.real(z)) for z in B.bg_state(x[j])]
        Abx[j] = 1.0 - A_ + 2.0*om_/(1.0 - V_*V_)*(1.0 + V_*V_/3.0)
        Nbx[j] = -2.0 + A_ - 2.0*om_/3.0
    eA = np.max(np.abs(Ab - Abx)); eN = np.max(np.abs(Nb - (Nbx + 1.0)))
    print(f"[gauge-slave] max|Abar_slaved - Abar_x| = {eA:.2e}   max|Nbar_slaved - (Nbar_x+1)| = {eN:.2e}")
    logn, s = 0.0, 0.0
    t0 = time.time(); nstep = int(round(s_end/ds))
    for i in range(nstep):
        ob, V = g.step(ob, V, ds, pen=2.0); s += ds
        nn = nrm(g, ob, V); logn += np.log(nn); ob /= nn; V /= nn
    shape = np.concatenate([ob, V]); shape /= np.linalg.norm(shape)
    sh_dev = min(np.linalg.norm(shape - shape0/np.linalg.norm(shape0)),
                 np.linalg.norm(shape + shape0/np.linalg.norm(shape0)))
    kap = logn/s
    print(f"[gauge] N={N} ds={ds:.2e} s={s:.2f}: growth kappa = {kap:.6f} (exact 1)   "
          f"shape drift = {sh_dev:.2e}   ({time.time()-t0:.1f}s)")
    return kap

def run_benettin(N=800, s_end=14.0, s_burn=4.0, nvec=3, cfl=0.4, seed=1, gs_every=10,
                 deflate_gauge=False, quiet=False):
    """Main run: nvec-frame Benettin/Gram-Schmidt from generic smooth data -> kappa_1..kappa_n.
    deflate_gauge: project the EXACT (analytically known) gauge mode (obar_x, V_x) out of every
    vector at each GS -> kappa_1 is then the largest NON-gauge exponent (the relevant-mode reader)."""
    g = Grid(N)
    ds = cfl*g.dxi/g.spd
    rng = np.random.default_rng(seed)
    gvec = np.concatenate([g.obx, g.Vx]); gvec = gvec/np.linalg.norm(gvec)
    W = []
    xiN = g.xi/g.xi[-1]
    for k in range(nvec):
        ob = np.sin((k+1)*np.pi*xiN)*np.exp(-8*(xiN-0.4-0.1*k)**2) + 0.3*rng.standard_normal(g.N)*xiN*(1-xiN)
        V = np.cos((k+2)*np.pi*xiN)*np.exp(-8*(xiN-0.5+0.1*k)**2) + 0.3*rng.standard_normal(g.N)*xiN*(1-xiN)
        W.append(np.concatenate([ob, V]))
    def gs(Wm):
        # NOTE: the SAME norm must be used for the log and the renormalization (a mismatched
        # weight injects a constant bias per GS event — caught by the first N=400 run).
        out = []; logs = []
        for k in range(len(Wm)):
            v = Wm[k].copy()
            if deflate_gauge: v -= (gvec @ v)*gvec
            for q in out: v -= (q @ v)*q
            nv = np.linalg.norm(v)
            logs.append(np.log(nv)); out.append(v/nv)
        return out, logs
    W, _ = gs(W)
    t0 = time.time(); nstep = int(round(s_end/ds)); s = 0.0
    acc = np.zeros(nvec); acc_s = 0.0
    win = np.zeros(nvec); win_s = 0.0; next_report = s_burn
    for i in range(nstep):
        for k in range(nvec):
            ob, V = W[k][:g.N], W[k][g.N:]
            ob, V = g.step(ob, V, ds, pen=2.0)
            W[k] = np.concatenate([ob, V])
        s += ds
        if (i+1) % gs_every == 0:
            W, logs = gs(W)
            if s > s_burn:
                acc += np.array(logs); acc_s += gs_every*ds
            win += np.array(logs); win_s += gs_every*ds
            if not quiet and s >= next_report:
                gov = abs(gvec @ W[0])
                print(f"    s={s:6.2f}  running kappa={['%+.4f' % k for k in win/win_s]}  |<w1,gauge>|={gov:.3f}", flush=True)
                win[:] = 0.0; win_s = 0.0; next_report += max(1.0, s_end/10)
    kaps = acc/acc_s
    if not quiet:
        gov = abs(gvec @ W[0])
        print(f"[benettin] N={N} ds={ds:.2e} s_end={s_end} burn={s_burn} deflate_gauge={deflate_gauge} "
              f"final |<w1,gauge>|={gov:.3f} ({time.time()-t0:.1f}s):")
        for k in range(nvec):
            extra = f"   beta = 1/kappa = {1.0/kaps[k]:.8f}  (ref 0.35580192, err {abs(1.0/kaps[k]-0.35580192):.2e})" if k == 0 else ""
            print(f"  kappa_{k+1} = {kaps[k]:+.6f}{extra}")
    return kaps

def run_gaugeB(N=400, s_end=3.0, cfl=0.4):
    """B CONTROL: init with the exact gauge mode (Abar_x, obar_x, V_x) -> growth EXACTLY 1."""
    g = Grid(N)
    ds = cfl*g.dxi/g.spd
    Ab, ob, V = g.Abx.copy(), g.obx.copy(), g.Vx.copy()
    n0 = np.sqrt(np.sum(Ab*Ab + ob*ob + V*V)); Ab /= n0; ob /= n0; V /= n0
    shape0 = np.concatenate([Ab, ob, V])
    Nb = g.slaveN(g.Abx, g.obx)
    eN = np.max(np.abs(Nb - (g.Nbx + 1.0)))
    print(f"[gaugeB-slave] max|Nbar_slaved - (Nbar_x+1)| = {eN:.2e}   "
          f"row0-residual(gauge) = {g.row0_res(g.Abx, g.obx, g.Vx):.2e}")
    logn, s = 0.0, 0.0
    t0 = time.time(); nstep = int(round(s_end/ds))
    for i in range(nstep):
        Ab, ob, V = g.stepB(Ab, ob, V, ds, pen=2.0); s += ds
        nn = np.sqrt(np.sum(Ab*Ab + ob*ob + V*V)); logn += np.log(nn); Ab /= nn; ob /= nn; V /= nn
    shape = np.concatenate([Ab, ob, V]); shape /= np.linalg.norm(shape)
    sh_dev = min(np.linalg.norm(shape - shape0/np.linalg.norm(shape0)),
                 np.linalg.norm(shape + shape0/np.linalg.norm(shape0)))
    print(f"[gaugeB] N={N} ds={ds:.2e} s={s:.2f}: growth kappa = {logn/s:.6f} (exact 1)   "
          f"shape drift = {sh_dev:.2e}   row0-res(final) = {g.row0_res(Ab, ob, V):.2e}   ({time.time()-t0:.1f}s)")
    return logn/s

def run_benettinB(N=400, s_end=14.0, s_burn=4.0, nvec=3, cfl=0.4, seed=1, gs_every=10,
                  deflate_gauge=False, quiet=False):
    """B main run: nvec-frame Benettin on the COMPLETE system (Abar dynamical via the fifth eq).
    Initial data starts ON the row-0 constraint surface (Abar slaved once at s=0); row-0 residual
    of w1 is monitored — physical modes keep it small (Bianchi), junk does not."""
    g = Grid(N)
    ds = cfl*g.dxi/g.spd
    rng = np.random.default_rng(seed)
    gv = np.concatenate([g.Abx, g.obx, g.Vx]); gv = gv/np.linalg.norm(gv)
    W = []
    xiN = g.xi/g.xi[-1]
    for k in range(nvec):
        ob = np.sin((k+1)*np.pi*xiN)*np.exp(-8*(xiN-0.4-0.1*k)**2) + 0.3*rng.standard_normal(g.N)*xiN*(1-xiN)
        V = np.cos((k+2)*np.pi*xiN)*np.exp(-8*(xiN-0.5+0.1*k)**2) + 0.3*rng.standard_normal(g.N)*xiN*(1-xiN)
        Ab0, _ = g.slave(ob, V)                 # constraint-consistent start (row 0 satisfied)
        W.append(np.concatenate([Ab0, ob, V]))
    def gs(Wm):
        out = []; logs = []
        for k in range(len(Wm)):
            v = Wm[k].copy()
            if deflate_gauge: v -= (gv @ v)*gv
            for q in out: v -= (q @ v)*q
            nv = np.linalg.norm(v)
            logs.append(np.log(nv)); out.append(v/nv)
        return out, logs
    W, _ = gs(W)
    t0 = time.time(); nstep = int(round(s_end/ds)); s = 0.0
    acc = np.zeros(nvec); acc_s = 0.0
    win = np.zeros(nvec); win_s = 0.0; next_report = s_burn
    for i in range(nstep):
        for k in range(nvec):
            Ab, ob, V = W[k][:g.N], W[k][g.N:2*g.N], W[k][2*g.N:]
            Ab, ob, V = g.stepB(Ab, ob, V, ds, pen=2.0)
            W[k] = np.concatenate([Ab, ob, V])
        s += ds
        if (i+1) % gs_every == 0:
            W, logs = gs(W)
            if s > s_burn:
                acc += np.array(logs); acc_s += gs_every*ds
            win += np.array(logs); win_s += gs_every*ds
            if not quiet and s >= next_report:
                gov = abs(gv @ W[0])
                r0 = g.row0_res(W[0][:g.N], W[0][g.N:2*g.N], W[0][2*g.N:])
                print(f"    s={s:6.2f}  running kappa={['%+.4f' % k for k in win/win_s]}  "
                      f"|<w1,gauge>|={gov:.3f}  row0-res(w1)={r0:.2e}", flush=True)
                win[:] = 0.0; win_s = 0.0; next_report += max(1.0, s_end/10)
    kaps = acc/acc_s
    if not quiet:
        gov = abs(gv @ W[0])
        r0 = g.row0_res(W[0][:g.N], W[0][g.N:2*g.N], W[0][2*g.N:])
        print(f"[benettinB] N={N} ds={ds:.2e} s_end={s_end} burn={s_burn} deflate_gauge={deflate_gauge} "
              f"final |<w1,gauge>|={gov:.3f} row0-res(w1)={r0:.2e} ({time.time()-t0:.1f}s):")
        for k in range(nvec):
            extra = f"   beta = 1/kappa = {1.0/kaps[k]:.8f}  (ref 0.35580192, err {abs(1.0/kaps[k]-0.35580192):.2e})" if k == 0 else ""
            print(f"  kappa_{k+1} = {kaps[k]:+.6f}{extra}")
    return kaps

def run_spec(N=200, cfl=0.4, nshow=14, pen=2.0, bg=None):
    """DECISIVE: full spectrum of the DISCRETE one-step map of formulation B (KO+penalty included
    — exactly the dynamics Benettin samples). kappa = log(lambda_step)/ds. If 2.81 is in the
    discrete dynamics it MUST appear here; if not, A/B the BC/domain variants at operator level.
    bg = (bg_state, xs) tuple to run on an alternative background (nr_ec2 = the true EC)."""
    g = Grid(N, bg_state=bg[0] if bg else None, xs=bg[1] if bg else None)
    ds = cfl*g.dxi/g.spd
    n3 = 3*N
    M = np.zeros((n3, n3))
    t0 = time.time()
    for j in range(n3):
        u = np.zeros(n3); u[j] = 1.0
        Ab, ob, V = g.stepB(u[:N], u[N:2*N], u[2*N:], ds, pen=pen)
        M[:, j] = np.concatenate([Ab, ob, V])
    lam = np.linalg.eigvals(M)
    kap = np.log(lam.astype(complex))/ds
    order = np.argsort(-kap.real)
    print(f"[spec] N={N} ds={ds:.2e} pen={pen} (one-step map, {time.time()-t0:.1f}s) — top {nshow} kappa:")
    gvec = np.concatenate([g.Abx, g.obx, g.Vx]); gvec /= np.linalg.norm(gvec)
    for i in order[:nshow]:
        tag = "  <== RELEVANT?" if abs(kap[i].real - 2.81055) < 0.3 and abs(kap[i].imag) < 0.5 else ""
        tag = tag or ("  <== gauge (exact 1)?" if abs(kap[i].real - 1.0) < 0.05 and abs(kap[i].imag) < 0.5 else "")
        print(f"  kappa = {kap[i].real:+9.5f} {kap[i].imag:+9.5f}i{tag}")
    return kap

if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "run"
    N = int(sys.argv[2]) if len(sys.argv) > 2 else 400
    if mode == "spec":
        run_spec(N if len(sys.argv) > 2 else 200)
    elif mode == "gauge":
        run_gaugeB(N)
    elif mode == "gaugeA":
        run_gauge(N)
    elif mode == "conv":
        for NN in (400, 800, 1600):
            run_benettinB(N=NN)
    elif mode == "nogauge":
        run_benettinB(N=N, deflate_gauge=True)
    elif mode == "runA":
        run_benettin(N=N)
    elif mode == "nogaugeA":
        run_benettin(N=N, deflate_gauge=True)
    else:
        run_benettinB(N=N)
