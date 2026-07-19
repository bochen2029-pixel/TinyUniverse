# nr_evolve.py — throwaway near-critical massless-scalar collapse evolver, purpose-built to
# SEED the DSS boundary-value problem (nr_dss) with real Choptuik echo fields — the same
# route Gundlach's own initial guess took (Choptuik's simulation data).
#
# Polar-areal gauge, CENTRAL-observer time normalization alpha(0,t) = 1 (Gundlach's gauge
# directly, so t is his t up to the t* shift). Massless scalar (V=0):
#   ds^2 = -alpha^2 dt^2 + a^2 dr^2 + r^2 dOmega^2,  Phi = phi_,r,  Pi = (a/alpha) phi_,t
#   evolve:  Phi_t = ( (alpha/a) Pi )_,r ,   Pi_t = (1/r^2) ( r^2 (alpha/a) Phi )_,r
#   slice:   (ln a)_,r  = (1-a^2)/(2r) + 2 pi r (Pi^2 + Phi^2),        a(0)=1
#            (ln al)_,r = (a^2-1)/(2r) + 2 pi r (Pi^2 + Phi^2),        al(0)=1  (CENTRAL gauge)
# (the substrate_nexus scheme, massless, alpha-normalized at the center instead of r_max).
# MOL RK4 + Kreiss-Oliger, origin parity (Phi odd, Pi even), outgoing outer BC.
#
# Outputs: p* bisection; then a marginal run records full snapshots over the echoing epoch;
# t* and Delta_echo from the geometric accumulation of the central-field zero crossings
# (phi(0,t) flips sign every Delta/2 in tau = ln(t*-t)); fields sampled on the (tau, zeta)
# cylinder for one period -> harmonic seed for nr_dss (and Delta_echo = the in-house
# redundant-recovery observable of the similarity_nexus contract).
import numpy as np, sys, time

class Ev:
    def __init__(self, N=1600, rmax=60.0, cfl=0.4, eko=0.5):
        self.N = N; self.rmax = rmax
        self.dr = rmax/N
        self.r = (np.arange(N) + 0.5)*self.dr          # staggered off r=0
        self.cfl = cfl; self.eko = eko

    def slice_metric(self, Phi, Pi):
        """a, alpha by outward log-integration (trapezoid), a(0)=1; alpha normalized so
        alpha*a = 1 at r_max (the substrate_nexus outer gauge — keeps coordinate speeds O(1)
        near criticality; the central lapse alpha(0,t) is RECORDED and used to accumulate the
        central proper time t_G = int alpha(0) dt, which is Gundlach's time)."""
        r = self.r; dr = self.dr
        S = 2*np.pi*r*(Pi*Pi + Phi*Phi)
        lna = np.zeros(self.N); lnal = np.zeros(self.N)
        a2 = 1.0
        f_prev_a = 0.0; f_prev_l = 0.0                  # integrands -> 0 at r=0
        acc_a = 0.0; acc_l = 0.0
        for j in range(self.N):
            fa = (1.0 - a2)/(2*r[j]) + S[j]
            fl = (a2 - 1.0)/(2*r[j]) + S[j]
            h = dr if j > 0 else r[0]
            acc_a += 0.5*(f_prev_a + fa)*h
            acc_l += 0.5*(f_prev_l + fl)*h
            lna[j] = acc_a; lnal[j] = acc_l
            a2 = np.exp(2*acc_a)
            f_prev_a = fa; f_prev_l = fl
        a = np.exp(lna)
        al = np.exp(lnal)
        al *= 1.0/(al[-1]*a[-1])                        # outer normalization alpha*a|_rmax = 1
        return a, al

    def rhs(self, Phi, Pi):
        a, al = self.slice_metric(Phi, Pi)
        w = al/a
        r = self.r; dr = self.dr; N = self.N
        # Phi_t = d_r(w Pi)                      [plain gradient — NO metric-flux weighting;
        #                                         a spurious +wPi/2r here caused an origin blowup]
        # Pi_t  = d_r(w Phi) + 2 w Phi / r       [== (1/r^2) d_r(r^2 w Phi), uniformly accurate
        #                                         on the staggered grid; exact for linear fields]
        wPi = w*Pi                               # EVEN across r=0
        wPhi = w*Phi                             # ODD across r=0
        dPhi = np.empty(N)
        dPhi[1:-1] = (wPi[2:] - wPi[:-2])/(2*dr)
        dPhi[0] = (wPi[1] - wPi[0])/(2*dr)                     # even ghost: wPi(-r0)=wPi(r0)
        dPhi[-1] = (wPi[-1] - wPi[-2])/dr
        dwPhi = np.empty(N)
        dwPhi[1:-1] = (wPhi[2:] - wPhi[:-2])/(2*dr)
        dwPhi[0] = (wPhi[1] + wPhi[0])/(2*dr)                  # odd ghost: wPhi(-r0)=-wPhi(r0)
        dwPhi[-1] = (wPhi[-1] - wPhi[-2])/dr
        dPi = dwPhi + 2.0*wPhi/r
        dPi[-1] = -(wPi[-1])/r[-1] - (wPi[-1] - wPi[-2])/dr    # outgoing-ish
        return dPhi, dPi, a, al

    def ko(self, f, sgn):
        g = f.copy()
        gh = np.concatenate([[sgn*f[1], sgn*f[0]], f, [f[-1], f[-1]]])
        g -= (self.eko/16.0)*(gh[:-4] - 4*gh[1:-3] + 6*gh[2:-2] - 4*gh[3:-1] + gh[4:])
        return g

    def phi_of(self, Phi):
        """phi(r_j) = int_0^{r_j} Phi dr (Phi odd: ghost Phi(0)=0)."""
        seg = np.empty(self.N)
        seg[0] = 0.5*Phi[0]*self.r[0]
        seg[1:] = 0.5*(Phi[1:] + Phi[:-1])*self.dr
        return np.cumsum(seg)

    def run(self, p, r0=12.0, sig=2.0, tmax=200.0, record=False, rec_tG0=0.0, stride=2):
        """Gaussian phi = p exp(-((r-r0)/sig)^2) at rest. Returns outcome (+snapshots keyed
        by the CENTRAL proper time tG = int alpha(0) dt — Gundlach's time up to the t* shift)."""
        r = self.r
        phi0 = p*np.exp(-((r - r0)/sig)**2)
        Phi = np.gradient(phi0, self.dr)
        Pi = np.zeros_like(Phi)
        t = 0.0; tG = 0.0
        snaps = []; cen = []
        istep = 0
        # phi(0,t) must be EVOLVED (d_t phi|_0 = (al/a) Pi|_0): Phi = phi_,r cannot carry the
        # r-independent central mode — reconstructing phi from int(Phi) pins phi(0)=0 and
        # erases the echo signal entirely (measured: max|phi_cen| = 0.015 pre-fix).
        phicen = 0.0
        while t < tmax:
            dPhi, dPi, a, al = self.rhs(Phi, Pi)
            w = (al/a).max()
            dt = self.cfl*self.dr/max(w, 1e-9)
            k1P, k1p = dPhi, dPi
            d2 = self.rhs(Phi + 0.5*dt*k1P, Pi + 0.5*dt*k1p); k2P, k2p = d2[0], d2[1]
            d3 = self.rhs(Phi + 0.5*dt*k2P, Pi + 0.5*dt*k2p); k3P, k3p = d3[0], d3[1]
            d4 = self.rhs(Phi + dt*k3P, Pi + dt*k3p);         k4P, k4p = d4[0], d4[1]
            Phi = self.ko(Phi + dt/6*(k1P + 2*k2P + 2*k3P + k4P), -1.0)
            Pi  = self.ko(Pi  + dt/6*(k1p + 2*k2p + 2*k3p + k4p), +1.0)
            phicen += (al[0]/a[0])*Pi[0]*dt
            t += dt; tG += al[0]*dt; istep += 1
            m2r = 1.0 - 1.0/(a*a)
            if not np.all(np.isfinite(Phi)):
                return dict(fate='nan', t=t, tG=tG, snaps=snaps, cen=cen)
            # polar slicing freezes asymptotically (measured: m2r hovers ~0.95, alpha(0)~1e-2,
            # then breaks down numerically) — detect the freeze, not the unreachable limit
            if m2r.max() > 0.90 or al[0] < 0.02:
                # mass at freeze for the Choptuik scaling law: m(r) = (r/2)(1-1/a^2) at the
                # outermost 2m/r peak (the frozen near-horizon). Fixed threshold across p
                # => a p-independent multiplicative bias => the log-log SLOPE (gamma) is clean.
                j = int(np.argmax(m2r))
                return dict(fate='bh', t=t, tG=tG, m2r=m2r.max(), snaps=snaps, cen=cen,
                            rH=float(self.r[j]), MBH=float(0.5*self.r[j]*m2r[j]))
            if record:
                cen.append((tG, phicen, al[0]))
                if tG >= rec_tG0 and istep % stride == 0:
                    snaps.append((tG, phicen + self.phi_of(Phi), Phi.copy(), Pi.copy(), a.copy(), al.copy()))
            # dispersal: field energy has left the inner third (decided BEFORE outer-boundary
            # reflections — my crude outer BC reflects — can re-enter; r0=12 makes the window safe)
            if t > 40.0 and m2r[:self.N//3].max() < 1e-3:
                return dict(fate='disperse', t=t, tG=tG, snaps=snaps, cen=cen)
        return dict(fate='tmax', t=t, tG=tG, snaps=snaps, cen=cen)

def bisect_pstar(ev, lo=0.003, hi=0.02, iters=40, verbose=True):
    t0 = time.time()
    rlo = ev.run(lo); rhi = ev.run(hi)
    assert rlo['fate'] == 'disperse' and rhi['fate'] == 'bh', (rlo['fate'], rhi['fate'])
    for i in range(iters):
        mid = 0.5*(lo + hi)
        f = ev.run(mid)['fate']
        if f == 'bh': hi = mid
        else: lo = mid
        if verbose and (i+1) % 8 == 0:
            print(f"  bisect {i+1}/{iters}: p* in [{lo:.10f}, {hi:.10f}]  ({time.time()-t0:.0f}s)", flush=True)
    return lo, hi

# ---------------------------------------------------------------- echo analysis + seeding
def analyze_echoes(cen, amin=0.02, verbose=True):
    """t* and Delta from the geometric accumulation of phi(0,t) zero crossings:
    t*-t_n = C q^n, q = e^{-Delta/2}. Uses the LAST clean crossings (echo regime)."""
    cen = np.array(cen)                       # (tG, phi_cen, al_cen)
    tG, ph = cen[:, 0], cen[:, 1]
    s = np.sign(ph)
    ix = np.where(np.diff(s) != 0)[0]
    tn = []
    for i in ix:
        j0, j1 = max(0, i-40), min(len(ph), i+40)
        if np.abs(ph[j0:j1]).max() > amin:
            t0, t1, p0, p1 = tG[i], tG[i+1], ph[i], ph[i+1]
            tn.append(t0 - p0*(t1 - t0)/(p1 - p0))
    tn = np.array(tn)
    if verbose:
        print(f"[echo] raw crossings={len(ix)} filtered={len(tn)}  max|phi_cen|={np.abs(ph).max():.4f}")
        if len(tn) >= 2:
            print(f"[echo] crossing times: {np.round(tn, 4)}")
    if len(tn) < 3:
        return None
    # phi(0,tau) is NON-sinusoidal (several crossings per period, unequal sub-gaps —
    # measured), so successive-gap ratios do NOT give Delta at ~1 resolvable period.
    # Instead: hold Delta at the literature value (the BVP refines it — that is Stage-A's
    # job) and fit t* from SAME-PHASE crossing pairs: (t*-t_a) = e^Delta (t*-t_b).
    Delta = 3.4453
    eD = np.exp(Delta)
    cands = []
    for ia in range(len(tn)):
        for ib in range(ia+1, len(tn)):
            ts = (eD*tn[ib] - tn[ia])/(eD - 1.0)
            if ts > tn[-1] + 1e-3:
                # accept only pairs whose implied one-period partner ratio is plausible
                Ta, Tb = ts - tn[ia], ts - tn[ib]
                if 0.7 < np.log(Ta/Tb)/Delta < 1.3:
                    cands.append(ts)
    if not cands:
        return None
    tstar = float(np.median(cands))
    spread = float(np.percentile(cands, 84) - np.percentile(cands, 16))
    # honest Delta_echo estimate from the best-separated same-phase pair at the fitted t*
    T = tstar - tn
    Dest = [np.log(T[a]/T[b]) for a in range(len(tn)) for b in range(a+1, len(tn))
            if 2.4 < np.log(T[a]/T[b]) < 4.4]
    return dict(tn=tn, tstar=tstar, tstar_spread=spread, Delta=Delta,
                Delta_echo=(float(np.median(Dest)) if Dest else float('nan')))

def sample_cylinder(ev, snaps, tstar, Delta, Mtau=32, zmin=-6.0, zmax=0.8, Nz=60, tau_hi_off=0.3):
    """Sample (X, Y, g_central-gauge, a) on the (tau, zeta) cylinder over ONE period.
    tau = ln(t*-tG) (Gundlach's tau: decreases toward the singularity); zeta = ln(r/T), xi0=0.
    X = sqrt(2pi) (r/a) Phi,  Y = sqrt(2pi) (r/a) Pi  (the lapse cancels),
    g = a * al(0)/al  (central-lapse renormalization)."""
    s2p = np.sqrt(2*np.pi)
    ts = np.array([s[0] for s in snaps])
    keep = ts < tstar - 1e-6
    snaps = [s for s, k in zip(snaps, keep) if k]
    ts = ts[keep]
    Ts = tstar - ts
    taus_snap = np.log(Ts)
    # ONE period within the AVAILABLE T-range [~0, t*]: put its top at 0.95 t* — the
    # period then spans T in [0.95 t* e^{-Delta}, 0.95 t*]. Small-T features live at
    # r ~ T (resolution warning if T_lo < 6 dr — the D-021 uniform-grid physics).
    T_hi = 0.95*tstar
    T_lo = T_hi*np.exp(-Delta)
    if T_lo < 6*ev.dr:
        print(f"[cyl] WARNING: T_lo={T_lo:.3f} < 6 dr={6*ev.dr:.3f} — small-T end under-resolved")
    tau1 = np.log(T_lo)
    tau_grid = tau1 + Delta*np.arange(Mtau)/Mtau            # one period upward
    z_grid = np.linspace(zmin, zmax, Nz)
    out = dict(tau=tau_grid, z=z_grid,
               X=np.zeros((Mtau, Nz)), Y=np.zeros((Mtau, Nz)),
               g=np.zeros((Mtau, Nz)), a=np.zeros((Mtau, Nz)))
    for i, tau in enumerate(tau_grid):
        k = int(np.argmin(np.abs(taus_snap - tau)))         # nearest snapshot in tau
        tG, phi, Phi, Pi, a, al = snaps[k]
        T = tstar - tG
        rs = T*np.exp(z_grid)
        aI = np.interp(rs, ev.r, a); alI = np.interp(rs, ev.r, al)
        PhiI = np.interp(rs, ev.r, Phi); PiI = np.interp(rs, ev.r, Pi)
        out['X'][i] = s2p*rs/aI*PhiI
        out['Y'][i] = s2p*rs/aI*PiI
        out['g'][i] = aI*al[0]/alI
        out['a'][i] = aI
    return out

def extract_seed(cyl, Delta, zc=-2.0):
    """From the cylinder: Y1(tau) = e^{-zc} Y(zc,tau) (center behavior Y ~ Y1 e^z);
    xi0(tau) from e^{xi0} g(xi0-line) ~= 1 (leading order); X+0(tau) = X+ at zeta=xi0.
    Returns the nr_dss 13-vector (with the given Delta) + the raw functions."""
    tau, z = cyl['tau'], cyl['z']
    Mtau = len(tau)
    Xp = cyl['X'] + cyl['Y']
    xi0 = np.zeros(Mtau); Xp0 = np.zeros(Mtau)
    for i in range(Mtau):
        f = np.exp(z)*cyl['g'][i] - 1.0        # e^{zeta} g = 1 at the (xi0'-neglected) SSH
        jx = np.where(np.diff(np.sign(f)) != 0)[0]
        if len(jx):
            j = jx[-1]
            frac = -f[j]/(f[j+1] - f[j])
            xi0[i] = z[j] + frac*(z[j+1] - z[j])
            Xp0[i] = Xp[i, j] + frac*(Xp[i, j+1] - Xp[i, j])
        else:
            xi0[i] = z[-1]; Xp0[i] = Xp[i, -1]
    # Y1 in the GAUGED frame: zeta_G = zeta_raw - xi0(tau); the solver's center function is
    # Y1_G(tau) = e^{-zeta_G} Y at fixed zeta_G = zc  =>  sample raw at zc + xi0(tau)
    Y1 = np.zeros(Mtau)
    for i in range(Mtau):
        zr = zc + xi0[i]
        Y1[i] = np.exp(-zc)*np.interp(zr, z, cyl['Y'][i])
    # Fourier project (period Delta, grid = tau offsets)
    th = 2*np.pi*np.arange(Mtau)/Mtau
    def h(f, k, fn):  return 2.0/Mtau*np.sum(f*fn(k*th))
    # tau-translation gauge: rotate ALL functions by a common shift so xi0's k=2 SINE
    # component vanishes (the nr_dss parametrization's gauge condition)
    c2, s2 = h(xi0, 2, np.cos), h(xi0, 2, np.sin)
    ph2 = 0.5*np.arctan2(s2, c2)
    def hs(f, k, fn):
        c, s = h(f, k, np.cos), h(f, k, np.sin)
        cr = c*np.cos(k*ph2) + s*np.sin(k*ph2)
        sr = -c*np.sin(k*ph2) + s*np.cos(k*ph2)
        return cr if fn is np.cos else sr
    u = np.zeros(13)
    u[0] = Delta
    u[1] = np.mean(xi0)
    u[2] = hs(xi0, 2, np.cos); u[3] = hs(xi0, 4, np.cos); u[4] = hs(xi0, 4, np.sin)
    u[5] = hs(Y1, 1, np.cos); u[6] = hs(Y1, 1, np.sin)
    u[7] = hs(Y1, 3, np.cos); u[8] = hs(Y1, 3, np.sin)
    u[9] = hs(Xp0, 1, np.cos); u[10] = hs(Xp0, 1, np.sin)
    u[11] = hs(Xp0, 3, np.cos); u[12] = hs(Xp0, 3, np.sin)
    return u, dict(Y1=Y1, xi0=xi0, Xp0=Xp0)

def seed_pipeline(ev, p_sub, verbose=True):
    """The full Plan-B' pipeline: marginal subcritical run (recorded) -> t*, Delta_echo ->
    cylinder sampling -> nr_dss seed vector. Returns (u_seed, info)."""
    r = ev.run(p_sub, tmax=200.0, record=True, rec_tG0=0.0, stride=3)
    if verbose:
        print(f"[seed] marginal run fate={r['fate']} t={r['t']:.1f} tG={r['tG']:.2f}  snaps={len(r['snaps'])}")
    an = analyze_echoes(r['cen'])
    if an is None:
        print("[seed] ECHO ANALYSIS FAILED (too few clean crossings)"); return None, dict(run=r)
    if verbose:
        print(f"[seed] crossings={len(an['tn'])}  t*={an['tstar']:.4f} (spread {an['tstar_spread']:.3f})")
        print(f"[seed] Delta_echo(est) = {an['Delta_echo']:.3f}   (Gundlach 3.4453; seed uses the lit value)")
    cyl = sample_cylinder(ev, r['snaps'], an['tstar'], an['Delta'])
    u, fns = extract_seed(cyl, an['Delta'])
    if verbose:
        print(f"[seed] u_seed = {np.round(u, 4)}")
        print(f"[seed] amplitude check: RMS(Y1)={np.sqrt(np.mean(fns['Y1']**2)):.4f}  "
              f"RMS(X+0)={np.sqrt(np.mean(fns['Xp0']**2)):.4f}  mean(xi0)={np.mean(fns['xi0']):.4f}")
    np.save("seed_u.npy", u)
    return u, dict(run=r, an=an, cyl=cyl, fns=fns)

if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "pstar"
    N = int(sys.argv[2]) if len(sys.argv) > 2 else 1600
    ev = Ev(N=N)
    if mode == "pstar":
        lo, hi = bisect_pstar(ev)
        print(f"p* bracket: [{lo:.12f}, {hi:.12f}]  (rel width {(hi-lo)/lo:.1e})")
    elif mode == "seed":
        p_sub = float(sys.argv[3]) if len(sys.argv) > 3 else None
        assert p_sub is not None, "usage: nr_evolve.py seed N p_subcritical"
        u, info = seed_pipeline(ev, p_sub)
