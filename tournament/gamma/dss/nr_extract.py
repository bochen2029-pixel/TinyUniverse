# nr_extract.py — EXTRACTION REBUILD (session-2 recipe #1): from near-critical evolver
# snapshots to a high-quality relaxation seed. Fixes the measured limiters of the v1
# extraction (nearest-snapshot aliasing; SSH located WITHOUT the xi0' term — but measured
# xi0' swings +-0.8..1.1, mislocating the line by up to ~0.8 in zeta and sampling X+ deep
# inside where it is small = the 3-5x amplitude deficit all six solver campaigns hit).
#
#   - 2-D interpolation: linear in tG between bracketing snapshots, after 1-D r-interp.
#   - xi0(tau) from the FULL coordinate condition (1+xi0') e^{xi0} g(zeta_raw=xi0) = 1,
#     iterated with spectral xi0' (fixed point, ~4 sweeps).
#   - boundary functions by BAND least squares: Y over zeta_G in [-4.5,-2] fit to
#     Y1 e^{z} + Y3 e^{3z}; X+ over [-0.6, 0] fit to X+0 + X+1 z.
#   - the full-field relaxation seed sampled in the GAUGED frame per tau-phase.
import numpy as np

S2P = np.sqrt(2*np.pi)

class Cyl2D:
    """2-D interpolant over the recorded snapshots (fields of tG and r)."""
    def __init__(self, ev, snaps, tstar):
        tG = np.array([s[0] for s in snaps])
        keep = tG < tstar - 1e-9
        self.tG = tG[keep]
        self.PHI = np.array([s[2] for s, k in zip(snaps, keep) if k])
        self.PI  = np.array([s[3] for s, k in zip(snaps, keep) if k])
        self.A   = np.array([s[4] for s, k in zip(snaps, keep) if k])
        self.AL  = np.array([s[5] for s, k in zip(snaps, keep) if k])
        self.ALC = self.AL[:, 0]
        self.r = ev.r
        self.tstar = tstar

    def raw(self, tau, rs):
        """(X, Y, g) at Gundlach-tau and radii rs (1-D array), 2-D interpolated."""
        T = np.exp(tau)
        tg = self.tstar - T
        k = np.searchsorted(self.tG, tg) - 1
        k = min(max(k, 0), len(self.tG) - 2)
        w = (tg - self.tG[k])/(self.tG[k+1] - self.tG[k])
        w = min(max(w, 0.0), 1.0)
        out = []
        for kk, ww in ((k, 1.0 - w), (k + 1, w)):
            a  = np.interp(rs, self.r, self.A[kk])
            al = np.interp(rs, self.r, self.AL[kk])
            Ph = np.interp(rs, self.r, self.PHI[kk])
            Pi = np.interp(rs, self.r, self.PI[kk])
            X = S2P*rs/a*Ph
            Y = S2P*rs/a*Pi
            g = a*self.ALC[kk]/al
            out.append((ww, X, Y, g))
        X = out[0][0]*out[0][1] + out[1][0]*out[1][1]
        Y = out[0][0]*out[0][2] + out[1][0]*out[1][2]
        g = out[0][0]*out[0][3] + out[1][0]*out[1][3]
        return X, Y, g

def dtau_spectral(f, Delta):
    F = np.fft.rfft(f)
    k = np.arange(len(F))
    return np.fft.irfft(1j*(2*np.pi*k/Delta)*F, len(f))

def extract(cyl, Delta, Mph=64, T_lo=0.55, zg=np.linspace(-6.0, 0.9, 140)):
    """Full extraction on Mph tau-phases over one period starting at T_lo."""
    tau0 = np.log(T_lo)
    taus = tau0 + Delta*np.arange(Mph)/Mph
    # raw fields on a (tau, zeta_raw) grid
    X = np.zeros((Mph, len(zg))); Y = np.zeros_like(X); G = np.zeros_like(X)
    for i, tau in enumerate(taus):
        rs = np.exp(tau)*np.exp(zg)
        X[i], Y[i], G[i] = cyl.raw(tau, rs)
    # xi0 from the FULL coordinate condition, iterated with spectral xi0'. Stabilized
    # (the naive iteration diverged, measured: xi0 -> +29, dxi0 -> +-1200): interpolation-
    # only crossings in a physical window, LOW-HARMONIC projection before differentiating,
    # damped fixed point.
    win = (zg >= -0.5) & (zg <= 0.9)                  # physical SSH window in raw zeta
    zw = zg[win]
    def locate(target_row, Grow, prev):
        f = np.exp(zw)*Grow[win] - target_row
        jx = np.where(np.diff(np.sign(f)) != 0)[0]
        if not len(jx):
            return prev
        j = jx[-1]
        fr = float(np.clip(-f[j]/(f[j+1] - f[j]), 0.0, 1.0))
        return zw[j] + fr*(zw[j+1] - zw[j])
    def smooth_even(v, kmax=6):
        F = np.fft.rfft(v); kk = np.arange(len(F))
        return np.fft.irfft(F*((kk % 2 == 0) & (kk <= kmax)), Mph)
    xi0 = np.zeros(Mph)
    for i in range(Mph):
        xi0[i] = locate(1.0, G[i], 0.5)
    xi0 = smooth_even(xi0)
    for _ in range(5):
        dx = dtau_spectral(xi0, Delta)
        target = 1.0/np.clip(1.0 + dx, 0.2, 5.0)
        new = np.array([locate(target[i], G[i], xi0[i]) for i in range(Mph)])
        xi0 = smooth_even(0.5*xi0 + 0.5*new)
    dxi0 = dtau_spectral(xi0, Delta)
    # gauged-frame field rows: per tau, resample at zeta_raw = zeta_G + xi0(tau)
    def gauged(zG_arr):
        Xg = np.zeros((Mph, len(zG_arr))); Yg = np.zeros_like(Xg); Gg = np.zeros_like(Xg)
        for i, tau in enumerate(taus):
            zr = zG_arr + xi0[i]
            rs = np.exp(tau)*np.exp(zr)
            Xg[i], Yg[i], Gg[i] = cyl.raw(tau, rs)
        return Xg, Yg, Gg
    # boundary functions by band least squares — RESOLUTION-MASKED per sample: the deep
    # center at small T samples r = T e^{z+xi0} BELOW the grid (measured: r ~ 0.008 vs
    # dr = 0.0375 inflated Y1 3x); only r >= rmin points enter the fits.
    rmin = 4.0*cyl.r[0]*2.0                           # ~4 dr (staggered r[0] = dr/2)
    zb = np.linspace(-4.5, -2.0, 24)
    Xb, Yb, Gb = gauged(zb)
    # per-phase Y1 estimate + resolvability weight; ~13% of phases (small T) have NO
    # resolvable center band (their one-period partners lie beyond t*) — their Y1 is
    # RECONSTRUCTED by the periodicity: weighted harmonic LSQ over the good phases.
    y1_est = np.zeros(Mph); wgt = np.zeros(Mph)
    for i, tau in enumerate(taus):
        rs = np.exp(tau)*np.exp(zb + xi0[i])
        m = rs >= rmin
        if m.sum() >= 5:
            A = np.column_stack([np.exp(zb[m]), np.exp(3*zb[m])])
            y1_est[i] = np.linalg.lstsq(A, Yb[i][m], rcond=None)[0][0]
            wgt[i] = m.sum()
    th = 2*np.pi*np.arange(Mph)/Mph
    H = np.column_stack([np.cos(th), np.sin(th), np.cos(3*th), np.sin(3*th),
                         np.cos(5*th), np.sin(5*th)])
    good = wgt > 0
    W = np.sqrt(wgt[good])
    coef = np.linalg.lstsq(H[good]*W[:, None], y1_est[good]*W, rcond=None)[0]
    Y1 = H@coef
    zs = np.linspace(-0.6, 0.0, 14)
    Xs, Ys, Gs = gauged(zs)
    As = np.column_stack([np.ones_like(zs), zs])
    coefXp = np.linalg.lstsq(As, (Xs + Ys).T, rcond=None)[0]
    Xp0 = coefXp[0]
    return dict(taus=taus, xi0=xi0, dxi0=dxi0, Y1=Y1, Xp0=Xp0, gauged=gauged, rmin=rmin,
                Delta=Delta,
                rmsY1=float(np.sqrt(np.mean(Y1**2))), rmsXp0=float(np.sqrt(np.mean(Xp0**2))))

def build_seed(rx, R, ext, Delta):
    """Full-field relaxation seed: EVERY node row fit by resolvability-WEIGHTED harmonic
    least squares over the tau-phases (holes at unresolvable phases are filled by the
    periodic parity basis itself — no per-phase blend seams). Deep center (z < -3.6,
    asymptotics exact) synthesized from the reconstructed Y1 harmonics."""
    Mph = len(ext['xi0'])
    taus = ext['taus']; xi0v = ext['xi0']; dxi0v = ext['dxi0']
    Y1v = ext['Y1']; dY1v = dtau_spectral(Y1v, Delta)
    rmin = ext['rmin']
    # gauge rotation angle (kills xi0's k=2 sine)
    th = 2*np.pi*np.arange(Mph)/Mph
    c2 = 2*np.mean(xi0v*np.cos(2*th)); s2 = 2*np.mean(xi0v*np.sin(2*th))
    ph2 = 0.5*np.arctan2(s2, c2)
    def rot(f):
        F = np.fft.rfft(f); kk = np.arange(len(F))
        return np.fft.irfft(F*np.exp(1j*kk*ph2), len(f))
    # parity bases evaluated ON THE PHASE GRID (theta already includes the rotation via
    # rotating the VALUES; bases in the standard frame)
    phi = th*Delta/(2*np.pi)                          # tau-offsets over one period
    w0 = 2*np.pi/Delta
    Bee = np.zeros((Mph, R.NE)); Bee[:, 0] = 1.0
    j = 1
    for k in range(2, R.KE+1, 2):
        Bee[:, j] = np.cos(k*w0*phi); Bee[:, j+1] = np.sin(k*w0*phi); j += 2
    Boo = np.zeros((Mph, R.NO))
    j = 0
    for k in range(1, R.KO+1, 2):
        Boo[:, j] = np.cos(k*w0*phi); Boo[:, j+1] = np.sin(k*w0*phi); j += 2
    def wfit(vals, wts, B):
        W = np.sqrt(np.maximum(wts, 1e-12))[:, None]
        return np.linalg.lstsq(B*W, vals*W[:, 0], rcond=None)[0]
    # gauged data rows + weights per node
    Xg, Yg, Gg = ext['gauged'](rx.z)
    u = np.zeros(rx.nu)
    u[0] = Delta
    Be, Bde, Bo, Bdo, Pe, Po = R.bases(Delta)
    xc = wfit(rot(xi0v), np.ones(Mph), Bee)
    u[1] = xc[0]; u[2] = xc[1]; u[3:R.NE] = xc[3:]
    F = np.zeros((rx.Nz, R.NPF))
    for k, zG in enumerate(rx.z):
        rs = np.exp(taus)*np.exp(zG + xi0v)
        wts = np.clip(np.log(np.maximum(rs, 1e-12)/rmin)/np.log(3.0), 0.0, 1.0)
        # per-PHASE blend: data where resolvable, asymptotic synthesis filling the holes —
        # the mixture fraction shifts smoothly with depth (a per-NODE binary switch put an
        # 18.7 g-eq spike at the z=-3.6 seam; per-phase blending removes the seam)
        Ysyn = Y1v*np.exp(zG)
        Xsyn = np.exp(zG)*(np.exp(xi0v)/3.0)*(dY1v - (1.0 + dxi0v)*Y1v)
        Gsyn = 1.0 - (Ysyn*Ysyn)/3.0
        gv = wts*Gg[:, k] + (1 - wts)*Gsyn
        xv = wts*Xg[:, k] + (1 - wts)*Xsyn
        yv = wts*Yg[:, k] + (1 - wts)*Ysyn
        gv, xv, yv = rot(gv), rot(xv), rot(yv)
        F[k, :R.NE] = wfit(gv, np.ones(Mph), Bee)
        F[k, R.NE:R.NE+R.NO] = wfit(xv + yv, np.ones(Mph), Boo)
        F[k, R.NE+R.NO:] = wfit(xv - yv, np.ones(Mph), Boo)
    u[R.NE:] = F.ravel()
    return u
