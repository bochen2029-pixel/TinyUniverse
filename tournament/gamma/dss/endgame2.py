# endgame2.py — session-6 LOW-K-pinned release ladder (the v2 endgame).
#
# Why v1 failed (endgame.log, all measured): the total-RMS pin constrains WHICH amplitude
# but not WHICH HARMONICS carry it — the (64,20,21) grind emptied k<7 while holding the
# pin (impostor #4 recurring at higher truncation: g collapsed [0.52,1.06]->[0.978,1.0],
# tail INVERTED to edge-piling). Pin-release then converged to the VACUUM (|r|=1e-13
# machine-zero, all-zero tails = the house tell), and the released Delta sat at its FROZEN
# INITIALIZATION (zero gradient in the vacuum) — "PRELIMINARY Delta = 3.445300" was the
# input echoed back to 7 digits, not a measurement.
#
# v2 discipline:
#   - LOW-K pin (nr_relax pin=('lowk',c,w)): RMS of the k<=5 odd SSH X+ COEFFICIENTS —
#     immune to high-k relocation.
#   - battery-GATED stages: abort the moment gdev collapses / the tail inverts /
#     machine-zero appears, instead of discovering the drain at the end.
#   - pin released by a WEIGHT RAMP (w: 30 -> 10 -> 3 -> 1 -> 0.3 -> 0.1 -> 0), battery at
#     each rung. Drain at every rung = the state is pin-supported = off-manifold verdict
#     (honest wall), not a Delta.
#   - Delta release only from a battery-green unpinned state; a released Delta that does
#     not MOVE from 3.4453 at all is flagged ZERO-INFORMATION (suspiciously-exact law).
#
# usage: python endgame2.py smoke   (wiring check, no solves, ~30 s)
#        python endgame2.py run     (the full ladder; logs to stdout)
import numpy as np, sys, time
import nr_relax as R
from newton_drive import lm

np.seterr(over='raise', invalid='raise', divide='raise', under='ignore')

LIT = 3.4453          # Gundlach 3.4453 +- 0.0005; the FROZEN Delta during amplitude work
FLOOR = 5e-3          # |r| a solution must reach before any Delta claim (run48 hole: the
                      # battery said GREEN on the 7.9e-2 PLATEAU — a non-solution — and the
                      # released Delta "measured" its own initialization to +9e-8. A Delta
                      # is only information when |r| is near the truncation floor.)

def pad_even(v, KE1, KE2):
    out = np.zeros(1 + KE2); out[:1 + KE1] = v; return out

def pad_odd(v, KO1, KO2):
    out = np.zeros(KO2 + 1); out[:KO1 + 1] = v; return out

def pack(rx, Delta, xi0c, gc, xpc, xmc):
    """Inverse of Relax.unpack under the CURRENT configuration (xi0 s2 slot dropped)."""
    u = np.zeros(rx.nu)
    u[0] = Delta
    u[1] = xi0c[0]; u[2] = xi0c[1]; u[3:R.NE] = xi0c[3:]
    u[R.NE:] = np.concatenate([gc, xpc, xmc], axis=1).ravel()
    return u

def lowk_amp(rx, u):
    _, _, _, xpc, _ = rx.unpack(u)
    return float(np.sqrt(np.mean(xpc[-1][:6]**2)))

def battery(rx, u, label, rn=None):
    """The which-solution battery, printed + gated. Healthy reference (v3_48):
    g=[0.521,1.064], tail [0.0499 ... 0.015] decaying, lowk-dominant."""
    D, xi0c, gc, xpc, xmc = rx.unpack(u)
    Be, Bde, Bo, Bdo, Pe, Po = R.bases(D)
    g = np.exp(gc@Be.T)
    mags = np.abs(xpc).max(axis=0)
    tail = [float(max(mags[2*j], mags[2*j + 1])) for j in range(R.NO//2)]
    lowk = max(tail[:3]); highk = max(tail[3:]) if len(tail) > 3 else 0.0
    checks = {
        'strong_field(g.min<0.90)': g.min() < 0.90,
        'lowk_dominant':            lowk > highk,
        'tail_decays':              (lowk > 0) and (tail[-1] < 0.7*lowk),
        'not_vacuum(maxtail>1e-8)': max(tail) > 1e-8,
        'Delta_band[3.40,3.49]':    3.40 <= D <= 3.49,
        'not_Delta_half':           abs(D - LIT/2) > 0.1,
    }
    if rn is not None and label.endswith('FINAL'):
        checks[f'at_floor(|r|<{FLOOR:g})'] = rn < FLOOR
    ok = all(checks.values())
    ts = " ".join(f"{t:.4f}" for t in tail)
    extra = f"  |r|={rn:.3e}" if rn is not None else ""
    print(f"  [battery:{label}] Delta={D:.7f}  g=[{g.min():.3f},{g.max():.3f}]  tail=[{ts}]{extra}", flush=True)
    print(f"  [battery:{label}] {'GREEN' if ok else 'FAIL: ' + ', '.join(k for k, v in checks.items() if not v)}", flush=True)
    return ok

def load48():
    R.configure(48, 14, 13)
    rx1 = R.Relax(Nz=40)
    u48 = np.load('v3_48.npy')
    return rx1, u48

def pad_to_64(rx1, u):
    D_, xi0c_, gc_, xpc_, xmc_ = rx1.unpack(u)
    R.configure(64, 20, 21)
    rx2 = R.Relax(Nz=40)
    u2 = pack(rx2, D_,
              pad_even(xi0c_, 14, 20),
              np.stack([pad_even(r_, 14, 20) for r_ in gc_]),
              np.stack([pad_odd(r_, 13, 21) for r_ in xpc_]),
              np.stack([pad_odd(r_, 13, 21) for r_ in xmc_]))
    return rx2, u2

def smoke():
    R.configure(48, 14, 13); R.vacuum_control()
    R.configure(64, 20, 21); R.vacuum_control()
    rx1, u48 = load48()
    c48 = lowk_amp(rx1, u48)
    r0 = rx1.residual(u48)
    rp = rx1.residual(u48, pin=('lowk', c48, 30.0))
    print(f"[smoke] v3_48 |r|={np.linalg.norm(r0):.4e}  lowk_amp={c48:.6f}  "
          f"pin_row@c={rp[-1]:.2e} (must be 0)", flush=True)
    battery(rx1, u48, 'v3_48-asis')
    rx2, u2 = pad_to_64(rx1, u48)
    print(f"[smoke] padded (64,20,21) |r|={np.linalg.norm(rx2.residual(u2)):.4f}  "
          f"(v1 endgame measured 1.1385 on the same state — must match)", flush=True)
    print(f"[smoke] lowk_amp after pad = {lowk_amp(rx2, u2):.6f} (must equal 48-value)", flush=True)

def run():
    t00 = time.time()
    R.configure(48, 14, 13); R.vacuum_control()
    R.configure(64, 20, 21); R.vacuum_control()

    # ---- stage 1: (48,14,13) lowk-pin sanity regrind from the healthy state
    rx1, u48 = load48()
    c48 = lowk_amp(rx1, u48)
    print(f"[stage1-48] lowk pin c={c48:.6f} (measured from v3_48)", flush=True)
    u, rn, st, _ = lm(rx1, u48, pin=('lowk', c48, 30.0), freeze_delta=LIT, max_iter=20, verbose=False)
    print(f"[stage1-48] {st}  |r|={rn:.4e}", flush=True)
    np.save('lowk48.npy', u)
    if not battery(rx1, u, 'post-48-lowk', rn):
        print("[ABORT] the 48-truncation lowk regrind failed the battery", flush=True); return

    # ---- stage 2: pad -> (64,20,21) grind, lowk pin ON, Delta frozen, gated per round
    rx2, u2 = pad_to_64(rx1, u)
    print(f"[stage2] padded |r| = {np.linalg.norm(rx2.residual(u2)):.4f}", flush=True)
    for rnd in range(6):
        u2, rn, st, _ = lm(rx2, u2, pin=('lowk', c48, 30.0), freeze_delta=LIT, max_iter=25, verbose=False)
        print(f"[stage2] round {rnd}: {st}  |r|={rn:.4e}", flush=True)
        np.save('lowk64.npy', u2)
        if not battery(rx2, u2, f'64-grind-r{rnd}', rn):
            print("[ABORT] the 64-grind left the healthy manifold (drain caught mid-flight)", flush=True); return
        if st == 'converged' or (st == 'slow' and rnd >= 1):
            break

    # ---- stages 3+4: ramp + release
    u2, verdict = ramp_release(rx2, u2, c48, 'lowk64')
    print(f"[done] verdict={verdict}  {time.time() - t00:.0f}s", flush=True)

def ramp_release(rx, u2, c, tag):
    """Pin-weight RAMP (Delta frozen) then Delta release, battery-gated per rung.
    Returns (u, 'drained' | 'released-fail' | 'released-green')."""
    for w in (10.0, 3.0, 1.0, 0.3, 0.1, 0.0):
        pin = ('lowk', c, w) if w > 0 else None
        u2, rn, st, _ = lm(rx, u2, pin=pin, freeze_delta=LIT, max_iter=25, verbose=False, lam0=1e-2)
        print(f"[{tag}-ramp] w={w:g}: {st}  |r|={rn:.4e}  lowk_amp={lowk_amp(rx, u2):.6f}", flush=True)
        np.save(f'{tag}_ramp.npy', u2)
        if not battery(rx, u2, f'{tag}-w{w:g}', rn):
            print(f"[VERDICT] {tag}: drained at pin weight w={w:g} — the state is PIN-SUPPORTED "
                  f"(off-manifold at this truncation). None faked.", flush=True)
            return u2, 'drained'
    np.save(f'{tag}_nopin.npy', u2)
    # Delta release (verbose — watch Delta actually MOVE)
    u2, rn, st, _ = lm(rx, u2, pin=None, freeze_delta=None, max_iter=40, verbose=True, lam0=1e-2)
    np.save(f'{tag}_final.npy', u2)
    ok = battery(rx, u2, f'{tag}-FINAL', rn)
    D = u2[0]
    print(f"[{tag}] Delta-released: {st}  |r|={rn:.4e}", flush=True)
    if ok:
        print(f"*** PRELIMINARY Delta[{tag}] = {D:.7f}  (drift from frozen {LIT}: {D - LIT:+.3e}; "
              f"lit 3.4453(5); impostor Delta/2 = {LIT/2:.4f}) ***", flush=True)
        if abs(D - LIT) < 1e-9:
            print("!!! Delta did NOT move from its initialization — ZERO-INFORMATION, "
                  "inspect the Delta gradient before believing anything", flush=True)
        return u2, 'released-green'
    print(f"[VERDICT] {tag}: survived the ramp but the FINAL state fails the battery — honest wall, none faked", flush=True)
    return u2, 'released-fail'

def run48():
    """Release ladder AT the native (48,14,13) — the cheapest decisive measurement.
    v1/v2-64 both padded BEFORE releasing; nobody has ever tried releasing the healthy
    state at its own truncation. If it is on-manifold there, the pin comes off and Delta
    moves to a preliminary (48,14,13) value; truncation refinement (M-first-then-K)
    comes after. If it drains even here, the state is pin-supported at its native
    truncation and seed/extraction is still the wall — either way, a verdict."""
    t00 = time.time()
    R.configure(48, 14, 13); R.vacuum_control()
    rx1, u48 = load48()
    c48 = lowk_amp(rx1, u48)
    print(f"[run48] lowk pin c={c48:.6f}", flush=True)
    u, rn, st, _ = lm(rx1, u48, pin=('lowk', c48, 30.0), freeze_delta=LIT, max_iter=20, verbose=False)
    print(f"[run48] sanity regrind: {st}  |r|={rn:.4e}", flush=True)
    if not battery(rx1, u, 'run48-sane', rn):
        print("[ABORT] the sanity regrind failed the battery", flush=True); return
    u, verdict = ramp_release(rx1, u, c48, 'r48')
    print(f"[run48] verdict: {verdict}  ({time.time() - t00:.0f}s)", flush=True)

def nzprobe():
    """Plateau vs Nz — the decisive discrimination for the 7.9e-2 plateau at (48,14,13):
    if the plateau DROPS as Nz rises, the z-discretization is inconsistent at Nz=40 (the
    discrete manifold sits ~0.08 from the physical one) and the path is an Nz climb; if it
    stays flat, the plateau is optimizer/valley flatness and the path is Psi-tc."""
    R.configure(48, 14, 13); R.vacuum_control()
    rx40 = R.Relax(Nz=40)
    u40 = np.load('r48_final.npy')
    c = lowk_amp(rx40, u40)
    print(f"[nzprobe] seed = r48_final (post-release, |r|40={np.linalg.norm(rx40.residual(u40)):.4e})", flush=True)
    for Nz in (40, 60, 80):
        rx = R.Relax(Nz=Nz)
        u = u40.copy() if Nz == 40 else R.upsample_u(rx40, u40, rx)
        r0 = np.linalg.norm(rx.residual(u))
        t0 = time.time()
        for rnd in range(2):
            u, rn, st, _ = lm(rx, u, pin=('lowk', c, 3.0), freeze_delta=LIT, max_iter=25, verbose=False)
        print(f"[nzprobe] Nz={Nz}: raw|r|={r0:.4e} -> plateau |r|={rn:.4e}  ({st}, {time.time()-t0:.0f}s)", flush=True)
        battery(rx, u, f'nz{Nz}', rn)
        np.save(f'nz{Nz}.npy', u)

if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else 'smoke'
    {'smoke': smoke, 'run': run, 'run48': run48, 'nzprobe': nzprobe}[mode]()
