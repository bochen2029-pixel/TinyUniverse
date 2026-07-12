# stageA.py — CROSS THE SONIC POINT. Use the exact 2nd-order local series at x=0 to step off both
# sides, integrate the center side inward (x -> -inf), verify A->1, V->0. This IS the sonic crossing.
#
# Exact sonic point:  V0=1/sqrt3, N0=2/sqrt3, A0=3/2, om0=3/4 ; A'=3, N'=-2/sqrt3.
# Two separatrix branches (from L'Hopital quadratic, verified sonic_exact.py):
#   B1: om'=9/2,  V'=2/sqrt3 ; 2nd order N2=2/sqrt3, A2=33, om2=99/2, V2=10*sqrt3/3
#   B2: om'=27/2, V'=4/sqrt3 ; 2nd order N2=-10*sqrt3/3, A2=87, om2=591/2, V2=12*sqrt3
import css_core as C
import math

S3 = math.sqrt(3.0)
PT = (2.0/S3, 1.5, 0.75, 1.0/S3)   # (N0,A0,om0,V0)
NP1, AP1 = -2.0/S3, 3.0            # N', A' (both branches share these)

BRANCHES = {
    'B1': dict(N1=-2.0/S3, A1=3.0, om1=4.5,  V1=2.0/S3,
               N2=2.0/S3,  A2=33.0, om2=49.5, V2=10.0*S3/3.0),
    'B2': dict(N1=-2.0/S3, A1=3.0, om1=13.5, V1=4.0/S3,
               N2=-10.0*S3/3.0, A2=87.0, om2=295.5, V2=12.0*S3),
}

def series_state(br, x):
    """2nd-order local series Y(x) = Y0 + Y1 x + Y2 x^2/2 about the sonic point."""
    b = BRANCHES[br]
    N = PT[0] + b['N1']*x + 0.5*b['N2']*x*x
    A = PT[1] + b['A1']*x + 0.5*b['A2']*x*x
    om = PT[2] + b['om1']*x + 0.5*b['om2']*x*x
    V = PT[3] + b['V1']*x + 0.5*b['V2']*x*x
    return (N, A, om, V)

def integrate_from(br, side, eps=1e-3, xstop=-16.0, h=2e-4):
    """Seed from series at x=side*eps, integrate toward xstop (center if xstop<0)."""
    x = side*eps
    Y = series_state(br, x)
    d = -1.0 if xstop < x else 1.0
    n = int(abs((xstop - x)/h))
    traj = [(x, *Y)]
    for i in range(n):
        try:
            Y = C.rk4_step(Y, d*h)
        except (ValueError, ZeroDivisionError, OverflowError):
            return traj, 'blow'
        x += d*h
        if not all(math.isfinite(v) for v in Y):
            return traj, 'blow'
        traj.append((x, *Y))
        if Y[1] <= 0 or Y[1] > 1e4:
            return traj, 'edge'
        # near center: V->0, A->1, N large
        if d < 0 and Y[0] > 5e3:
            return traj, 'center'
    return traj, 'done'

print("STAGE A — sonic crossing via exact 2nd-order local series. Step off x=0, integrate inward.")
print("Target center: A->1, V->0, N->inf.\n")
for br in ['B1', 'B2']:
    for side in [-1.0, +1.0]:
        traj, status = integrate_from(br, side)
        e = traj[-1]
        tag = ''
        if side < 0 and abs(e[4]) < 5e-3 and abs(e[2]-1.0) < 0.05:
            tag = '  <=== REGULAR CENTER (A->1, V->0)'
        print(f"  {br} side={side:+.0f}: end x={e[0]:+.3f} N={e[1]:.4f} A={e[2]:.5f} om={e[3]:.6f} V={e[4]:+.6f}  [{status}]{tag}")
    print()
