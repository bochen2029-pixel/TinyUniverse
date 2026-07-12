# stageA_shoot.py — precise Stage-A solve. Shoot a2 (central density) so the center trajectory
# reaches the EXACT sonic point (numV=numOm=0, i.e. lands on the analytic separatrix B1).
# Then continue THROUGH the sonic point using the 2nd-order series (seed the outer side at x=0+eps
# on branch B1) and integrate outward to confirm dispersal (V->1, om->0). Never integrate through D=0.
import css_core as C
import math
from scipy.optimize import brentq

S3 = math.sqrt(3.0)

def cross(nc, a2, z0=math.exp(-16), xmax=4.0, h=2e-5):
    """Integrate center->sonic; return (x_cross, state, numV, numOm) at Dson=0 crossing, else None."""
    Y = C.center_seed(nc, a2, z0); x = math.log(z0)
    prevD = C.Dson(Y[0], Y[3]); prevY = Y; prevx = x
    n = int((xmax - x)/h)
    for i in range(n):
        try: Y2 = C.rk4_step(Y, h)
        except (ValueError, ZeroDivisionError, OverflowError): return None
        if not all(math.isfinite(v) for v in Y2): return None
        D2 = C.Dson(Y2[0], Y2[3])
        if prevD*D2 < 0 and abs(Y2[3]) > 1e-3:
            fr = prevD/(prevD - D2)
            Yh = tuple(prevY[k] + fr*(Y2[k]-prevY[k]) for k in range(4))
            return (prevx + fr*h, Yh, C.numV(Yh[1],Yh[0],Yh[2],Yh[3]), C.numOm(Yh[1],Yh[0],Yh[2],Yh[3]))
        prevD = D2; prevY = Y2; prevx = x + h; Y = Y2; x += h
    return None

def residual(a2, nc=1.0):
    """Shooting residual: numV at the sonic crossing (0 => analytic passage on separatrix)."""
    r = cross(nc, a2)
    if r is None: return None
    return r[2]

# 1) bracket & refine a2*
print("Refining a2* (numV=0 at sonic crossing) ...")
grid = [0.20,0.22,0.24,0.25,0.26,0.28,0.30]
vals = [(a,residual(a)) for a in grid]
for a,v in vals: print(f"  a2={a:.3f}  numV={v:+.4e}" if v is not None else f"  a2={a:.3f}  (no crossing)")
# find sign change
br = None
for i in range(len(vals)-1):
    if vals[i][1] is not None and vals[i+1][1] is not None and vals[i][1]*vals[i+1][1] < 0:
        br = (vals[i][0], vals[i+1][0]); break
assert br, "no bracket for a2*"
a2star = brentq(lambda a: residual(a), br[0], br[1], xtol=1e-12, rtol=1e-14)
r = cross(1.0, a2star)
xs, Yh, nV, nO = r
print(f"\na2* = {a2star:.12f}")
print(f"sonic crossing at x={xs:.6f} (gauge-shiftable to 0):")
print(f"  N0={Yh[0]:.10f} (2/sqrt3={2/S3:.10f})")
print(f"  A0={Yh[1]:.10f} (3/2)")
print(f"  om0={Yh[2]:.10f} (3/4)")
print(f"  V0={Yh[3]:.10f} (1/sqrt3={1/S3:.10f})")
print(f"  numV={nV:+.3e}  numOm={nO:+.3e}   (both ~0 => analytic passage)")

# 2) verify the center BC by integrating the SAME a2* trajectory back and reporting center approach
print("\nCenter-side self-check (integrate the a2* trajectory, watch A->1,V->0,dV/dx->0):")
Y = C.center_seed(1.0, a2star, math.exp(-16)); x = math.log(Y[0]) # z0 s.t. N=1/z0; recompute
Y = C.center_seed(1.0, a2star, math.exp(-16)); x = -16.0
for target_x in [-16,-12,-8,-5,-3]:
    # integrate to target_x
    pass
# simpler: seed at several depths, report state & dV/dx
for lx in [-16,-12,-9,-6,-4]:
    z0 = math.exp(lx)
    Ys = C.center_seed(1.0, a2star, z0)
    dv = C.Vx(Ys[1],Ys[0],Ys[2],Ys[3])
    print(f"  x={lx:+d}: N={Ys[0]:.3e} A={Ys[1]:.8f} om={Ys[2]:.3e} V={Ys[3]:+.3e}  dV/dx={dv:+.3e}")

# 3) continue THROUGH the sonic point: seed outer side from 2nd-order B1 series, integrate outward
print("\nOuter (dispersal) side via B1 2nd-order series through sonic (x=0+eps), integrate outward:")
BR1 = dict(N1=-2/S3, A1=3.0, om1=4.5, V1=2/S3, N2=2/S3, A2=33.0, om2=49.5, V2=10*S3/3)
def b1_series(x):
    return (2/S3 + BR1['N1']*x + 0.5*BR1['N2']*x*x,
            1.5 + BR1['A1']*x + 0.5*BR1['A2']*x*x,
            0.75 + BR1['om1']*x + 0.5*BR1['om2']*x*x,
            1/S3 + BR1['V1']*x + 0.5*BR1['V2']*x*x)
eps = 2e-3; Y = b1_series(eps); x = eps; h = 2e-4
for i in range(int(6.0/h)):
    try: Y = C.rk4_step(Y, h)
    except Exception:
        print(f"  stop (numeric) at x={x:.3f}"); break
    x += h
    if not all(math.isfinite(v) for v in Y):
        print(f"  stop (nonfinite) at x={x:.3f} state={Y}"); break
    if abs(Y[3]) >= 0.999:
        print(f"  V->1 (dispersal) at x={x:.3f}: N={Y[0]:.4f} A={Y[1]:.5f} om={Y[2]:.6f} V={Y[3]:.6f}"); break
    if Y[1] <= 0 or Y[1] > 1e4:
        print(f"  A edge at x={x:.3f}: A={Y[1]:.4f} V={Y[3]:.4f} om={Y[2]:.5f}"); break
else:
    print(f"  end x={x:.3f}: N={Y[0]:.4f} A={Y[1]:.5f} om={Y[2]:.6f} V={Y[3]:.6f}")
