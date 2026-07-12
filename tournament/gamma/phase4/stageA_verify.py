# stageA_verify.py — Stage-A self-checks: a2* convergence (G-CONVERGE), physical anchors.
# mass aspect m/r: for metric ds^2=-alpha^2 dt^2 + a^2 dr^2 + r^2 dOmega^2, the Misner-Sharp mass
#   2m/r = 1 - 1/a^2 = 1 - 1/A.   At sonic (A=3/2): 2m/r = 1-2/3 = 1/3 -> m/r = 1/6.
# The similarity exponent n: EC use xi=r/t^n; KHA use x=ln(-r/t). The homothetic structure gives a
# relation; the robust invariant we can quote is the sonic-point m/r and the exact algebraic state.
import css_core as C
import math
from scipy.optimize import brentq
S3 = math.sqrt(3.0)

def cross(nc, a2, z0, h):
    Y = C.center_seed(nc, a2, z0); x = math.log(z0)
    prevD = C.Dson(Y[0], Y[3]); prevY = Y; prevx = x
    n = int((4.0 - x)/h)
    for i in range(n):
        try: Y2 = C.rk4_step(Y, h)
        except Exception: return None
        if not all(math.isfinite(v) for v in Y2): return None
        D2 = C.Dson(Y2[0], Y2[3])
        if prevD*D2 < 0 and abs(Y2[3]) > 1e-3:
            fr = prevD/(prevD - D2)
            Yh = tuple(prevY[k] + fr*(Y2[k]-prevY[k]) for k in range(4))
            return (prevx+fr*h, Yh, C.numV(Yh[1],Yh[0],Yh[2],Yh[3]))
        prevD = D2; prevY = Y2; prevx = x+h; Y = Y2; x += h
    return None

def a2star(z0, h):
    f = lambda a: cross(1.0, a, z0, h)[2]
    return brentq(f, 0.23, 0.27, xtol=1e-13, rtol=1e-15)

print("G-CONVERGE: a2* under refinement of (step h, seed depth z0):")
print(f"{'h':>10} {'z0=e^':>8} {'a2*':>16}   {'V0':>12} {'A0':>10}")
base = None
rows = []
for lz in [-14, -16, -18]:
    for h in [4e-5, 2e-5, 1e-5]:
        z0 = math.exp(lz)
        a = a2star(z0, h)
        r = cross(1.0, a, z0, h)
        rows.append((h, lz, a, r[1][3], r[1][1]))
        print(f"{h:10.1e} {lz:8d} {a:16.12f}   {r[1][3]:12.8f} {r[1][1]:10.6f}")
avals = [r[2] for r in rows]
spread = max(avals) - min(avals)
print(f"\na2* spread across refinement ladder = {spread:.3e}")

# best (finest) solution
h = 1e-5; z0 = math.exp(-18)
a = a2star(z0, h); r = cross(1.0, a, z0, h)
N0,A0,om0,V0 = r[1]
print(f"\nFinest Stage-A critical solution (a2*={a:.12f}):")
print(f"  sonic state: N0={N0:.10f}, A0={A0:.10f}, om0={om0:.10f}, V0={V0:.10f}")
print(f"  exact:       N0={2/S3:.10f}, A0=1.5, om0=0.75, V0={1/S3:.10f}")
print(f"  |dev| from exact: dN={abs(N0-2/S3):.2e} dA={abs(A0-1.5):.2e} dom={abs(om0-0.75):.2e} dV={abs(V0-1/S3):.2e}")
print(f"\nPhysical anchors:")
print(f"  Misner-Sharp 2m/r at sonic = 1 - 1/A0 = {1-1/A0:.8f}  (exact 1/3={1/3:.8f}) -> m/r = {(1-1/A0)/2:.8f} (exact 1/6={1/6:.8f})")
print(f"  sound speed check: V0 = {V0:.8f}, c_s=1/sqrt3={1/S3:.8f}  (flow speed = sound speed at sonic)")
# central density: om=4 pi r^2 a^2 rho ; near center om ~ 1.5 a2 z^2, A->1 so rho_c ~ (1.5 a2)/(4 pi z^2 r^2)
# with r=e^{x-s}=z e^{-s}: 4 pi r^2 rho = om/A -> at center finite central 4pi rho r^2/z^2 relation.
print(f"  central density param a2* = {a:.10f} (sets the homothetic amplitude; nc={1.0} is x-gauge)")
