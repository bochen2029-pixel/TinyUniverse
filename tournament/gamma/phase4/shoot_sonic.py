# shoot_sonic.py — the analytic-passage shooting condition.
# A generic outgoing trajectory hits Dson=0 with numV,numOm != 0, so V'=numV/det diverges (turns
# singular). The CRITICAL a2* reaches the sonic point along a separatrix: numV=0 AND numOm=0 there.
# Since det->0 too, the natural finite residual is  R(a2) := (numV/|grad terms|) at the crossing, or
# more robustly the value of V' extrapolated. We monitor numV and numOm at the crossing vs a2 and
# find the a2 that zeroes them (they should share the zero if the passage is genuinely analytic).
import css_core as C
import math

def cross_residual(nc, a2, z0=math.exp(-14), xmax=4.0, h=5e-5):
    Y = C.center_seed(nc, a2, z0)
    x = math.log(z0)
    prevD = C.Dson(Y[0], Y[3]); prevY = Y; prevx = x
    n = int((xmax - x)/h)
    for i in range(n):
        try:
            Y2 = C.rk4_step(Y, h)
        except (ValueError, ZeroDivisionError, OverflowError):
            return None
        if not all(math.isfinite(v) for v in Y2):
            return None
        D2 = C.Dson(Y2[0], Y2[3])
        if prevD*D2 < 0 and abs(Y2[3]) > 1e-3:
            fr = prevD/(prevD - D2)
            Yh = tuple(prevY[k] + fr*(Y2[k]-prevY[k]) for k in range(4))
            N, A, om, V = Yh[0], Yh[1], Yh[2], Yh[3]
            nV = C.numV(A, N, om, V); nO = C.numOm(A, N, om, V)
            return (prevx + fr*h, Yh, nV, nO)
        prevD = D2; prevY = Y2; prevx = x + h
        Y = Y2; x += h
    return None

print("Residual of analytic passage (numV, numOm at the sonic crossing) vs a2:")
print(f"{'a2':>8} {'x_son':>8} {'V0':>8} {'N0':>8} {'A0':>8} {'om0':>8} {'numV':>12} {'numOm':>12}")
nc = 1.0
prev = None
for a2 in [0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.45,0.50,0.60,0.70,0.80,1.0,1.2,1.5,2.0,3.0]:
    r = cross_residual(nc, a2)
    if r is None:
        print(f"{a2:8.3f}   (no clean crossing)")
        continue
    xs, Yh, nV, nO = r
    print(f"{a2:8.3f} {xs:8.3f} {Yh[3]:8.4f} {Yh[0]:8.4f} {Yh[1]:8.4f} {Yh[2]:8.4f} {nV:12.4e} {nO:12.4e}")
