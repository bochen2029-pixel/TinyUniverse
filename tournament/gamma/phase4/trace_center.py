# trace_center.py — integrate the CORRECT-seed outgoing trajectory from the regular center and see
# exactly how/where it meets the singular locus det=0 (Dson=0). Scan the shooting parameter a2.
import css_core as C
import math

def integrate(nc, a2, z0=math.exp(-14), xmax=6.0, h=1e-4):
    Y = C.center_seed(nc, a2, z0)
    x = math.log(z0)
    prevD = C.Dson(Y[0], Y[3])
    prevY = Y; prevx = x
    hit = None
    n = int((xmax - x)/h)
    for i in range(n):
        try:
            Y2 = C.rk4_step(Y, h)
        except (ValueError, ZeroDivisionError, OverflowError):
            return ('blow', x, Y, hit)
        x2 = x + h
        if not all(math.isfinite(v) for v in Y2):
            return ('blow', x, Y, hit)
        D2 = C.Dson(Y2[0], Y2[3])
        # detect crossing of det=0 (i.e. Dson sign change) once we're off the center (V not tiny)
        if hit is None and prevD*D2 < 0 and abs(Y2[3]) > 1e-3:
            fr = prevD/(prevD - D2)
            Yh = tuple(prevY[k] + fr*(Y2[k]-prevY[k]) for k in range(4))
            hit = (prevx + fr*h, Yh)
        prevD = D2; prevY = Y2; prevx = x2
        Y = Y2; x = x2
        # stop shortly after the sonic crossing to inspect
        if hit is not None and x > hit[0] + 0.02:
            return ('sonic', x, Y, hit)
        if Y[1] <= 0 or Y[1] > 1e4 or abs(Y[3]) >= 0.999:
            return ('edge', x, Y, hit)
    return ('end', x, Y, hit)

print("Outgoing trajectory from regular center (correct seed v1=1/(2nc), w2=1.5 a2).")
print("Scan a2 (central density). Report sonic crossing (Dson=0) state:")
print(f"{'nc':>5} {'a2':>7} | {'status':>7} {'x_son':>8} {'N':>8} {'A':>8} {'om':>8} {'V':>8}  Dresid")
for nc in [1.0]:
    for a2 in [0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0]:
        st, x, Y, hit = integrate(nc, a2)
        if hit:
            xs, Yh = hit
            print(f"{nc:5.1f} {a2:7.3f} | {st:>7} {xs:8.3f} {Yh[0]:8.4f} {Yh[1]:8.4f} {Yh[2]:8.4f} {Yh[3]:8.4f}  {C.Dson(Yh[0],Yh[3]):+.1e}")
        else:
            print(f"{nc:5.1f} {a2:7.3f} | {st:>7}   no-sonic  end x={x:.2f} N={Y[0]:.3f} A={Y[1]:.3f} om={Y[2]:.4f} V={Y[3]:.4f}")
