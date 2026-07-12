# diag_bg.py — where does the Stage-B background (center_seed, A2STAR) actually reach the sonic locus?
# Critical consistency: stageB_v2 builds R at the EXACT algebraic state (V0=1/sqrt3). Does the
# integrated background actually pass through THAT state, or a different one?
import numpy as np, math
import css_core as C

S3 = math.sqrt(3.0)
A2STAR = 0.250067905275

def run(nc, z0, h, xmax=4.0):
    Y = C.center_seed(nc, A2STAR, z0); x = math.log(z0)
    prevD = C.Dson(Y[0], Y[3]); prevY = Y; prevx = x
    n = int((xmax - x)/h)
    for i in range(n):
        Y2 = C.rk4_step(Y, h)
        if not all(math.isfinite(v) for v in Y2): return ("blew", x, Y)
        D2 = C.Dson(Y2[0], Y2[3])
        if prevD*D2 < 0 and abs(Y2[3]) > 1e-3:
            fr = prevD/(prevD - D2)
            Yh = tuple(prevY[j] + fr*(Y2[j]-prevY[j]) for j in range(4))
            xh = prevx + fr*h
            return ("sonic", xh, Yh)
        prevD = D2; prevY = Y2; prevx = x+h; Y = Y2; x += h
    return ("nocross", x, Y)

print("exact algebraic sonic state: N0=%.6f A0=1.5 om0=0.75 V0=%.6f" % (2/S3, 1/S3))
print()
print(f"{'nc':>5} {'z0=e^':>7} {'h':>8} {'status':>8} {'x_s':>9} {'N':>9} {'A':>9} {'om':>9} {'V':>9}")
for nc in [1.0, 1.5]:
    for lz in [-11, -14]:
        for h in [1e-4, 2e-5]:
            st, xh, Y = run(nc, math.exp(lz), h)
            print(f"{nc:5.1f} {lz:7d} {h:8.0e} {st:>8} {xh:9.4f} {Y[0]:9.5f} {Y[1]:9.5f} {Y[2]:9.5f} {Y[3]:9.5f}")
