# Refine the EC central density oi in the bracket where lhop(first sonic crossing)=0, and CHECK the
# result against the analytic sonic data: at the converged oi the crossing (N,om,V,A) should equal the
# HKA sonic data (4.7-4.9) for that V, and the trajectory should graze the locus (single tangent).
import numpy as np, math
from scipy.optimize import brentq
import hka_center as C
import hka_shoot_oi as S

G=C.G; GM1=C.GM1; CS=C.CS

def sonic_data(V0):
    N0=(1-CS*V0)/(CS-V0)
    A0=(G**2+4*G-4+8*GM1**1.5*V0-(3*G-2)*(2-G)*V0**2)/(G**2*(1-V0**2))
    om0=2*CS*(CS-V0)*(1+CS*V0)/(G**2*(1-V0**2))
    return A0,N0,om0

def f(oi, N_inf=1.0, h=5e-5):
    val,r=S.F(round(oi,10), N_inf, h=h)
    return val if val is not None else 1e3

if __name__=="__main__":
    import sys
    N_inf=1.0
    for h in (1e-4, 5e-5, 2.5e-5):
        try:
            oi_star=brentq(lambda o: f(o,N_inf,h), 0.35, 0.38, xtol=1e-8, rtol=1e-10)
        except Exception as e:
            print(f"h={h}: brentq failed {e}"); continue
        val,r=S.F(oi_star,N_inf,h=h)
        V=r['V']; A0,N0,om0=sonic_data(V)
        print(f"h={h:.1e}: oi*={oi_star:.7f}  crossing x={r['x']:.4f} V={V:.6f} N={r['N']:.6f} om={r['om']:.6f} A={r['A']:.6f}")
        print(f"          HKA sonic data at that V: A0={A0:.6f} N0={N0:.6f} om0={om0:.6f}")
        print(f"          match: dN={r['N']-N0:+.2e} dom={r['om']-om0:+.2e} dA={r['A']-A0:+.2e}  lhop={val:.2e}")
