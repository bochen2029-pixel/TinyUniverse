# hka_stageA_convergence.py — the Stage-A convergence RECEIPT: oi* and the sonic point vs (a) solver
# tolerance, (b) center launch depth x0, (c) integrator. Demonstrates the EC background is converged.
import numpy as np, math
from scipy.optimize import brentq
from scipy.integrate import solve_ivp
import hka_ec as E
G=E.G; GM1=E.GM1; CS=E.CS

def ec_oi(N_inf, rtol, x0, method='RK45'):
    def shoot(oi):
        z0=math.exp(x0); Y0=E.center_state3(N_inf,oi,z0)
        sol=solve_ivp(E.rhs3,[x0,3.0],Y0,method=method,rtol=rtol,atol=rtol*1e-2,
                      events=E._event_Dson,max_step=0.05)
        if sol.status==1 and len(sol.t_events[0])>0:
            Ye=sol.y_events[0][0]; return Ye[2]+1/math.sqrt(3)
        return 1.0
    return brentq(shoot,0.36,0.39,xtol=1e-13,rtol=1e-14)

if __name__=="__main__":
    N_inf=1.0
    print("STAGE-A CONVERGENCE RECEIPT")
    print("\n(a) oi* vs solver tolerance (x0=-12, RK45):")
    prev=None
    for rtol in (1e-8,1e-9,1e-10,1e-11,1e-12,1e-13):
        oi=ec_oi(N_inf,rtol,-12.0)
        drow=f"  rtol={rtol:.0e}: oi*={oi:.12f}"
        if prev is not None: drow+=f"   d(prev)={oi-prev:+.2e}"
        print(drow); prev=oi
    print(f"  -> oi* = 3/8 exactly; |oi*-0.375|={abs(oi-0.375):.2e}")

    print("\n(b) oi* vs center launch depth x0 (rtol=1e-12, RK45):")
    for x0 in (-8.0,-10.0,-12.0,-14.0,-16.0):
        oi=ec_oi(N_inf,1e-12,x0)
        print(f"  x0={x0:6.1f}: oi*={oi:.12f}   |oi*-3/8|={abs(oi-0.375):.2e}")

    print("\n(c) oi* vs integrator (x0=-12, rtol=1e-12):")
    for m in ('RK45','DOP853','Radau'):
        try:
            oi=ec_oi(N_inf,1e-12,-12.0,method=m)
            print(f"  {m:8s}: oi*={oi:.12f}   |oi*-3/8|={abs(oi-0.375):.2e}")
        except Exception as ex:
            print(f"  {m:8s}: {ex}")

    print("\n(d) sonic point at converged oi* (should be (3/2, 2/sqrt3, 3/4, -1/sqrt3) exactly):")
    oi=ec_oi(N_inf,1e-12,-12.0)
    r=E.shoot_to_sonic(N_inf,oi,rtol=1e-12,atol=1e-14)
    print(f"  A0={r['A']:.12f} (3/2={1.5})")
    print(f"  N0={r['N']:.12f} (2/sqrt3={2/math.sqrt(3):.12f})")
    print(f"  om0={r['om']:.12f} (3/4={0.75})")
    print(f"  V0={r['V']:.12f} (-1/sqrt3={-1/math.sqrt(3):.12f})")
    print(f"  2m/r=1-1/A0={1-1/r['A']:.12f} (1/3={1/3:.12f})")
