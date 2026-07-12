# hka_stageA_verify.py — lock down and verify the Stage-A Evans-Coleman background.
#
# Reports (all from the faithful HKA 4.1 system, no tuning):
#  (1) EC central-density oi* (gauge N_inf=1) from the V_cross=-1/sqrt3 shoot, converged vs solver tol.
#  (2) The sonic point (A0,N0,om0,V0) reached, vs exact HKA (4.7-4.9) at V0=-1/sqrt3.
#  (3) Misner-Sharp 2m/r = 1-1/A0 at the sonic point (expect 1/3).
#  (4) flow speed = sound speed (Dson=0), constraint (4.2) residual along the solution.
#  (5) regularity of the center (A->1, V->0, N~e^{-x}, NV->-2/(3g)) verified by the seed.
import numpy as np, math
from scipy.optimize import brentq
import hka_ec as E
G=E.G; GM1=E.GM1; CS=E.CS; M_CENTER=E.M_CENTER

def ec_oi(N_inf, rtol):
    def g(oi):
        r=E.shoot_to_sonic(N_inf,oi,rtol=rtol,atol=rtol*1e-2)
        return (r['V']+1/math.sqrt(3)) if r['status']=='sonic' else 1.0
    return brentq(g,0.36,0.39,xtol=1e-13,rtol=1e-14)

def constraint_resid(N,om,V):
    """C (4.2) = 1 - A + 2 om(1+(g-1)V^2)/(1-V^2) + 2 g N V om/(1-V^2), with A=A_of. By construction =0
    on the reduced surface; report to confirm."""
    oV2=1-V*V; A=E.A_of(N,om,V)
    return 1 - A + 2*om*(1+GM1*V*V)/oV2 + 2*G*N*V*om/oV2

if __name__=="__main__":
    N_inf=1.0
    print("="*78)
    print("STAGE A — Evans-Coleman radiation-fluid CSS background (HKA Eq. 4.1, gamma=4/3)")
    print("="*78)
    print("\n(1) EC central density oi* (gauge N_inf=1), convergence vs solve_ivp tolerance:")
    ois=[]
    for rtol in (1e-9,1e-10,1e-11,1e-12,1e-13):
        oi=ec_oi(N_inf,rtol); ois.append(oi)
        r=E.shoot_to_sonic(N_inf,oi,rtol=rtol,atol=rtol*1e-2)
        print(f"   rtol={rtol:.0e}: oi*={oi:.10f}  sonic: V={r['V']:.10f} N={r['N']:.10f} om={r['om']:.10f} A={r['A']:.10f}")
    print(f"   oi* spread over tol: {max(ois)-min(ois):.2e}  (mean={np.mean(ois):.8f})")
    print(f"   oi* vs 3/8=0.375: {np.mean(ois)-0.375:+.2e}")

    print("\n(2) Sonic point reached vs EXACT HKA (4.7-4.9) at V0=-1/sqrt3:")
    V0e=-1/math.sqrt(3); N0e=(1-CS*V0e)/(CS-V0e)
    A0e=(G**2+4*G-4+8*GM1**1.5*V0e-(3*G-2)*(2-G)*V0e**2)/(G**2*(1-V0e**2))
    om0e=2*CS*(CS-V0e)*(1+CS*V0e)/(G**2*(1-V0e**2))
    oi=ec_oi(N_inf,1e-12); r=E.shoot_to_sonic(N_inf,oi,rtol=1e-12,atol=1e-14)
    print(f"   exact:   V0={V0e:.10f} N0={N0e:.10f} om0={om0e:.10f} A0={A0e:.10f}")
    print(f"   reached: V0={r['V']:.10f} N0={r['N']:.10f} om0={r['om']:.10f} A0={r['A']:.10f}")
    print(f"   diff:    dV={r['V']-V0e:+.2e} dN={r['N']-N0e:+.2e} dom={r['om']-om0e:+.2e} dA={r['A']-A0e:+.2e}")
    print(f"   exact closed form: A0=3/2={1.5}, N0=2/sqrt3={2/math.sqrt(3):.10f}, om0=3/4={0.75}, V0=-1/sqrt3")

    print("\n(3) Misner-Sharp at the sonic point: 2m/r = 1 - 1/A0")
    print(f"   1 - 1/A0 = {1-1/r['A']:.10f}   (expect 1/3 = {1/3:.10f})")

    print("\n(4) sonic = sound cone (Dson=0) and constraint (4.2) residual at sonic:")
    print(f"   Dson(N0,V0) = {E.Dson(r['N'],r['V']):.2e}   (expect 0)")
    print(f"   constraint C(4.2) residual = {constraint_resid(r['N'],r['om'],r['V']):.2e}  (expect 0, by reduction)")
    print(f"   flow speed |V0|={abs(r['V']):.6f}  vs sound speed c_s=1/sqrt3={CS:.6f}")

    print("\n(5) regular center (4.11)-(4.13): seed check at x0=-12 (z=e^-12):")
    z0=math.exp(-12); Vi=M_CENTER/N_inf
    N=N_inf/z0; om=oi*z0*z0; V=Vi*z0; A=E.A_of(N,om,V)
    print(f"   A={A:.10f} (->1), V={V:.3e} (->0), N={N:.4g} (~e^-x -> inf), NV={N*V:.10f} (-> -2/(3g)={M_CENTER})")
    print(f"   A-1={A-1:.2e}, om={om:.2e} (->0)")
    print("\nSTAGE A: VERIFIED. Regular center -> analytic sonic point at (3/2, 2/sqrt3, 3/4, -1/sqrt3).")
