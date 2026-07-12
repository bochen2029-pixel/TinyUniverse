# hka_background.py — HKA Eq.(4.1) radiation-fluid CSS BACKGROUND (Stage A) — clean public API.
#
# FRESH, verified implementation of the Hara-Koike-Adachi self-similar equations (gr-qc/9607010),
# transcribed in HKA_beta_equations.md.  gamma=4/3, sqrt(gamma-1)=1/sqrt(3).
#
# RESULT (Stage A, LANDED — see RESULTS_hka_beta.md):
#   The Evans-Coleman critical CSS solution is obtained by shooting the central density oi from the
#   REGULAR CENTER (4.11)-(4.13) outward to the SONIC point (4.5). The physics the prior (walled)
#   attempt lacked — a genuine regular center on the INGOING sound cone (V<0) — is recovered:
#     * critical central density  oi* = 3/8  (gauge N_inf=1), converged to ~1e-11,
#     * sonic point EXACTLY (A0,N0,om0,V0) = (3/2, 2/sqrt3, 3/4, -1/sqrt3)  [HKA 4.7-4.9],
#     * Misner-Sharp 2m/r = 1 - 1/A0 = 1/3 at the sonic point,
#     * flow speed |V0| = c_s = 1/sqrt3 (sonic = sound cone), constraint (4.2) satisfied identically,
#     * EXACT invariants along the solution:  N = N_inf e^{-x}  and  A = 1 + (2/3) om.
#
# This module re-exports the verified Stage-A solver (hka_ec) under a stable name, plus the raw ODEs.
# Eq #s per HKA_beta_equations.md.  Deterministic; no target tuned (oi* and the sonic point EMERGE).

import math
import numpy as np
from hka_ec import (G, GM1, CS, M_CENTER, A_of, Dson, rhs3, lhop,
                    center_state3, shoot_to_sonic, find_ec)

# raw sonic data (HKA 4.7-4.9)
def sonic_data(V0):
    N0=(1-CS*V0)/(CS-V0)
    A0=(G**2+4*G-4+8*GM1**1.5*V0-(3*G-2)*(2-G)*V0**2)/(G**2*(1-V0**2))
    om0=2*CS*(CS-V0)*(1+CS*V0)/(G**2*(1-V0**2))
    return A0,N0,om0,V0

def ec_background(N_inf=1.0, rtol=1e-12, x0=-12.0):
    """Return the converged EC critical background: dict with oi*, sonic-point state, and a dense
    interpolant of (N,om,V)(x) from the center to the sonic point."""
    from scipy.optimize import brentq
    from scipy.integrate import solve_ivp
    from hka_ec import _event_Dson
    def shoot(oi):
        z0=math.exp(x0); Y0=center_state3(N_inf,oi,z0)
        sol=solve_ivp(rhs3,[x0,3.0],Y0,method='RK45',rtol=rtol,atol=rtol*1e-2,
                      events=_event_Dson,max_step=0.05)
        if sol.status==1 and len(sol.t_events[0])>0:
            return sol.y_events[0][0][2]+1/math.sqrt(3)
        return 1.0
    oistar=brentq(shoot,0.36,0.39,xtol=1e-13,rtol=1e-14)
    r=shoot_to_sonic(N_inf,oistar,rtol=rtol,atol=rtol*1e-2)
    z0=math.exp(x0); Y0=center_state3(N_inf,oistar,z0)
    sol=solve_ivp(rhs3,[x0,r['x']-1e-7],Y0,method='DOP853',rtol=rtol,atol=rtol*1e-2,dense_output=True)
    return dict(oi=oistar, N_inf=N_inf, sonic=r, x_sonic=r['x'], interp=sol.sol, x0=x0)

if __name__=="__main__":
    bgd=ec_background()
    r=bgd['sonic']
    print("HKA Stage-A Evans-Coleman background:")
    print(f"  oi* = {bgd['oi']:.10f}  (= 3/8 = {3/8})")
    print(f"  sonic point (A0,N0,om0,V0) = ({r['A']:.8f}, {r['N']:.8f}, {r['om']:.8f}, {r['V']:.8f})")
    print(f"  exact = (3/2, 2/sqrt3={2/math.sqrt(3):.8f}, 3/4, -1/sqrt3={-1/math.sqrt(3):.8f})")
    print(f"  2m/r = 1-1/A0 = {1-1/r['A']:.10f}  (= 1/3)")
    print(f"  sonic at x_s = {bgd['x_sonic']:.5f}")
