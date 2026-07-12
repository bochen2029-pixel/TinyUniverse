# Look at the shape of lhop(oi) and V_crossing(oi) at the sonic crossing, cleanly, with solve_ivp.
import numpy as np, math
import hka_ec as E

if __name__=="__main__":
    N_inf=1.0
    print(f"{'oi':>8} {'x_cross':>9} {'V_cross':>11} {'V+1/sqrt3':>11} {'N_cross':>10} {'om_cross':>10} {'A_cross':>10} {'lhop':>12}")
    for oi in np.arange(0.30,0.451,0.0075):
        r=E.shoot_to_sonic(N_inf, round(oi,5), rtol=1e-11, atol=1e-13)
        if r['status']=='sonic':
            print(f"{oi:8.4f} {r['x']:9.4f} {r['V']:11.7f} {r['V']+1/math.sqrt(3):11.3e} {r['N']:10.6f} {r['om']:10.6f} {r['A']:10.6f} {r['lhop']:12.3e}")
        else:
            print(f"{oi:8.4f}  status={r['status']}")
