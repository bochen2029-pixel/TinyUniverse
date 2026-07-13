# nr_check2.py — decide whether nr_sonic.bg_series is correct, REFERENCE-FREE, by plugging the
# series back into the ODE and checking the residual Taylor coefficients (should be ~0 up to
# truncation). Then a CLEAN reference: fresh integrations that STOP at each t (never entering the
# singular 0/0 region), so the reference itself is trustworthy.
import numpy as np, math
import nr_sonic as NS, hka_ec as E
from scipy.integrate import solve_ivp

g = 4.0/3.0

def ode_residual(order=12):
    A, N, om, V = NS.bg_series(order)
    K = order+1
    def ev(a, t): return sum(a[k]*t**k for k in range(len(a)))
    # residual of metric N' - N(-2+A-(2-g)om) and fluid M w' - b, as series coeffs
    NS_ = NS
    (M11,M12,M21,M22),(b1,b2) = NS_._fluid_MB(A,N,om,V)
    wpo = NS_.sderiv_shift(om); wpv = NS_.sderiv_shift(V)
    Nprime = NS_.sderiv_shift(N)
    metric_res = Nprime - NS_.smul(N, NS_.sconst(K,-2.0)+A-(2.0-g)*om)
    fluid1 = NS_.smul(M11,wpo)+NS_.smul(M12,wpv)-b1
    fluid2 = NS_.smul(M21,wpo)+NS_.smul(M22,wpv)-b2
    # coefficients 0..order-1 should be ~0 (order-1 because a derivative loses one order)
    print("ODE residual Taylor coeffs (should be ~0 for k=0..order-2):")
    for name,res in [("metricN",metric_res),("fluid_om",fluid1),("fluid_V",fluid2)]:
        mx = max(abs(res[k]) for k in range(order-1))
        print(f"  {name}: max|coeff[0..{order-2}]| = {mx:.2e}   first few: {[f'{res[k].real:+.1e}' for k in range(5)]}")

def clean_ref(order=12):
    A,N,om,V = NS.bg_series(order)
    coeffs=(A,N,om,V)
    r = E.shoot_to_sonic(1.0, 0.375, rtol=1e-12, atol=1e-14); xs=r['x']
    z0=math.exp(-16.0); Y0=E.center_state3(1.0,0.375,z0)
    print(f"\nCLEAN reference (fresh integration stopping AT each x_s+t, tol 1e-13):")
    worst=0.0
    for t in [-0.02,-0.05,-0.08,-0.12,-0.18]:
        sol=solve_ivp(E.rhs3,[-16.0, xs+t], Y0, method='DOP853', rtol=1e-13, atol=1e-15)
        Nr,omr,Vr=sol.y[:,-1]; ref=np.array([E.A_of(Nr,omr,Vr),Nr,omr,Vr])
        ser=NS.eval_bg(coeffs,t).real
        err=np.abs(ser-ref); worst=max(worst,err.max())
        print(f"  t={t:+.3f}: max|series-ref|={err.max():.2e}   ser_om={ser[2]:.6f} ref_om={ref[2]:.6f}")
    print(f"  >>> worst {worst:.2e}")

if __name__=="__main__":
    ode_residual(12)
    clean_ref(12)
