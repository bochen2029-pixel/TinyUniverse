# POSITIVE CONTROL: the GAUGE mode (5.20) is an exact perturbation solution regular at both ends.
# In the sonic gauge Nbar_p(x_s)=0 it has kappa=kbar=0.35699 (md fn15). So if we take v_center = the
# GAUGE-mode direction (background log-derivative vector) and match against the sonic analytic subspace,
# the rejection ||c - QQ^H c|| must -> 0 at kappa=0.35699. This validates the center-integration + sonic
# Frobenius + match machinery independently of the physical-mode search.
import numpy as np, math
import hka_beta4 as B
import hka_frobenius as FR
import hka_pert_sonic as PS

def rej_gauge(kappa, xc=-13.0, xm=None, h=1e-3, order=8):
    B.bg(); B.bg_path()
    xs=B.bg()['xs']
    if xm is None: xm=xs-0.03
    vc=B.gauge_dir(xc)                    # use the GAUGE direction as the center seed (NOT projected out)
    c=B.integ(vc,complex(kappa),xc,xm,h=h)
    if c is None: return None
    modes,R=FR.analytic_modes(complex(kappa),B.bg(),order=order)
    if len(modes)<3: return None
    tm=xm-xs
    scols=[B.eval_series(a,tm) for a in modes]
    Q,_=np.linalg.qr(np.column_stack(scols)); Q=Q[:,:3]
    rej=c-Q.dot(Q.conj().T.dot(c))
    return np.linalg.norm(rej)

if __name__=="__main__":
    print("POSITIVE CONTROL: rejection of the GAUGE-mode center seed vs sonic analytic subspace.")
    print("Expect a ZERO at kappa=kbar~0.35699 (the sonic-gauge gauge mode).")
    print(f"{'kappa':>9} {'||rej||':>12}")
    for k in np.arange(0.20,0.55,0.02):
        r=rej_gauge(round(k,4))
        print(f"{k:9.4f} {r:12.4e}")
    # fine near 0.357
    print("fine:")
    for k in np.arange(0.34,0.38,0.005):
        r=rej_gauge(round(k,5))
        print(f"{k:9.5f} {r:12.4e}")
