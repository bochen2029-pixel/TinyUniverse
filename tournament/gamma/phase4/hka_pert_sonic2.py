# Is L(x;kappa) actually singular at the sonic point along the EC background, or does the singular
# part cancel (row-proportionality / 5.15)? Evaluate L at small t=(x-x_s) using the bg series and
# watch the largest |L entry| and the eigenvalues as t->0.
import numpy as np, math
import hka_pert_core as PC
import hka_pert_sonic as PS
import hka_ec as E

bgser,_,_=PS.bg_series_near_sonic()
xs=bgser['xs']
def fld_at(t):
    A=np.polyval(bgser['A'][::-1],t); N=np.polyval(bgser['N'][::-1],t)
    om=np.polyval(bgser['om'][::-1],t); V=np.polyval(bgser['V'][::-1],t)
    omx,Vx=PS._fluid_slopes_c(A,N,om,V)
    return (A,N,om,V,omx/om,Vx)

for kappa in [2.81055255]:
    print(f"kappa={kappa}: L along EC bg near sonic (t=x-x_s):")
    print(f"{'t':>10} {'max|L|':>12} {'|L*t| max':>12}   eigenvalues(L)")
    for t in [-1e-1,-1e-2,-1e-3,-1e-4,-1e-5]:
        fld=fld_at(t)
        L=PC.Lnum(fld,complex(kappa))
        ev=np.linalg.eigvals(L)
        print(f"{t:10.0e} {np.max(np.abs(L)):12.3e} {np.max(np.abs(L*t)):12.3e}   {np.sort_complex(ev).round(3)}")
    # also: is det([[Ax,Bx],[Cx,Dx]]) ~ t?  Ax=1+NV, Dx=g(1+NV)/oV2, Bx=g(N+V)/oV2, Cx=(g-1)(N+V)
    print("\n check det(Mx)=AxDx-BxCx ~ t (the sonic singularity source):")
    G=PC.G; GM1=PC.GM1
    for t in [-1e-2,-1e-3,-1e-4]:
        A,N,om,V,ox,Vx=fld_at(t); oV2=1-V*V
        Ax=1+N*V; Bx=G*(N+V)/oV2; Cx=(G-1)*(N+V); Dx=G*(1+N*V)/oV2
        detMx=Ax*Dx-Bx*Cx
        print(f"  t={t:.0e}: det(Mx)={detMx:.4e}  (ratio det/t={detMx/t:.4f})")
