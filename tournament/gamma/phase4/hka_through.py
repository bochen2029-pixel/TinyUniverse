# hka_through.py — integrate the EC background THROUGH the sonic point into the outer region, to
# confirm the EC selection: exactly ONE zero of V (md line 77) and correct outer behavior.
#
# The sonic point (Dson=0) makes the resolved om',V' blow up. To cross it smoothly, switch to the
# DESINGULARIZED flow dY/dxi (polynomial, finite at Dson=0) in a small window around the crossing,
# then resume the resolved flow. We reuse the constraint-reduced 3D system.
#
# For the EC solution the sonic point is at (A0,N0,om0,V0)=(3/2, 2/sqrt3, 3/4, -1/sqrt3). We verify:
#  - starting from the center at the EC oi, the trajectory reaches this sonic point,
#  - crossing it (via desing flow) and continuing, V has exactly ONE zero in the outer region,
#  - the outer region approaches the EC stationary behavior.
import numpy as np, math
import hka_ec as E

G=E.G; GM1=E.GM1; CS=E.CS

def desing3(Y3):
    """dY3/dxi = d*(N',om',V') with uniform scale d = (fluid-matrix determinant)*om/? Use the same
    numerators as hka_reduced. Smooth at the sonic point."""
    N,om,V=Y3; g=G; oV2=1-V*V
    A=E.A_of(N,om,V)
    Nx=N*(-2+A-(2-g)*om)
    RHS_d=(3*(2-g)/2)*N*V-((2+g)/2)*A*N*V+(2-g)*N*V*om
    RHS_e=(2-g)*(g-1)*N*om+((7*g-6)/2)*N+((2-3*g)/2)*A*N
    # Fluid pair M'.(om_x,V_x) = om*(RHS_d,RHS_e), M'=[[1+NV, g(N+V)om/oV2],[(g-1)(N+V), g(1+NV)om/oV2]].
    a=(1+N*V); b=g*(N+V)/oV2*om; cc=(g-1)*(N+V); d=g*(1+N*V)/oV2*om
    det=a*d-b*cc                             # ∝ om*Dson (finite, ->0 at sonic)
    nOm=(om*RHS_d)*d - b*(om*RHS_e)          # = om_x*det
    nV = a*(om*RHS_e) - (om*RHS_d)*cc        # = V_x *det
    return np.array([Nx*det, nOm, nV]), det

def rk4_res(Y3,h):
    f=lambda Y: np.array(E.rhs3(0,Y))
    k1=f(Y3);k2=f(Y3+0.5*h*k1);k3=f(Y3+0.5*h*k2);k4=f(Y3+h*k3)
    return Y3+(h/6)*(k1+2*k2+2*k3+k4)

def rk4_des(Y3,dxi):
    f=lambda Y: desing3(Y)[0]
    k1=f(Y3);k2=f(Y3+0.5*dxi*k1);k3=f(Y3+0.5*dxi*k2);k4=f(Y3+dxi*k3)
    return Y3+(dxi/6)*(k1+2*k2+2*k3+k4)

def integrate_through(N_inf, oi, x0=-12.0, h=2e-5, Dswitch=1e-3, nmax=2000000, xmax=8.0):
    """Center -> through sonic -> outer. Resolved flow when |Dson|>Dswitch; desingularized flow when
    close to the sonic point. Track x (via dx=det/... for desing leg) and count V zeros."""
    z0=math.exp(x0); Y=E.center_state3(N_inf,oi,z0); x=x0
    Vprev=Y[2]; vz=0
    xs=[x]; Ys=[Y.copy()]
    mode='res'
    for i in range(nmax):
        D=E.Dson(Y[0],Y[2])
        if abs(D)>Dswitch:
            Yn=rk4_res(Y,abs(h)); x+=abs(h)
        else:
            # desingularized step; map dx = (dx/dxi) dxi. dx/dxi = det (the scale used in desing3).
            _,det=desing3(Y)
            dxi=abs(h)/max(abs(det),1e-6)*0.3   # keep dx-ish steps modest
            dxi=math.copysign(min(abs(dxi),0.01), 1.0)
            Yn=rk4_des(Y,dxi)
            _,det2=desing3(Y)
            x+=det2*dxi
        if not np.all(np.isfinite(Yn)): return dict(status='nan',x=x,vz=vz,xs=np.array(xs),Ys=np.array(Ys))
        if Yn[2]*Vprev<0: vz+=1
        Vprev=Yn[2]
        Y=Yn
        if i%20==0: xs.append(x); Ys.append(Y.copy())
        if x>xmax or abs(Y[2])>3 or Y[1]<-1 or Y[1]>60 or Y[0]<1e-6 or Y[0]>1e9:
            break
    return dict(status='done',x=x,vz=vz,Y=Y,xs=np.array(xs),Ys=np.array(Ys))

if __name__=="__main__":
    import sys,time
    N_inf=1.0
    # get a well-converged oi first (V_cross=-1/sqrt3 condition)
    from scipy.optimize import brentq
    def g(oi):
        rr=E.shoot_to_sonic(N_inf,oi,rtol=1e-11,atol=1e-13)
        return (rr['V']+1/math.sqrt(3)) if rr['status']=='sonic' else 1.0
    oistar=brentq(g,0.36,0.39,xtol=1e-12)
    print(f"EC oi* (V_cross=-1/sqrt3): {oistar:.9f}")
    t0=time.time()
    out=integrate_through(N_inf,oistar)
    print(f"through-sonic integration: status={out['status']} x_end={out['x']:.3f} Vzeros={out['vz']}")
    if 'Y' in out:
        N,om,V=out['Y']; print(f"  end state: N={N:.5g} om={om:.5g} V={V:.5f} A={E.A_of(N,om,V):.5f}")
    # print V profile summary
    Ys=out['Ys']; xs=out['xs']
    print("  V profile (x, V) sampled:")
    for j in range(0,len(xs),max(1,len(xs)//20)):
        print(f"    x={xs[j]:+.3f}  N={Ys[j][0]:.4g} om={Ys[j][1]:.4f} V={Ys[j][2]:+.4f}")
    print(f"# {time.time()-t0:.1f}s")
