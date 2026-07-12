# hka_pert_sonic.py — perturbation structure AT the sonic point (regular singular point).
#
# L(x;kappa) is singular at the sonic point because Minv=[[Ax,Bx],[Cx,Dx]]^{-1} blows up (AxDx-BxCx=0
# = the sonic condition 4.5). Near x_s, L(x) ~ R/(x-x_s) + regular, with R the RESIDUE. The indicial
# exponents (eigenvalues of R) classify modes: analytic (exponent 0) vs non-analytic. md line 103-104:
# "indicial exponents {0,0,0,1-2kappa}" -> 3 analytic + 1 non-analytic (exponent 1-2kappa).
# The analytic perturbation at the sonic point is fixed up to scale by kappa via row-proportionality
# + identity (5.15) + gauge Nbar_p(0)=0. We build R and the sonic Frobenius series here.
#
# We compute R by expanding the background about the sonic point x_s. Along the EC background near
# x_s, use the analytic sonic series of (A,N,om,V) in t=(x-x_s). Then L(t) = R/t + L0 + L1 t + ...
# We extract R = lim t*L(t) via the desingularized background series.
#
# Eq #s per HKA_beta_equations.md.
import numpy as np, sympy as sp, math, pickle
import hka_pert_core as PC
import hka_ec as E
import hka_desing as D

G=4.0/3.0; GM1=G-1.0; CS=math.sqrt(GM1)
N0=2/math.sqrt(3); OM0=0.75; V0=-1/math.sqrt(3); A0=1.5

# background sonic series in t=(x-x_s), from the desingularized flow eigenvector (analytic dir).
# We get it from a short high-accuracy integration of the RESOLVED background near the sonic point,
# on BOTH sides, fitting a Taylor series in t. Simpler: use hka_ec background integrated to just
# before the sonic point, then Taylor-fit (A,N,om,V)(t) locally.

def bg_series_near_sonic(order=6, dt=2e-3, npts=25, N_inf=1.0, oi=0.375):
    """Fit Taylor coeffs of (A,N,om,V) in t=(x-x_s) by sampling the EC background near the sonic point
    on the center side (t<0). x_s is where the shoot hits the sonic point."""
    from scipy.integrate import solve_ivp
    # integrate center->sonic, dense, get x_s and a dense sample just inside
    z0=math.exp(-12.0); Y0=E.center_state3(N_inf,oi,z0)
    r=E.shoot_to_sonic(N_inf,oi,rtol=1e-12,atol=1e-14)
    xs=r['x']
    # sample the trajectory on [xs-dt*npts, xs) via dense solve_ivp
    sol=solve_ivp(E.rhs3,[-12.0,xs-1e-9],Y0,method='DOP853',rtol=1e-12,atol=1e-14,dense_output=True)
    ts=np.linspace(-dt*npts, -dt*0.5, npts)
    xsamp=xs+ts
    Ysamp=sol.sol(xsamp)      # 3 x npts (N,om,V)
    Nn,omm,Vv=Ysamp
    Aa=np.array([E.A_of(Nn[i],omm[i],Vv[i]) for i in range(npts)])
    # fit polynomials in t
    coefsA=np.polyfit(ts,Aa,order)[::-1]
    coefsN=np.polyfit(ts,Nn,order)[::-1]
    coefsO=np.polyfit(ts,omm,order)[::-1]
    coefsV=np.polyfit(ts,Vv,order)[::-1]
    return dict(A=coefsA,N=coefsN,om=coefsO,V=coefsV,xs=xs), ts, (Aa,Nn,omm,Vv)

def L_laurent(kappa, bgser, order=6, eps=1e-4):
    """Laurent coefficients of L(x_s+t;kappa) = R/t + L0 + L1 t + ... via DFT on |t|=eps.
    Returns list [R, L0, L1, ...] (R is the t^{-1} coeff). t*L = R + L0 t + L1 t^2 + ...; the p-th
    Taylor coeff of g(t)=t*L(t) is (t*L)_p, and L's Laurent coeff of t^{p-1} = g_p."""
    M=64
    w=np.exp(2j*np.pi*np.arange(M)/M); tt=eps*w
    G=np.empty((M,4,4),complex)
    for m in range(M):
        t=tt[m]
        A=np.polyval(bgser['A'][::-1],t); N=np.polyval(bgser['N'][::-1],t)
        om=np.polyval(bgser['om'][::-1],t); V=np.polyval(bgser['V'][::-1],t)
        omx,Vx=_fluid_slopes_c(A,N,om,V); ombar_x=omx/om
        L=PC.Lnum((A,N,om,V,ombar_x,Vx),complex(kappa))
        G[m]=t*L                              # analytic: g(t)=t*L(t)=R + L0 t + ...
    # Taylor coeffs of g: g_p = (1/M) sum_m G[m] w^{-p m} / eps^p
    coeffs=[]
    for p in range(order+1):
        c=np.zeros((4,4),complex)
        for m in range(M):
            c+=G[m]*np.conj(w[m])**p
        coeffs.append(c/M/(eps**p))
    # L Laurent: L = sum_p g_p t^{p-1}; so R=g_0 (t^-1 coeff), L0=g_1, L1=g_2,...
    return coeffs   # coeffs[0]=R, coeffs[1]=L0, ...

def residue_matrix(kappa, bgser, eps=1e-4):
    return L_laurent(kappa,bgser,order=1,eps=eps)[0]

def _fluid_slopes_c(A,N,om,V):
    g=G; oV2=1-V*V
    RHS_d=(3*(2-g)/2)*N*V-((2+g)/2)*A*N*V+(2-g)*N*V*om
    RHS_e=(2-g)*(g-1)*N*om+((7*g-6)/2)*N+((2-3*g)/2)*A*N
    a=(1+N*V)/om; b=g*(N+V)/oV2; c=(g-1)*(N+V)/om; d=g*(1+N*V)/oV2
    det=a*d-b*c
    return (d*RHS_d-b*RHS_e)/det,(a*RHS_e-c*RHS_d)/det

if __name__=="__main__":
    bgser,ts,samp=bg_series_near_sonic()
    print(f"sonic point x_s={bgser['xs']:.5f}")
    print("bg series at sonic (t^0 coeffs A,N,om,V):", bgser['A'][0],bgser['N'][0],bgser['om'][0],bgser['V'][0])
    print("  (exact 1.5, 1.1547, 0.75, -0.57735)")
    print("bg series slopes (t^1 coeffs):", bgser['A'][1],bgser['N'][1],bgser['om'][1],bgser['V'][1])
    for kappa in [2.81055255, 1.0, 0.35699]:
        R=residue_matrix(kappa,bgser)
        ev=np.linalg.eigvals(R)
        print(f"\nkappa={kappa}: residue R indicial exponents (eig R) = {np.sort_complex(ev).round(5)}")
        print(f"   expect {{0,0,0,1-2kappa}}: 1-2kappa={1-2*kappa:.5f}")
