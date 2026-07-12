# Numerical prototype: integrate the KHA autonomous ODE system from the CENTER outward,
# using the exact transcribed equations, to discover the true sonic-point state and to
# validate the ODE RHS (esp. the sign of omega). This is ground truth; no closed forms.
#
# ODE state y=(N,A,omega,V) as functions of x.
#   A' = A*(1 - A + (2 omega/(1-V^2))(1+V^2/3))          [eq a1]
#   N' = N*(-2 + A - 2 omega/3)                            [eq a2]
#   fluid eqs (b) with _,s=0 give 2x2 for (omega', V'):
#     eq1:  (1+NV) (omega'/omega) + (4(N+V)/(3(1-V^2))) V' = RHS1
#     eq2:  (4V+N+3NV^2)(omega'/omega) + (4(1+V^2+2NV)/(1-V^2)) V' = RHS2
#   where RHS1 = -[ -NV A'/(3A) + 4V N'/3 + 2N(1+4omega/(9(1-V^2))) ]
#         RHS2 = -[ N(1-V^2)A'/A + 4(1+V^2)N' + 2N(1+3V^2) ]
#   (A',N' substituted from a1,a2). Solve 2x2 for (u=omega'/omega, V'); omega'=u*omega.
#
# Center regularity (x->-inf): A=1, V=0. Expand: near center V ~ v1 e^{x}? Standard CSS
# center expansion has V linear in the local radius. We seed at x0 very negative with
# A=1, V=0 and a chosen omega_c; N from the center relation. Then the *value* of omega_c
# is the shoot parameter, adjusted so the solution passes ANALYTICALLY through the sonic
# point (numerators finite where det->0).

import numpy as np

def metric_slopes(N,A,om,V):
    Ap = A*(1 - A + (2*om/(1-V*V))*(1 + V*V/3.0))
    Np = N*(-2 + A - 2.0*om/3.0)
    return Ap, Np

def fluid_slopes(N,A,om,V):
    Ap, Np = metric_slopes(N,A,om,V)
    omV2 = 1 - V*V
    # 2x2 matrix for (u, Vp), u=omega'/omega
    m11 = (1 + N*V)
    m12 = 4*(N+V)/(3*omV2)
    m21 = (4*V + N + 3*N*V*V)
    m22 = 4*(1 + V*V + 2*N*V)/omV2
    rhs1 = -( -N*V*Ap/(3*A) + 4*V*Np/3.0 + 2*N*(1 + 4*om/(9*omV2)) )
    rhs2 = -(  N*omV2*Ap/A   + 4*(1+V*V)*Np + 2*N*(1+3*V*V) )
    det = m11*m22 - m12*m21
    u  = (rhs1*m22 - m12*rhs2)/det
    Vp = (m11*rhs2 - rhs1*m21)/det
    omp = u*om
    return Ap, Np, omp, Vp, det

def Dloc(N,V):  # sonic locus 3N^2V^2 - N^2 + 4NV - V^2 + 3 (=0 at sonic); note this ~ det*const
    return 3*N*N*V*V - N*N + 4*N*V - V*V + 3

def center_N(om_c):
    # at center A=1,V=0: from a2, N'/N = -1 - 2om/3. N stays finite; the center value of N
    # is fixed by requiring N->? Actually N = alpha/a * e^{-x} -> as x->-inf, e^{-x}->+inf.
    # The regular center has alpha/a -> const, so N ~ const*e^{-x} DIVERGES at center?!
    # Re-examine: N = alpha a^{-1} e^{-x}. As x->-inf (center), e^{-x}->inf. For N finite we
    # need alpha/a -> 0 like e^{x}. Hmm. Let's instead integrate from the sonic point.
    return None

# The center is a SINGULAR point too (x->-inf). Cleaner to integrate FROM the sonic point
# (x=0) in BOTH directions, using a regular sonic-point seed. So we need the sonic state.
# Let's find it empirically: at x=0, det=0 must hold => D(N0,V0)=0. Also the numerators
# rhs1,rhs2 must be in the range of the singular matrix (Fredholm) for u,Vp finite.
# Scan (V0) on branch1 N0(V0); for each, compute the Fredholm defect with omega0 chosen to
# kill it, then read A0 from the OTHER solvability... Let's just scan V0,A0,omega0 space for
# the point where BOTH det=0 and the solution is finite, by requiring the two numerators
# (Cramer) vanish proportionally. Instead: pick V0, set N0 from D=0 (branch1). Then require
# rhs1*m22 - m12*rhs2 = 0 AND m11*rhs2 - rhs1*m21 = 0 at det=0 (both numerators vanish) for
# a FINITE slope (0/0). Two equations, two unknowns (A0,omega0). Solve numerically.

from scipy.optimize import fsolve

def sonic_residual(p, V0):
    A0, om0 = p
    N0 = (-2*V0 + np.sqrt(3)*(V0*V0-1))/(3*V0*V0-1)   # branch1
    Ap, Np = metric_slopes(N0,A0,om0,V0)
    omV2 = 1-V0*V0
    m11=(1+N0*V0); m12=4*(N0+V0)/(3*omV2)
    m21=(4*V0+N0+3*N0*V0*V0); m22=4*(1+V0*V0+2*N0*V0)/omV2
    rhs1=-(-N0*V0*Ap/(3*A0)+4*V0*Np/3.0+2*N0*(1+4*om0/(9*omV2)))
    rhs2=-( N0*omV2*Ap/A0+4*(1+V0*V0)*Np+2*N0*(1+3*V0*V0))
    numU = rhs1*m22 - m12*rhs2   # must ->0
    numV = m11*rhs2 - rhs1*m21   # must ->0
    return [numU, numV]

print("scan V0: solve (A0,omega0) so BOTH sonic numerators vanish (regular node):")
print(" V0       N0        A0        omega0     resid")
for V0 in np.linspace(-0.6,0.6,25):
    try:
        N0 = (-2*V0 + np.sqrt(3)*(V0*V0-1))/(3*V0*V0-1)
        sol,info,ier,msg = fsolve(sonic_residual, [1.0, 1.0], args=(V0,), full_output=True)
        r = sonic_residual(sol,V0)
        tag = "OK" if (abs(r[0])+abs(r[1])<1e-8 and sol[1]>0) else ""
        print(f" {V0:+.3f}  {N0:+.5f}  {sol[0]:+.5f}  {sol[1]:+.5f}   {abs(r[0])+abs(r[1]):.1e} {tag}")
    except Exception as e:
        print(f" {V0:+.3f}  err {e}")
