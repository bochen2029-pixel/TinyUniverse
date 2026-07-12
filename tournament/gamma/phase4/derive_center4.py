# Determine the leading center behavior by directly balancing the ODEs at small z=e^x.
# Try N = n_c / z (leading), A = 1 + a2 z^p, V = v1 z^q, omega = w z^m and find (p,q,m) and
# the coefficient RELATIONS by plugging into the 4 ODEs and matching lowest orders NUMERICALLY.
import numpy as np

def metric_slopes(N,A,om,V):
    Ap=A*(1-A+(2*om/(1-V*V))*(1+V*V/3.0)); Np=N*(-2+A-2.0*om/3.0); return Ap,Np
def full_slopes(N,A,om,V):
    Ap,Np=metric_slopes(N,A,om,V); omV2=1-V*V
    m11=(1+N*V);m12=4*(N+V)/(3*omV2);m21=(4*V+N+3*N*V*V);m22=4*(1+V*V+2*N*V)/omV2
    rhs1=-(-N*V*Ap/(3*A)+4*V*Np/3.0+2*N*(1+4*om/(9*omV2)))
    rhs2=-( N*omV2*Ap/A+4*(1+V*V)*Np+2*N*(1+3*V*V))
    det=m11*m22-m12*m21
    return Np,Ap,(rhs1*m22-m12*rhs2)/det*om/om if False else None  # placeholder

# Numerically probe the leading behavior: assume at small z, N=nc/z. Then N'=z dN/dz=-nc/z=-N.
# ODE N'/N=-2+A-2om/3 -> -1 = -2+A-2om/3 -> A - 2om/3 = 1. At center A->1 => om->0 consistent (om~o(1)).
# So leading: A=1+..., om-> (3/2)(A-1) at leading? A-2om/3=1 => om=(3/2)(A-1). Interesting: om tied to A!
# Let's test the V equation leading order. Suppose V=v1 z^q. Compute Vp and see the power balance.
# Do it numerically: pick tiny z, set A=1+a2 z^2, om=(3/2)(A-1)=1.5 a2 z^2, V=v1 z^q, N=nc/z with
# nc TBD, and require the ODEs dx(A)=Ap, dx(V)=Vp be consistent to leading order for some (nc,a2,v1,q).

import sympy as sp
z=sp.symbols('z',positive=True)
nc,a2,v1,q,a3,v2,w2b=sp.symbols('nc a2 v1 q a3 v2 w2b',positive=False)
# ansatz to O(z^2) with unknown powers fixed to integers guess q=1
for qguess in [1,2]:
    N=nc/z; A=1+a2*z**2; V=v1*z**qguess; om=sp.Rational(3,2)*a2*z**2  # om tied to A at leading
    def dx(F): return z*sp.diff(F,z)
    Ap=A*(1-A+(2*om/(1-V**2))*(1+V**2/3)); Np=N*(-2+A-sp.Rational(2,3)*om)
    omV2=1-V**2
    m11=(1+N*V);m12=4*(N+V)/(3*omV2);m21=(4*V+N+3*N*V**2);m22=4*(1+V**2+2*N*V)/omV2
    rhs1=-(-N*V*Ap/(3*A)+sp.Rational(4,3)*V*Np+2*N*(1+4*om/(9*omV2)))
    rhs2=-( N*omV2*Ap/A+4*(1+V**2)*Np+2*N*(1+3*V**2))
    det=sp.simplify(m11*m22-m12*m21)
    Vp=sp.simplify((m11*rhs2-rhs1*m21)/det)
    # leading order of dx(V)-Vp
    r=sp.series(dx(V)-Vp, z, 0, qguess+2).removeO()
    print(f"q={qguess}: leading V-eq residual series =", sp.simplify(r))
    # leading order of dx(A)-Ap
    rA=sp.series(dx(A)-Ap, z, 0, 4).removeO()
    print(f"       A-eq residual series =", sp.simplify(rA))
    print()
