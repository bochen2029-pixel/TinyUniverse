# hka_desing.py — DESINGULARIZED HKA (4.1) flow, for the sonic-point passage.
#
# The resolved RHS om_x=numOm/det, V_x=numV/det blows up at the sonic point det(M)=0. Introduce a
# new independent variable xi with dx/dxi = det(M) (the fluid-matrix determinant, ∝ Dson). Then
#   dA/dxi = det*A_x,  dN/dxi = det*N_x,  dom/dxi = numOm,  dV/dxi = numV
# is a smooth POLYNOMIAL autonomous system, and the sonic point is a genuine CRITICAL POINT.
# Its Jacobian eigenstructure gives the analytic eigendirections; the regular-center (Evans-Coleman)
# trajectory leaves the sonic point along a specific eigenvector. This is the standard, robust
# desingularization (Maison gr-qc/9504008; Gundlach review) — no fragile series through 0/0.
#
# Derive the polynomial RHS ONCE symbolically, cache lambdified numeric functions. Eq #s per md.

import sympy as sp, numpy as np, pickle, os, math

_PKL = "hka_desing.pkl"
G = 4.0/3.0; GM1 = G-1.0; CS = math.sqrt(GM1)

def _derive():
    g = sp.Rational(4,3)
    A,N,om,V = sp.symbols('A N om V', real=True)
    oV2 = 1 - V**2
    # metric slopes (4.1a,b)
    Ax = A*(1 - A + 2*om*(1+(g-1)*V**2)/oV2)
    Nx = N*(-2 + A - (2-g)*om)
    # fluid pair (4.1d,e): coeff matrix in (om_x, V_x) and RHS (multiply eq by om to clear om_x/om)
    #  (1+NV) om_x + g(N+V)/(1-V^2) * om * V_x = om*RHS_d
    #  (g-1)(N+V) om_x + g(1+NV)/(1-V^2) * om * V_x = om*RHS_e
    RHS_d = (3*(2-g)/2)*N*V - ((2+g)/2)*A*N*V + (2-g)*N*V*om
    RHS_e = (2-g)*(g-1)*N*om + ((7*g-6)/2)*N + ((2-3*g)/2)*A*N
    a11 = (1+N*V);        a12 = g*(N+V)/oV2*om
    a21 = (g-1)*(N+V);    a22 = g*(1+N*V)/oV2*om
    b1 = om*RHS_d;        b2 = om*RHS_e
    det = sp.simplify(a11*a22 - a12*a21)               # ∝ om*Dson/(1-V^2)
    numOm = sp.simplify(b1*a22 - a12*b2)               # Cramer numerator for om_x
    numV = sp.simplify(a11*b2 - b1*a21)                # Cramer numerator for V_x
    # desingularized: multiply by (1-V^2) to also clear the metric denominators uniformly.
    #   dx/dxi = det ; we use DET = det (keeps sign). Metric slopes carry 1/(1-V^2) in Ax only.
    #   To keep everything polynomial, use scale factor S = det (which already has no denominator
    #   except through a12,a22 ~ 1/(1-V^2)). Compute det*Ax etc. and simplify.
    dA = sp.simplify(det*Ax); dN = sp.simplify(det*Nx)
    dom = sp.simplify(numOm); dV = sp.simplify(numV)
    # common denominator check
    syms = (A,N,om,V)
    fns = {k: sp.lambdify(syms, sp.simplify(v), 'numpy') for k,v in
           dict(dA=dA,dN=dN,dom=dom,dV=dV,det=det,Ax=Ax,Nx=Nx,numOm=numOm,numV=numV).items()}
    # also the analytic slopes om_x=numOm/det, V_x=numV/det for cross-check
    data = dict(dA=sp.srepr(dA), dN=sp.srepr(dN), dom=sp.srepr(dom), dV=sp.srepr(dV),
                det=sp.srepr(det), numOm=sp.srepr(numOm), numV=sp.srepr(numV),
                Ax=sp.srepr(Ax), Nx=sp.srepr(Nx))
    return data

def _load():
    if not os.path.exists(_PKL):
        data = _derive()
        pickle.dump(data, open(_PKL, 'wb'))
    else:
        data = pickle.load(open(_PKL, 'rb'))
    A,N,om,V = sp.symbols('A N om V', real=True)
    syms=(A,N,om,V)
    F = {}
    for k in ('dA','dN','dom','dV','det','numOm','numV','Ax','Nx'):
        F[k] = sp.lambdify(syms, sp.sympify(data[k]), 'numpy')
    return F, data

_F, _DATA = _load()

def desing_rhs(Y):
    """dY/dxi (polynomial, smooth at the sonic point). Y=(A,N,om,V)."""
    A,N,om,V = Y
    return np.array([_F['dA'](A,N,om,V), _F['dN'](A,N,om,V),
                     _F['dom'](A,N,om,V), _F['dV'](A,N,om,V)], float)

def det_fluid(Y):
    A,N,om,V = Y; return float(_F['det'](A,N,om,V))

def resolved_rhs(Y):
    """dY/dx = desing_rhs / det (finite away from sonic). For cross-check vs hka_background.rhs."""
    A,N,om,V = Y; d = _F['det'](A,N,om,V)
    return np.array([_F['Ax'](A,N,om,V), _F['Nx'](A,N,om,V),
                     _F['numOm'](A,N,om,V)/d, _F['numV'](A,N,om,V)/d], float)

def sonic_point(V0):
    """(A0,N0,om0,V0) at the sonic point (HKA 4.7-4.9)."""
    N0 = (1-CS*V0)/(CS-V0)
    A0 = (G**2+4*G-4+8*GM1**1.5*V0-(3*G-2)*(2-G)*V0**2)/(G**2*(1-V0**2))
    om0 = 2*CS*(CS-V0)*(1+CS*V0)/(G**2*(1-V0**2))
    return np.array([A0,N0,om0,V0])

def sonic_jacobian(V0, eps=1e-6):
    """Numeric Jacobian of desing_rhs at the sonic point. Eigen-decompose for analytic directions."""
    Y0 = sonic_point(V0)
    J = np.zeros((4,4))
    f0 = desing_rhs(Y0)
    for j in range(4):
        Yp = Y0.copy(); Yp[j]+=eps
        J[:,j] = (desing_rhs(Yp)-f0)/eps
    w, Vc = np.linalg.eig(J)
    return Y0, J, w, Vc


if __name__ == "__main__":
    import sys
    V0 = -0.25
    if len(sys.argv)>1: V0=float(sys.argv[1])
    Y0 = sonic_point(V0)
    print(f"V0={V0}  sonic Y0=(A,N,om,V)={Y0.round(6)}")
    print(f"det_fluid at sonic = {det_fluid(Y0):.3e}  (expect ~0)")
    print(f"desing_rhs at sonic = {desing_rhs(Y0).round(6)}  (expect ~0, it's a critical point)")
    Y0,J,w,Vc = sonic_jacobian(V0)
    print("\nJacobian eigenvalues (of desingularized flow at sonic point):")
    for i in range(4):
        print(f"  lam={w[i].real:+.5f}{w[i].imag:+.5f}j   eigvec={np.real_if_close(Vc[:,i]).round(4)}")
    # cross-check resolved vs hka_background.rhs away from sonic
    import hka_background as HB
    Yt = np.array([1.4, 1.6, 0.4, -0.3])
    print("\ncross-check resolved_rhs vs hka_background.rhs at Y=", Yt)
    print("  desing/det:", resolved_rhs(Yt).round(6))
    print("  hka_background.rhs:", HB.rhs(Yt).round(6))
