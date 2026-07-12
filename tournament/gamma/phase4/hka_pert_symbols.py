# hka_pert_symbols.py — symbolic assembly of the HKA perturbation operator (5.5)-(5.13).
#
# Builds the linear perturbation system for Psi=(Abar_p, Nbar_p, ombar_p, V_p)^T, mode ~ e^{kappa s}.
# From HKA_beta_equations.md §4:
#   As=1, Bs=g V/(1-V^2), Cs=(g-1)V, Ds=g/(1-V^2)                                  (5.5)
#   Ax=1+NV, Bx=g(N+V)/(1-V^2), Cx=(g-1)(N+V), Dx=g(1+NV)/(1-V^2)                  (5.6)
#   E1..E4, F1..F4                                                                  (5.7),(5.8)
#   G1,G3,G4 ; H1,H3                                                                (5.9),(5.10)
# Matrix form (5.13):
#   Abar_p' = G1 Abar_p + G3 ombar_p + G4 V_p                (first-order, from 4.1a linearization)
#   Nbar_p' = H1 Abar_p + H3 ombar_p                          (from 4.1b linearization)
#   [[Ax,Bx],[Cx,Dx]] . (ombar_p', V_p') = E.Psi + kappa (As ombar_p + Bs V_p ; ...)
#                                        = (E1 Abar+E2 Nbar+E3 ombar+E4 V) + kappa(As ombar_p+Bs V_p)  [row d]
#                                          (F1 Abar+F2 Nbar+F3 ombar+F4 V) + kappa(Cs ombar_p+Ds V_p)  [row e]
# so  (ombar_p', V_p') = [[Ax,Bx],[Cx,Dx]]^{-1} ( E.Psi + kappa Ms.(ombar_p,V_p) ),  Ms=[[As,Bs],[Cs,Ds]].
# The full operator: Psi' = L(x;kappa) Psi,  L AFFINE in kappa (verified: kappa enters only via +kappa Ms).
#
# The background fields (N,A,om,V) and their x-derivatives (ombar_x=om'/om, V_x, etc.) enter the E,F
# coefficients. We lambdify L(x;kappa) = P(fields) + kappa Q(fields) entrywise for fast numeric use.
#
# Aux (5.14) identity for cross-check:
#   kappa Abar_p + Abar_p' = -(2 g N V om/(1-V^2))(Nbar_p+ombar_p) - (2 g N om(1+V^2)/(1-V^2)^2) V_p
# Eq #s per HKA_beta_equations.md.

import sympy as sp, pickle, os

_PKL="hka_pert_L.pkl"

def build():
    g=sp.Rational(4,3)
    A,N,om,V=sp.symbols('A N om V', real=True)
    omx,Vx=sp.symbols('om_x V_x', real=True)     # background ombar_x := d(ln om)/dx = om'/om ; V_x
    # NOTE HKA use ombar_x = (ln om)_x. In (5.7)-(5.8) the terms "N ombar_x", "-(g-1) N ombar_x" etc.
    # use ombar_x = om'/om (log derivative). We pass it as an independent symbol omx := om'/om.
    k=sp.symbols('kappa')
    oV2=1-V**2
    # (5.5)
    As=sp.Integer(1); Bs=g*V/oV2; Cs=(g-1)*V; Ds=g/oV2
    # (5.6)
    Ax=1+N*V; Bx=g*(N+V)/oV2; Cx=(g-1)*(N+V); Dx=g*(1+N*V)/oV2
    # (5.7) E1..E4
    E1=-((g+2)/2)*A*N*V
    E2=((6-3*g)/2)*N*V-((2+g)/2)*A*N*V+(2-g)*om*N*V-N*V*omx-g*N*Vx/oV2
    E3=(2-g)*om*N*V
    E4=((6-3*g)/2)*N-((2+g)/2)*A*N+(2-g)*om*N-N*omx-g*(1+2*N*V+V**2)*Vx/oV2**2
    # (5.8) F1..F4
    F1=((2-3*g)/2)*A*N
    F2=(2-g)*(g-1)*om*N+((7*g-6)/2)*N+((2-3*g)/2)*A*N-(g-1)*N*omx-g*N*V*Vx/oV2
    F3=(2-g)*(g-1)*om*N
    F4=-(g-1)*omx-g*(N+2*V+N*V**2)*Vx/oV2**2
    # (5.9),(5.10)
    G1=-A; G3=2*(1+(g-1)*V**2)*om/oV2; G4=4*g*om*V/oV2**2
    H1=A; H3=(g-2)*om
    # assemble L: Psi=(Abar,Nbar,ombar,V) columns 0..3
    # Abar' = G1 Abar + 0 Nbar + G3 ombar + G4 V
    row_Abar=[G1, 0, G3, G4]
    # Nbar' = H1 Abar + 0 + H3 ombar + 0
    row_Nbar=[H1, 0, H3, 0]
    # (ombar',V') = Minv ( Evec.Psi + k Ms.(ombar,V) )
    Mx=sp.Matrix([[Ax,Bx],[Cx,Dx]])
    Minv=Mx.inv()
    Ms=sp.Matrix([[As,Bs],[Cs,Ds]])
    Evec=sp.Matrix([[E1,E2,E3,E4],[F1,F2,F3,F4]])   # 2x4 acting on Psi
    # RHS as linear in Psi: source_from_E = Evec (2x4). kappa part: k*Ms.(ombar,V) = k*Ms * [Psi_2;Psi_3]
    # Build 2x4 operator: (ombar',V') = Minv*Evec*Psi + k*Minv*Ms*[[0,0,1,0],[0,0,0,1]]*Psi
    sel=sp.Matrix([[0,0,1,0],[0,0,0,1]])
    L_ov = Minv*Evec + k*(Minv*Ms*sel)              # 2x4
    L=sp.Matrix([row_Abar, row_Nbar, list(L_ov.row(0)), list(L_ov.row(1))])  # 4x4
    L=sp.simplify(L)
    data=dict(L=sp.srepr(L))
    pickle.dump(data, open(_PKL,'wb'))
    return L, (A,N,om,V,omx,Vx,k)

if __name__=="__main__":
    import time
    t0=time.time()
    L,syms=build()
    A,N,om,V,omx,Vx,k=syms
    print(f"assembled L(x;kappa) 4x4 in {time.time()-t0:.1f}s")
    # verify affine in kappa: dL/dk should be kappa-free
    dLdk=sp.diff(L,k)
    print("L is affine in kappa:", not any(dLdk[i,j].has(k) for i in range(4) for j in range(4)))
    # check the Abar row is exactly (5.9) (first order, k-free)
    print("Abar' row [G1,0,G3,G4]:", [sp.simplify(L[0,j]-x) for j,x in enumerate([-A, 0, 2*(1+(sp.Rational(1,3))*V**2)*om/(1-V**2), 4*sp.Rational(4,3)*om*V/(1-V**2)**2])])
