# hka_pert_derive.py — derive the Stage-B perturbation operator L(x;kappa) by DIRECTLY linearizing the
# full HKA background PDEs (3.6a,b,d,e), instead of transcribing the paper's precomputed (5.5)-(5.13).
#
# WHY: the transcribed operator failed the gauge-mode exactness gate (RESULTS_hka_beta.md). The clean
# fix: the perturbation equations ARE the linearization of the full self-similar PDEs with d/ds -> kappa.
# The gauge transformation (s,x)->(f^-1(s), x-s+f^-1(s)), N->fdot N is an EXACT symmetry of (3.6), so a
# correctly-linearized L annihilates the pure-gauge mode (5.20) automatically => it passes the gate.
#
# Perturbation: lnA,lnN,lnom perturbed multiplicatively, V additively; mode ~ eps*pert(x)*e^{kappa s}:
#   A -> A(1+eps aa),  N -> N(1+eps nn),  om -> om(1+eps oo),  V -> V+eps vv
#   A_s/A -> eps kappa aa,  om_s/om -> eps kappa oo,  V_s -> eps kappa vv   (background is s-independent)
#   A_x/A -> lAx+eps aap,  om_x/om -> om_x+eps oop,  V_x -> V_x+eps vvp,  N_x/N -> lNx+eps nnp
# Psi = (aa,nn,oo,vv) = (Abar_p, Nbar_p, ombar_p, V_p);  output L is 4x4, AFFINE in kappa.
# Symbols match hka_pert_core.py exactly: A,N,om,V (real), om_x,V_x (real), kappa. Eq #s per HKA_beta_equations.md.
import sympy as sp, pickle

_PKL="hka_pert_L.pkl"

def build():
    g=sp.Rational(4,3)
    A,N,om,V=sp.symbols('A N om V', real=True)
    om_x,V_x=sp.symbols('om_x V_x', real=True)      # background (ln om)_x and V_x
    kappa=sp.symbols('kappa')
    eps=sp.symbols('eps')
    aa,nn,oo,vv=sp.symbols('aa nn oo vv')           # Abar_p, Nbar_p, ombar_p, V_p
    aap,nnp,oop,vvp=sp.symbols('aap nnp oop vvp')   # their x-derivatives
    lAx,lNx=sp.symbols('lAx lNx')                   # background (lnA)_x, (lnN)_x  (drop out)

    # perturbed fields
    Af=A*(1+eps*aa); Nf=N*(1+eps*nn); omf=om*(1+eps*oo); Vf=V+eps*vv
    oV2f=1-Vf**2
    # perturbed derivative quantities
    AxoA=lAx+eps*aap; NxoN=lNx+eps*nnp; oxoom=om_x+eps*oop; Vxf=V_x+eps*vvp
    oso=eps*kappa*oo; Vsf=eps*kappa*vv               # om_s/om, V_s  (A_s/A=eps kappa aa unused: 3.6a has none)

    # ---- full PDE residuals (3.6), LHS - RHS ----
    Ra=AxoA-(1-Af+2*omf*(1+(g-1)*Vf**2)/oV2f)                                   # 3.6a (no s-term)
    Rb=NxoN-(-2+Af-(2-g)*omf)                                                    # 3.6b (no s-term)
    LHSd=oso+g*Vf*Vsf/oV2f+(1+Nf*Vf)*oxoom+g*(Nf+Vf)*Vxf/oV2f                    # 3.6d LHS (s-terms first)
    RHSd=(3*(2-g)/2)*Nf*Vf-((2+g)/2)*Af*Nf*Vf+(2-g)*Nf*Vf*omf
    Rd=LHSd-RHSd
    LHSe=(g-1)*Vf*oso+g*Vsf/oV2f+(g-1)*(Nf+Vf)*oxoom+g*(1+Nf*Vf)*Vxf/oV2f        # 3.6e LHS (s-terms first)
    RHSe=-(g-2)*(g-1)*Nf*omf+((7*g-6)/2)*Nf+((2-3*g)/2)*Af*Nf
    Re=LHSe-RHSe

    lin=lambda R: sp.diff(R,eps).subs(eps,0)         # O(eps) = linearized equation
    La,Lb,Ld,Le=lin(Ra),lin(Rb),lin(Rd),lin(Re)

    aap_s=sp.solve(La,aap)[0]                          # Abar row (kappa-free)
    nnp_s=sp.solve(Lb,nnp)[0]                          # Nbar row (kappa-free)
    ov=sp.solve([Ld,Le],[oop,vvp],dict=True)[0]        # invert the 2x2 -> ombar',V' rows
    rows=[aap_s,nnp_s,ov[oop],ov[vvp]]

    Psi=[aa,nn,oo,vv]
    L=sp.zeros(4,4)
    for i,e in enumerate(rows):
        e=sp.expand(e)
        for j,p in enumerate(Psi):
            L[i,j]=sp.simplify(e.coeff(p))
    aux={eps,aa,nn,oo,vv,aap,nnp,oop,vvp,lAx,lNx}
    leftover=L.free_symbols & aux
    assert not leftover, f"aux symbols leaked into L: {leftover}"
    pickle.dump(dict(L=sp.srepr(L)),open(_PKL,'wb'))
    return L

if __name__=="__main__":
    import time; t0=time.time()
    L=build()
    k=sp.symbols('kappa')
    print(f"derived L in {time.time()-t0:.1f}s")
    print("free symbols:", sorted(map(str,L.free_symbols)))
    print("affine in kappa:", not any(sp.diff(L[i,j],k).has(k) for i in range(4) for j in range(4)))
    # Abar row must match the known-correct (G1,0,G3,G4)
    A,N,om,V=sp.symbols('A N om V', real=True)
    want=[-A,0,2*(1+sp.Rational(1,3)*V**2)*om/(1-V**2),4*sp.Rational(4,3)*om*V/(1-V**2)**2]
    print("Abar row matches (G1,0,G3,G4):",[sp.simplify(L[0,j]-want[j]) for j in range(4)])
