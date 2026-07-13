# hka_pert_hka99.py — the AUTHORITATIVE Stage-B perturbation operator, transcribed VERBATIM from
# the HKA long paper (gr-qc/9607010, rflanl.tex, eq:EOM-var / eq:EOM-eigenmodes).
#
# The paper writes:   M_x . Psi' = (Gmat - kappa M_s) . Psi        (eq:EOM-eigenmodes)
#   M_x  = [[1,0,0,0],[0,1,0,0],[0,0,Ax,Bx],[0,0,Cx,Dx]]     (NOT identity! fluid block = [[Ax,Bx],[Cx,Dx]])
#   Gmat = [[G1,0,G3,G4],[H1,0,H3,0],[E1,E2,E3,E4],[F1,F2,F3,F4]]
#   M_s  = [[0,0,0,0],[0,0,0,0],[0,0,As,Bs],[0,0,Cs,Ds]]
# => L(kappa) = M_x^{-1} (Gmat - kappa M_s).   Psi = (Abar_p, Nbar_p, obar_p, V_p).
#
# The coefficients (5.5-5.10) MATCH the prior transcription verbatim; the wall was the ASSEMBLY
# (M_x != I). Gauge mode (eq:CT.var): Psi_g = (Abar_x, Nbar_x + kbar, obar_x, V_x); exact for all kbar.
import numpy as np, math
G = 4.0/3.0

def coeffs(fld):
    """fld = (A,N,om,V, obar_x, V_x) as PC.bg_fields returns (obar_x = om'/om, V_x = V')."""
    A, N, om, V, obx, Vx = [float(np.real(z)) for z in fld]
    g = G; oV2 = 1.0 - V*V
    As, Bs, Cs, Ds = 1.0, g*V/oV2, (g-1.0)*V, g/oV2
    Ax, Bx, Cx, Dx = 1.0+N*V, g*(N+V)/oV2, (g-1.0)*(N+V), g*(1.0+N*V)/oV2
    E1 = -((g+2.0)/2.0)*A*N*V
    E2 = ((6.0-3.0*g)/2.0)*N*V - ((2.0+g)/2.0)*A*N*V + (2.0-g)*om*N*V - N*V*obx - g*N*Vx/oV2
    E3 = (2.0-g)*om*N*V
    E4 = ((6.0-3.0*g)/2.0)*N - ((2.0+g)/2.0)*A*N + (2.0-g)*om*N - N*obx - g*(1.0+2.0*N*V+V*V)*Vx/oV2**2
    F1 = ((2.0-3.0*g)/2.0)*A*N
    F2 = (2.0-g)*(g-1.0)*om*N + ((7.0*g-6.0)/2.0)*N + ((2.0-3.0*g)/2.0)*A*N - (g-1.0)*N*obx - g*N*V*Vx/oV2
    F3 = (2.0-g)*(g-1.0)*om*N
    F4 = -(g-1.0)*obx - g*(N+2.0*V+N*V*V)*Vx/oV2**2
    G1, G3, G4 = -A, 2.0*(1.0+(g-1.0)*V*V)*om/oV2, 4.0*g*om*V/oV2**2
    H1, H3 = A, (g-2.0)*om
    return dict(As=As,Bs=Bs,Cs=Cs,Ds=Ds,Ax=Ax,Bx=Bx,Cx=Cx,Dx=Dx,
                E1=E1,E2=E2,E3=E3,E4=E4,F1=F1,F2=F2,F3=F3,F4=F4,G1=G1,G3=G3,G4=G4,H1=H1,H3=H3)

def Lnum(fld, kappa):
    c = coeffs(fld); k = complex(kappa)
    Mx = np.array([[1,0,0,0],[0,1,0,0],[0,0,c['Ax'],c['Bx']],[0,0,c['Cx'],c['Dx']]], complex)
    Gm = np.array([[c['G1'],0,c['G3'],c['G4']],
                   [c['H1'],0,c['H3'],0],
                   [c['E1'],c['E2'],c['E3'],c['E4']],
                   [c['F1'],c['F2'],c['F3'],c['F4']]], complex)
    Ms = np.array([[0,0,0,0],[0,0,0,0],[0,0,c['As'],c['Bs']],[0,0,c['Cs'],c['Ds']]], complex)
    return np.linalg.solve(Mx, Gm - k*Ms)

def _bg_derivs(A, N, om, V, obx, Vx):
    ABx = 1.0 - A + 2.0*om/(1.0 - V*V)*(1.0 + V*V/3.0)   # dlnA/dx (Eq1)
    NBx = -2.0 + A - 2.0*om/3.0                          # dlnN/dx (Eq2)
    return ABx, NBx, obx, Vx

if __name__ == "__main__":
    import hka_beta4 as B
    B.bg(); B.bg_path(); xs = B.bg()['xs']

    # sanity: kappa^1 gauge condition  E2 = As*obar_x + Bs*V_x  and  F2 = Cs*obar_x + Ds*V_x  (holds on bg)
    print("kappa^1 gauge condition (should be ~0 on the background):")
    for x in [-4.0, -1.0, xs-0.1]:
        A,N,om,V,obx,Vx = [float(np.real(z)) for z in B.bg_state(x)]
        c = coeffs(B.bg_state(x))
        r_E2 = c['E2'] - (c['As']*obx + c['Bs']*Vx)
        r_F2 = c['F2'] - (c['Cs']*obx + c['Ds']*Vx)
        print(f"  x={x:+.3f}:  E2-(As obx+Bs Vx)={r_E2:+.2e}   F2-(Cs obx+Ds Vx)={r_F2:+.2e}")

    print("\nGAUGE-MODE EXACTNESS GATE (authoritative HKA99 operator):")
    def gauge(x, kbar):
        A,N,om,V,obx,Vx = [float(np.real(z)) for z in B.bg_state(x)]
        ABx,NBx,OBx,VBx = _bg_derivs(A,N,om,V,obx,Vx)
        return np.array([ABx, NBx+kbar, OBx, VBx], complex)
    xlist = [-8.0,-4.0,-1.0,xs-0.05]
    for kbar in [0.35699, 1.0, 2.81055255]:
        rel_dP, rel_psi = [], []
        for xx in xlist:
            Psi = gauge(xx,kbar)
            dP = (gauge(xx+1e-6,kbar)-gauge(xx-1e-6,kbar))/2e-6
            res = dP - Lnum(B.bg_state(xx),kbar).dot(Psi)
            rel_dP.append(np.linalg.norm(res)/max(np.linalg.norm(dP),1e-12))
            rel_psi.append(np.linalg.norm(res)/np.linalg.norm(Psi))
        print(f"  kbar={kbar:>10}:  |res|/|dP| = {[f'{r:.1e}' for r in rel_dP]}   |res|/|Psi| = {[f'{r:.1e}' for r in rel_psi]}")
