# hka_beta_lyap.py — fluid-beta eigenvalue by the LYAPUNOV / time-evolution method (HKA99 sec V's
# "second method", "very regular"): discretize the s-EVOLUTION generator Q and read its spectrum.
# Unlike the eigenmode ODE (singular at the sonic point), the s-evolution is REGULAR everywhere
# (the s-block [[As,Bs],[Cs,Ds]] is invertible: det = g[1-(g-1)V^2]/(1-V^2) > 0).
#
# From eq:EOM-var:  [M_s d_s + M_x d_x] Psi = Gmat Psi,  Psi=(Abar,Nbar,obar,V).
# The metric rows have NO d_s -> spatial CONSTRAINTS:  d_x Abar = G.Psi ,  d_x Nbar = H.Psi .
# The fluid rows evolve:  d_s(obar,V) = Msinv[ (E.Psi, F.Psi) - Mx_fl d_x(obar,V) ],
#   Msinv=[[As,Bs],[Cs,Ds]]^{-1},  Mx_fl=[[Ax,Bx],[Cx,Dx]].
# Eliminate (Abar,Nbar) via the constraints (linear maps of (obar,V)) -> Q acts on (obar,V).
# Eigenmode obar,V ~ e^{kappa s} => Q u = kappa u.  Relevant mode = largest Re kappa (unique) => beta=1/Re.
import numpy as np, sys
import hka_pert_core as PC, hka_pert_hka99 as H99
PC.Lnum = H99.Lnum
import hka_beta4 as B
B.bg(); B.bg_path(); XS = B.bg()['xs']

def cheb(N):
    x = np.cos(np.pi*np.arange(N+1)/N)
    c = np.hstack([2., np.ones(N-1), 2.])*(-1.)**np.arange(N+1)
    X = np.tile(x, (N+1, 1)).T; dX = X - X.T
    D = np.outer(c, 1./c)/(dX + np.eye(N+1)); D -= np.diag(D.sum(axis=1))
    return D, x

def spectrum(N=100, xc=-10.0, xmax=6.0):
    Dc, xche = cheb(N); n = N+1
    xp = xc + (xmax - xc)*(xche + 1)/2.0            # j=0 -> xmax ; j=N -> xc (center)
    D = Dc*2.0/(xmax - xc)
    cof = {kk: np.zeros(n) for kk in ['G1','G3','G4','H1','H3','E1','E2','E3','E4','F1','F2','F3','F4']}
    Msi = np.zeros((n, 2, 2)); MMx = np.zeros((n, 2, 2))
    for j in range(n):
        c = H99.coeffs(B.bg_state(xp[j]))
        for kk in cof: cof[kk][j] = c[kk]
        Ms = np.array([[c['As'], c['Bs']], [c['Cs'], c['Ds']]])
        Mx = np.array([[c['Ax'], c['Bx']], [c['Cx'], c['Dx']]])
        Msi[j] = np.linalg.inv(Ms); MMx[j] = Msi[j].dot(Mx)
    dg = np.diag
    # constraint Abar: (D - diag(G1)) Abar = diag(G3) ob + diag(G4) V ; BC Abar(xc)=0 at j=N
    MA = D - dg(cof['G1']); MA[N, :] = 0; MA[N, N] = 1
    Rom = dg(cof['G3']).copy(); Rom[N, :] = 0
    RV = dg(cof['G4']).copy(); RV[N, :] = 0
    MAi = np.linalg.inv(MA); SAom = MAi.dot(Rom); SAV = MAi.dot(RV)     # Abar = SAom ob + SAV V
    # constraint Nbar: D Nbar = diag(H1) Abar + diag(H3) ob ; BC Nbar(x_s)=0 at nearest grid pt
    js = int(np.argmin(np.abs(xp - XS)))
    MN = D.copy(); MN[js, :] = 0; MN[js, js] = 1
    Ra = dg(cof['H1']).copy(); Ra[js, :] = 0
    Ro = dg(cof['H3']).copy(); Ro[js, :] = 0
    MNi = np.linalg.inv(MN)
    SNom = MNi.dot(Ra.dot(SAom) + Ro); SNV = MNi.dot(Ra.dot(SAV))       # Nbar = SNom ob + SNV V
    # E.Psi, F.Psi as linear maps of (ob,V)
    EPo = dg(cof['E1']).dot(SAom) + dg(cof['E2']).dot(SNom) + dg(cof['E3'])
    EPv = dg(cof['E1']).dot(SAV) + dg(cof['E2']).dot(SNV) + dg(cof['E4'])
    FPo = dg(cof['F1']).dot(SAom) + dg(cof['F2']).dot(SNom) + dg(cof['F3'])
    FPv = dg(cof['F1']).dot(SAV) + dg(cof['F2']).dot(SNV) + dg(cof['F4'])
    # d_s ob = Msi00 EP + Msi01 FP - MMx00 D ob - MMx01 D V ; similarly d_s V
    dob_o = dg(Msi[:,0,0]).dot(EPo) + dg(Msi[:,0,1]).dot(FPo) - dg(MMx[:,0,0]).dot(D)
    dob_v = dg(Msi[:,0,0]).dot(EPv) + dg(Msi[:,0,1]).dot(FPv) - dg(MMx[:,0,1]).dot(D)
    dV_o  = dg(Msi[:,1,0]).dot(EPo) + dg(Msi[:,1,1]).dot(FPo) - dg(MMx[:,1,0]).dot(D)
    dV_v  = dg(Msi[:,1,0]).dot(EPv) + dg(Msi[:,1,1]).dot(FPv) - dg(MMx[:,1,1]).dot(D)
    Q = np.block([[dob_o, dob_v], [dV_o, dV_v]])
    return np.linalg.eigvals(Q)

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 100
    w = spectrum(N)
    w = w[np.isfinite(w)]
    real = w[np.abs(w.imag) < 1e-4*np.maximum(np.abs(w.real), 1)].real
    real = np.sort(real[np.abs(real) < 40])[::-1]
    print(f"N={N}: largest real growth rates kappa (relevant = largest Re>0; expect 2.81055 -> beta 0.35580):")
    for k in real[:12]:
        tag = "  <== RELEVANT?" if abs(k-2.81055) < 0.25 else ""
        print(f"   kappa = {k:+.5f}   beta = 1/k = {1.0/k:+.5f}{tag}" if abs(k) > 1e-6 else f"   kappa = {k:+.5f}")
