# hka_beta_spec.py — fluid-beta eigenvalue by SPECTRAL (Chebyshev) COLLOCATION of the authoritative
# operator (hka_pert_hka99). The eigenmode system  M_x Psi' - Gmat Psi = -kappa M_s Psi  becomes a
# GENERALIZED matrix eigenvalue problem  A Psi = kappa B Psi  on a Chebyshev grid -> the WHOLE spectrum
# at once (no shooting, no initial guess, no mode collapse). The physical fluid-beta mode is the unique
# real eigenvalue with Re kappa > 0 (paper: kappa=2.81055255, beta=0.35580192).
#
# BCs: center (x_c): Abar=0, V=0 (regular).  sonic (x_m ~ x_s): Nbar=0 (gauge) + sonic solvability
# Cx*r3 - Ax*r4 = 0 with (r3,r4) the fluid rows of (Gmat - kappa M_s)Psi (kappa-linear -> stays in A,B).
import numpy as np, sys
from scipy.linalg import eig
import hka_pert_core as PC, hka_pert_hka99 as H99
PC.Lnum = H99.Lnum
import hka_beta4 as B
B.bg(); B.bg_path(); XS = B.bg()['xs']

def cheb(N):
    x = np.cos(np.pi*np.arange(N+1)/N)
    c = np.hstack([2., np.ones(N-1), 2.])*(-1.)**np.arange(N+1)
    X = np.tile(x, (N+1, 1)).T
    dX = X - X.T
    D = np.outer(c, 1./c)/(dX + np.eye(N+1))
    D = D - np.diag(D.sum(axis=1))
    return D, x

def mats(x, kappa=0.0):
    c = H99.coeffs(B.bg_state(x))
    Mx = np.array([[1,0,0,0],[0,1,0,0],[0,0,c['Ax'],c['Bx']],[0,0,c['Cx'],c['Dx']]], float)
    Gm = np.array([[c['G1'],0,c['G3'],c['G4']],[c['H1'],0,c['H3'],0],
                   [c['E1'],c['E2'],c['E3'],c['E4']],[c['F1'],c['F2'],c['F3'],c['F4']]], float)
    Ms = np.array([[0,0,0,0],[0,0,0,0],[0,0,c['As'],c['Bs']],[0,0,c['Cs'],c['Ds']]], float)
    return Mx, Gm, Ms, c

def idcoef(x):
    """kappa-INDEPENDENT identity coeffs (cN,com,cV) of eq:alg-PP; the Abar coeff is (kappa - A)."""
    A, N, om, V, obx, Vx = [float(np.real(z)) for z in B.bg_state(x)]
    g = 4.0/3.0; oV2 = 1.0 - V*V
    cN = 2.0*g*N*V*om/oV2
    com = 2.0*om*(1.0 + (g-1.0)*V*V + g*N*V)/oV2
    cV = 2.0*g*om*(N*(1.0 + V*V) + 2.0*V)/oV2**2
    return A, cN, com, cV

def solve_spectrum(N=120, xc=-14.0, xm=None):
    if xm is None: xm = XS - 0.01
    Dc, xche = cheb(N)
    xp = xc + (xm - xc)*(xche + 1)/2.0      # map [-1,1]->[xc,xm]; j=0 -> xm (sonic), j=N -> xc (center)
    D = Dc*2.0/(xm - xc)
    n = N + 1
    A = np.zeros((4*n, 4*n)); Bm = np.zeros((4*n, 4*n))
    for j in range(n):
        Mx, Gm, Ms, c = mats(xp[j])
        # rows 1,2,3 = the (Nbar,obar,V) ODE:  (Mx D - Gm) Psi = -kappa Ms Psi
        for m in range(n):
            A[4*j+1:4*j+4, 4*m:4*m+4] = Mx[1:4, :]*D[j, m]
        A[4*j+1:4*j+4, 4*j:4*j+4] -= Gm[1:4, :]
        Bm[4*j+1:4*j+4, 4*j:4*j+4] = -Ms[1:4, :]
        # row 0 = the ALGEBRAIC identity eq:alg-PP (replaces the redundant Abar ODE; builds the
        # identity in so the physical modes are captured):  (kappa - A)Abar + cN Nbar + com obar + cV V = 0
        #   =>  (-A Abar + cN Nbar + com obar + cV V) = kappa (-Abar)
        Aid, cN, com, cV = idcoef(xp[j])
        A[4*j+0, 4*j:4*j+4] = [-Aid, cN, com, cV]
        Bm[4*j+0, 4*j:4*j+4] = [-1.0, 0.0, 0.0, 0.0]
    # --- BCs ---
    def setrow(row, Arow, Brow):
        A[row, :] = 0; Bm[row, :] = 0; A[row, :] = Arow; Bm[row, :] = Brow
    # center j=N: Abar=0, V=0
    jc = N
    r = np.zeros(4*n); r[4*jc+0] = 1; setrow(4*jc+0, r, np.zeros(4*n))
    r = np.zeros(4*n); r[4*jc+3] = 1; setrow(4*jc+3, r, np.zeros(4*n))
    # sonic j=0: gauge Nbar=0 ; solvability Cx*r3 - Ax*r4 = 0 (kappa-linear)
    j0 = 0
    _, Gm0, Ms0, c0 = mats(XS - 1e-7)
    r = np.zeros(4*n); r[4*j0+1] = 1; setrow(4*j0+1, r, np.zeros(4*n))     # gauge
    # solvability replaces the V-row at j0:
    Arow = np.zeros(4*n); Brow = np.zeros(4*n)
    lhsA = c0['Cx']*Gm0[2, :] - c0['Ax']*Gm0[3, :]      # Cx*(Gm row3) - Ax*(Gm row4)
    lhsB = c0['Cx']*Ms0[2, :] - c0['Ax']*Ms0[3, :]
    Arow[4*j0:4*j0+4] = lhsA; Brow[4*j0:4*j0+4] = lhsB   # A - kappa B form: (Gm..) - kappa(Ms..) => B row = +lhsB? see below
    setrow(4*j0+3, Arow, Brow)
    # A Psi = kappa B Psi ; ODE rows: (Mx D - Gm) = -kappa Ms => (Mx D - Gm)Psi = kappa(-Ms)Psi  OK.
    # solvability: [Cx r3 - Ax r4] = 0 with r = (Gm - kappa Ms)Psi => lhsA Psi - kappa lhsB Psi = 0
    #   => A-row = lhsA, B-row = lhsB.  (already set)
    w, V = eig(A, Bm)
    return w, V, xp

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 120
    w, V, xp = solve_spectrum(N)
    w = w[np.isfinite(w)]
    real = w[np.abs(w.imag) < 1e-6*np.maximum(np.abs(w.real), 1)].real
    real = np.sort(real[np.abs(real) < 50])
    pos = real[real > 0.05]
    print(f"N={N}: finite real eigenvalues in (0,50):")
    for k in pos[:25]:
        tag = "  <== RELEVANT (target 2.81055)" if abs(k-2.81055) < 0.2 else ("  (gauge ~0.357)" if abs(k-0.357)<0.05 else "")
        print(f"   kappa = {k:.6f}   beta=1/k = {1.0/k:.6f}{tag}")
