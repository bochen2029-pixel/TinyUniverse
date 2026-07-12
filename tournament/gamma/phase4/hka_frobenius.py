# hka_frobenius.py — analytic Frobenius solutions of Psi' = L(x;kappa) Psi at the sonic point.
#
# L(x_s+t) = R/t + L0 + L1 t + ...  (Laurent, from hka_pert_sonic.L_laurent). Indicial exponents =
# eig(R) = {0,0,0, mu} with mu=1+2kappa (our t=x-x_s convention). The exponent-0 (ANALYTIC) modes are
# Taylor series Psi(t) = sum_{n>=0} a_n t^n. Plugging in: t Psi' = R Psi + t L0 Psi + ...  matching
# powers gives the recursion  (n I - R) a_n = sum_{j>=0} L_j a_{n-1-j}  (with a_{-1}=0 handled).
# For n=0: R a_0 = 0  -> a_0 in ker(R) (the 3D analytic kernel). For n>=1: (nI-R) a_n = RHS (solvable
# since n is not an indicial exponent, as long as n != mu; mu is generally non-integer).
#
# This gives 3 independent analytic solutions (one per ker(R) basis vector). The physical perturbation
# is the combination fixed by: gauge Nbar_p(0)=0 (a_0 has Nbar-comp 0) + solvability/identity (5.15).
# We build the 3 analytic series and expose them; the shoot (hka_beta.py) matches against the
# center-regular subspace.
#
# Eq #s per HKA_beta_equations.md.
import numpy as np, math
import hka_pert_sonic as PS
import hka_pert_core as PC

def analytic_modes(kappa, bgser, order=8, eps=1e-4):
    """Return list of up to 3 analytic Taylor solutions [a_0,a_1,...,a_order] (each a_n is length-4).
    a_0 spans ker(R)."""
    lau=PS.L_laurent(kappa,bgser,order=order+1,eps=eps)
    R=lau[0]; Ls=lau[1:]        # Ls[j] = L_j (coeff of t^j in L's regular part), j=0,1,...
    # kernel of R (3D)
    U,S,Vh=np.linalg.svd(R)
    ker=[np.conj(Vh[i]) for i in range(4) if abs(S[i])<1e-6*max(S.max(),1)]
    modes=[]
    for a0 in ker:
        a=[np.array(a0,complex)]
        ok=True
        for n in range(1,order+1):
            rhs=np.zeros(4,complex)
            for j in range(0,n):           # L_j a_{n-1-j}, j=0..n-1
                if j<len(Ls):
                    rhs+=Ls[j].dot(a[n-1-j])
            Mn=n*np.eye(4)-R
            try:
                a.append(np.linalg.solve(Mn,rhs))
            except np.linalg.LinAlgError:
                ok=False; break
        if ok: modes.append(a)
    return modes, R

def eval_mode(a, t):
    s=np.zeros(4,complex); tp=1.0
    for n in range(len(a)):
        s+=a[n]*tp; tp*=t
    return s

if __name__=="__main__":
    bgser,_,_=PS.bg_series_near_sonic()
    for kappa in [2.81055255, 1.0, complex(2.8,0.1)]:
        modes,R=analytic_modes(kappa,bgser)
        print(f"kappa={kappa}: {len(modes)} analytic modes; ker(R) dim={len(modes)}")
        # show a_0 (the kernel basis) and Nbar component (gauge Nbar_p(0)=0 picks combos)
        for i,a in enumerate(modes):
            a0=a[0]/a[0][np.argmax(np.abs(a[0]))]
            print(f"   mode{i} a_0(Abar,Nbar,ombar,V)={np.round(a0,4)}")
