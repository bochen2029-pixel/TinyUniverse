# diag_center.py — indicial structure of the perturbation ODE at the CENTER (x->-inf).
# We need to know how many center-regular modes there are: that sets whether the eigenvalue
# residual is a single projection (1 regular mode) or a Wronskian (2+ regular modes).
import numpy as np, math
import css_core as C, pert_core as P

S3 = math.sqrt(3.0)
A2STAR = 0.250067905275
NC = 1.0

def bg_slopes(Y):
    N,A,om,V = Y
    return (C.Nx(A,N,om,V), C.Ax(A,N,om,V), C.omx(A,N,om,V), C.Vx(A,N,om,V))

print("Center indicial exponents = eig(L) as x->-inf (z->0), for several k:")
print("(regular modes: Re(eig) > 0, since hp ~ z^eig and z->0)")
for k in [0.0, 0.35699, 0.357, 2.81055, 1.5, complex(2.0,1.0)]:
    z = 1e-9
    Y = C.center_seed(NC, A2STAR, z)
    Nb,Ab,ob,Vb = bg_slopes(Y)
    L = P.Lmat(Y[0],Y[1],Y[2],Y[3], Nb,Ab,ob,Vb, complex(k))
    ev = np.linalg.eigvals(L)
    ev = ev[np.argsort(ev.real)]
    reg = np.sum(ev.real > 1e-6)
    print(f"  k={str(k):>16}:  eig(L) = {np.round(ev,4)}   -> #regular(Re>0)={reg}")

print("\nInterpretation: #regular modes at center = dimension of the shooting space.")
print("If 1 -> single projection residual is correct. If 2 -> need a 2x2 determinant of two")
print("independent center-regular solutions projected on the non-analytic sonic dual.")
