# nr_diag.py — diagnose the order-1 degeneracy: is c_1 fixed by the desingularized eigenvector?
import numpy as np, math
import nr_sonic as NS, hka_desing as D

g=4.0/3.0
A0,N0,om0,V0 = 1.5, 2/math.sqrt(3), 0.75, -1/math.sqrt(3)
V0f=-1/math.sqrt(3)

# desingularized Jacobian eigenstructure at the sonic point
Y0,J,w,Vc = D.sonic_jacobian(V0f)
print("desing eigenvalues:", [f"{ww.real:+.4f}{ww.imag:+.4f}j" for ww in w])
# gradient of det_fluid at Y0
Y0a=D.sonic_point(V0f); eps=1e-7
detgrad=np.zeros(4)
for j in range(4):
    Yp=Y0a.copy(); Yp[j]+=eps; Ym=Y0a.copy(); Ym[j]-=eps
    detgrad[j]=(D.det_fluid(Yp)-D.det_fluid(Ym))/(2*eps)
print("grad(det_fluid) at Y0:", detgrad.round(5))

# c_1 = lam v / (grad det . v) for the largest real eigenvalue
idx=np.argmax(w.real)
lam=w[idx].real; v=np.real(Vc[:,idx])
gdv=detgrad.dot(v)
c1_4d = lam*v/gdv
print(f"\nlargest real lam={lam:.6f}, eigvec v(A,N,om,V)={v.round(4)}, grad.v={gdv:.5f}")
print(f"c_1 = lam v/(grad.v) = (A1,N1,om1,V1) = {c1_4d.round(6)}")

# check order-0 fluid: does (om1,V1)_eig satisfy M0 (om1,V1) = b0 ?
K=6
N=NS.sconst(K,N0); om=NS.sconst(K,om0); V=NS.sconst(K,V0)
A=NS._A_series(N,om,V)
(M11,M12,M21,M22),(b1,b2)=NS._fluid_MB(A,N,om,V)
M0=np.array([[M11[0],M12[0]],[M21[0],M22[0]]])
b0=np.array([b1[0],b2[0]])
w1=np.array([c1_4d[2],c1_4d[3]])
print(f"\nM0 = {M0.round(5).tolist()}, det={np.linalg.det(M0):.2e}")
print(f"b0 = {b0.round(5)}")
print(f"M0.w1_eig = {M0.dot(w1).round(6)}   (should equal b0)  resid={np.abs(M0.dot(w1)-b0).max():.2e}")

# metric N1 from the ODE vs eigenvector N-component
N1_metric = N0*(-2+A0-(2-g)*om0)
print(f"\nN1 (metric ODE) = {N1_metric:.6f}   vs eigvec N1 = {c1_4d[1]:.6f}")
# also A1 algebraic vs eigvec
A1_alg = NS._A_series(np.array([N0,c1_4d[1],0,0,0,0],complex),
                      np.array([om0,c1_4d[2],0,0,0,0],complex),
                      np.array([V0,c1_4d[3],0,0,0,0],complex))[1]
print(f"A1 (algebraic from c1) = {A1_alg.real:.6f}   vs eigvec A1 = {c1_4d[0]:.6f}")
