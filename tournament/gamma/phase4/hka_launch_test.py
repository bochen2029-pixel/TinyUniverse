# Launch off the sonic point along the desingularized unstable manifold; integrate in xi and watch
# where the trajectory goes (toward regular center A->1,V->0 as x->-inf?). Map x via dx=det dxi.
import numpy as np
import hka_desing as D

def launch(V0, sign_manifold=+1, eps=1e-5, dxi=-2e-4, nmax=400000, xstop=-16.0):
    """Perturb off the sonic point along Re(eigvec) of the leading complex pair; integrate desing
    flow in xi. Track x = integral of det dxi. sign_manifold picks +/- along the real eigenvector;
    dxi sign chosen so x decreases (toward center). Returns arrays."""
    Y0, J, w, Vc = D.sonic_jacobian(V0)
    # leading complex pair (max real part)
    order = np.argsort(-w.real)
    lam = w[order[0]]; vec = Vc[:, order[0]]
    vreal = np.real(vec); vreal = vreal/np.linalg.norm(vreal)
    Y = Y0 + sign_manifold*eps*vreal
    x = 0.0
    xs=[0.0]; Xs=[Y0.copy()]; xis=[0.0]
    xi = 0.0
    for i in range(nmax):
        f = D.desing_rhs(Y)
        Yn = Y + dxi*f
        d = D.det_fluid(Y)
        dx = d*dxi
        x += dx; xi += dxi
        Y = Yn
        if not np.all(np.isfinite(Y)): break
        if i % 200 == 0:
            xs.append(x); Xs.append(Y.copy()); xis.append(xi)
        if x < xstop: break
        # stop if V hits 0 crossing or A diverges
        if abs(Y[0]) > 50 or abs(Y[1]) > 1e6: break
    return np.array(xs), np.array(Xs), np.array(xis)

if __name__ == "__main__":
    import sys
    V0 = -0.25
    if len(sys.argv)>1: V0=float(sys.argv[1])
    for sm in (+1,-1):
        for dxi in (-2e-4, +2e-4):
            xs, Xs, xis = launch(V0, sign_manifold=sm, dxi=dxi)
            if len(Xs)<2:
                print(f"sm={sm:+d} dxi={dxi:+.0e}: died immediately"); continue
            A,N,om,V = Xs[-1]
            print(f"sm={sm:+d} dxi={dxi:+.0e}: reached x={xs[-1]:+.3f} (n={len(xs)}) "
                  f"A={A:.4f} N={N:.4g} om={om:.4f} V={V:.4f}  Vsign_changes="
                  f"{int(np.sum(np.diff(np.sign(Xs[:,3]))!=0))}")
