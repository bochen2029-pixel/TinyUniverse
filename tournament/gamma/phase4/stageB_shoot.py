# stageB_shoot.py — the eigenvalue shoot for beta = 1/Re(k0).
# Co-integrate (background H, perturbation hp) from center outward to just before the sonic point.
# Perturbation seeded on the center-REGULAR eigendirection (indicial lambda>=0 -> lambda=+2 mode).
# At the sonic point (regular singular point, detf ~ Dson -> 0), regularity requires the numerator
# numVP -> 0. That complex residual R(k) has zeros at the eigenvalues. Scan/Newton on complex k.
#
# Because hp is linear & homogeneous and the center-regular space is (here) 1-D, the solution is
# fixed up to scale; R(k) := numVP at the sonic crossing (normalized). Discard gauge mode ~0.357;
# the relevant physical mode is ~2.81.
import numpy as np, math, css_core as C, pert_core as P

S3 = math.sqrt(3.0)
A2STAR = 0.250067905275     # from Stage A
NC = 1.0

def bg_slopes(Y):
    N,A,om,V = Y
    return (C.Nx(A,N,om,V), C.Ax(A,N,om,V), C.omx(A,N,om,V), C.Vx(A,N,om,V))

def center_regular_seed(k, z0):
    """Perturbation hp at x0=ln(z0) on the center-regular (lambda=+2) eigendirection of L."""
    Y = C.center_seed(NC, A2STAR, z0)
    Nb,Ab,ob,Vb = bg_slopes(Y)
    L = P.Lmat(Y[0],Y[1],Y[2],Y[3], Nb,Ab,ob,Vb, k)
    ev, V = np.linalg.eig(L)
    # pick eigenvector with the most-positive real eigenvalue (lambda=+2 -> regular z^2)
    idx = int(np.argmax(ev.real))
    vec = V[:, idx]
    # normalize (real gauge: make Ap component real positive if nonzero, else first nonzero)
    for c in [1,0,2,3]:
        if abs(vec[c]) > 1e-9:
            vec = vec/vec[c]; break
    return Y, vec, ev[idx]

def rk4_joint(Y, hp, k, h):
    """One RK4 step of the joint (background 4 real, perturbation 4 complex) system."""
    def deriv(Yr, hpc):
        Nb,Ab,ob,Vb = bg_slopes(Yr)
        dY = (Nb,Ab,ob,Vb)
        L = P.Lmat(Yr[0],Yr[1],Yr[2],Yr[3], Nb,Ab,ob,Vb, k)
        dh = L.dot(hpc)
        return np.array(dY, float), dh
    dY1,dh1 = deriv(Y, hp)
    dY2,dh2 = deriv(tuple(Y[i]+0.5*h*dY1[i] for i in range(4)), hp+0.5*h*dh1)
    dY3,dh3 = deriv(tuple(Y[i]+0.5*h*dY2[i] for i in range(4)), hp+0.5*h*dh2)
    dY4,dh4 = deriv(tuple(Y[i]+h*dY3[i] for i in range(4)), hp+h*dh3)
    Yn = tuple(Y[i]+(h/6)*(dY1[i]+2*dY2[i]+2*dY3[i]+dY4[i]) for i in range(4))
    hpn = hp+(h/6)*(dh1+2*dh2+2*dh3+dh4)
    return Yn, hpn

def residual(k, z0=math.exp(-12), h=2e-4, xmax=2.0):
    """numVP at the sonic crossing, normalized by |hp|, as the analyticity residual."""
    Y, hp, lam = center_regular_seed(k, z0)
    x = math.log(z0)
    prevD = C.Dson(Y[0], Y[3]); prevY = Y; prevhp = hp; prevx = x
    n = int((xmax - x)/h)
    for i in range(n):
        Y, hp = rk4_joint(Y, hp, k, h)
        x += h
        if not all(math.isfinite(v) for v in Y) or not np.all(np.isfinite(hp)):
            return None
        # renormalize hp to avoid overflow (linear system)
        nrm = np.max(np.abs(hp))
        if nrm > 1e6:
            hp = hp/nrm; prevhp = prevhp/nrm
        D2 = C.Dson(Y[0], Y[3])
        if prevD*D2 < 0 and abs(Y[3]) > 1e-3:
            fr = prevD/(prevD - D2)
            Yh = tuple(prevY[j] + fr*(Y[j]-prevY[j]) for j in range(4))
            hph = prevhp + fr*(hp-prevhp)
            Nb,Ab,ob,Vb = bg_slopes(Yh)
            nV = P.numVP(Yh[0],Yh[1],Yh[2],Yh[3], Nb,Ab,ob,Vb, k, hph)
            nO = P.numOmP(Yh[0],Yh[1],Yh[2],Yh[3], Nb,Ab,ob,Vb, k, hph)
            # residual: the (om,V)-fluid analyticity numerators must vanish. Use both, normalized.
            scale = np.max(np.abs(hph)) + 1e-30
            return (nV/scale, nO/scale, Yh, hph)
        prevD = D2; prevY = Y; prevhp = hp; prevx = x
    return None

if __name__ == "__main__":
    print("Residual R(k)=numVP/|hp| at sonic vs real k (scan for zeros):")
    print(f"{'k':>7} {'|numVP|':>12} {'|numOmP|':>12} {'arg(numVP)':>12}")
    for k in np.arange(0.0, 6.01, 0.25):
        r = residual(complex(k))
        if r is None:
            print(f"{k:7.2f}   (no crossing/blow)"); continue
        nV,nO,Yh,hph = r
        print(f"{k:7.2f} {abs(nV):12.4e} {abs(nO):12.4e} {np.angle(nV):12.4f}")
