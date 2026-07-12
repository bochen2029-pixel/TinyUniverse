# stageB_fix.py — CORRECTED eigenvalue extraction for the fluid-CSS critical exponent beta.
#
# Physics (verified, diag_residue.py):
#   At the sonic point the perturbation ODE hp' = L(x)hp is Fuchsian: L ~ R/(x - x_s) + regular
#   (Dson ~ (4/3)(x - x_s), so the pole is in x). Indicial exponents (eigs of R):
#       mu in {0, 0, 0, 1-2k}.
#   Analytic branches: mu=0 (x3). Non-analytic branch: mu = 1-2k, right-eigenvector
#       v_na = [0, 0, 3*sqrt(3)/2, 1]   (k-independent).
#   Near the sonic point:  hp(x) ~ c_a*(analytic) + c_na*(x_s - x)^(1-2k) * v_na.
#   Center-regular solution is the physical one; EIGENVALUE condition = c_na(k) = 0.
#
# Extraction of c_na (finite, pole-free):
#   Left dual of v_na is wL(k); it has a spurious pole at k=3/2. Pole-normalize: wL_n=(2k-3)*wL.
#   wL_n . hp isolates the non-analytic amplitude times the divergent power, so
#       c_na(k) = [wL_n . hp] / [ (x_s - x)^(1-2k) ]   (finite; the divergent power divides out).
#   Track it as a signed complex number vs k. Zeros = eigenvalues (gauge k~0.357, physical k~2.81).
#
# Robust cross-check: at an eigenvalue the center-regular solution stays BOUNDED at x_s; off it,
#   |hp| blows up like (x_s - x)^(1-2k). So log|hp(near x_s)| has sharp MINIMA at eigenvalues.
import numpy as np, math, sympy as sp, pickle
import css_core as C, pert_core as P

S3 = math.sqrt(3.0)
A2STAR = 0.250067905275
NC = 1.0
DSON_SLOPE = 4.0/3.0     # Dson ~ DSON_SLOPE*(x - x_s) near sonic (verified diag_residue)

# --- pole-normalized left dual wL_n(k) = (2k-3)*wL, and the non-analytic right vector v_na ---
_d = pickle.load(open("Lmat.pkl","rb"))
_k = sp.symbols('k')
# from diag_residue (exact, verified):
#   wL   = [ -5/(2k-3), -7*sqrt3/(18k-27), 2*sqrt3*(2k+1)/(9(2k-3)), 1 ]
#   wL_n = (2k-3)*wL = [ -5, -7*sqrt3/9, 2*sqrt3*(2k+1)/9, 2k-3 ]
wL_n = sp.Matrix([[ -5,
                    -7*sp.sqrt(3)/9,
                    2*sp.sqrt(3)*(2*_k+1)/9,
                    2*_k-3 ]])
f_wL_n = sp.lambdify(_k, wL_n, 'numpy')
v_na = np.array([0.0, 0.0, 3*S3/2.0, 1.0], dtype=complex)

def bg_slopes(Y):
    N,A,om,V = Y
    return (C.Nx(A,N,om,V), C.Ax(A,N,om,V), C.omx(A,N,om,V), C.Vx(A,N,om,V))

def center_regular_seed(k, z0):
    Y = C.center_seed(NC, A2STAR, z0)
    Nb,Ab,ob,Vb = bg_slopes(Y)
    L = P.Lmat(Y[0],Y[1],Y[2],Y[3], Nb,Ab,ob,Vb, k)
    ev, Vec = np.linalg.eig(L)
    idx = int(np.argmax(ev.real))          # lambda=+2 regular mode at center
    vec = Vec[:, idx]
    for c in [1,0,2,3]:
        if abs(vec[c]) > 1e-9:
            vec = vec/vec[c]; break
    return Y, vec

def rk4_joint(Y, hp, k, h):
    def deriv(Yr, hpc):
        Nb,Ab,ob,Vb = bg_slopes(Yr); dY = (Nb,Ab,ob,Vb)
        L = P.Lmat(Yr[0],Yr[1],Yr[2],Yr[3], Nb,Ab,ob,Vb, k)
        return np.array(dY,float), L.dot(hpc)
    dY1,dh1 = deriv(Y, hp)
    dY2,dh2 = deriv(tuple(Y[i]+0.5*h*dY1[i] for i in range(4)), hp+0.5*h*dh1)
    dY3,dh3 = deriv(tuple(Y[i]+0.5*h*dY2[i] for i in range(4)), hp+0.5*h*dh2)
    dY4,dh4 = deriv(tuple(Y[i]+h*dY3[i] for i in range(4)), hp+h*dh3)
    Yn = tuple(Y[i]+(h/6)*(dY1[i]+2*dY2[i]+2*dY3[i]+dY4[i]) for i in range(4))
    return Yn, hp+(h/6)*(dh1+2*dh2+2*dh3+dh4)

def integrate_to_sonic(k, dstop, z0, h):
    """Co-integrate to the shell |Dson|=dstop just BEFORE the sonic point. Return (hp, Dson, logabs, xs_minus_x).
    Normalize the perturbation to the ANALYTIC scale by dividing hp by max|analytic-projected| ... instead
    we keep hp seeded at O(1) at the center and DO NOT renormalize destructively (track a logscale factor)."""
    Y, hp = center_regular_seed(k, z0)
    x = math.log(z0)
    logscale = 0.0                 # accumulated log-magnitude removed by renormalization
    prevD = C.Dson(Y[0], Y[3]); prevY = Y; prevhp = hp.copy(); prevx = x
    n = int((3.0 - x)/h)
    for i in range(n):
        Y, hp = rk4_joint(Y, hp, k, h); x += h
        if not (all(math.isfinite(v) for v in Y) and np.all(np.isfinite(hp))):
            return None
        nrm = np.max(np.abs(hp))
        if nrm > 1e12:
            hp = hp/nrm; prevhp = prevhp/nrm; logscale += math.log(nrm)
        D = C.Dson(Y[0], Y[3])
        # stop shell: just BEFORE sonic (Dson<0 side), crossing |Dson|=dstop from below
        if D < 0 and -D <= dstop and -prevD > dstop:
            # linear interp in Dson to the exact shell -D=dstop  => D=-dstop
            fr = (prevD - (-dstop))/(prevD - D)      # prevD>-dstop(more neg? check) fraction to D=-dstop
            fr = min(max(fr,0.0),1.0)
            hph = prevhp + fr*(hp-prevhp)
            Yh = tuple(prevY[j] + fr*(Y[j]-prevY[j]) for j in range(4))
            Dh = C.Dson(Yh[0], Yh[3])
            xs_minus_x = -Dh/DSON_SLOPE     # (x_s - x) = -Dson/slope  (Dson<0 -> positive)
            logabs = math.log(np.max(np.abs(hph))) + logscale
            return hph, Dh, logabs, xs_minus_x, logscale
        prevD = D; prevY = Y; prevhp = hp.copy(); prevx = x
    return None

def c_na(k, dstop=0.02, z0=math.exp(-12), h=5e-5):
    """The non-analytic amplitude c_na(k). Zero at eigenvalues."""
    r = integrate_to_sonic(k, dstop, z0, h)
    if r is None: return None
    hph, Dh, logabs, dxm, logscale = r
    w = np.array(f_wL_n(k), dtype=complex).flatten()
    proj = w.dot(hph)                                  # pole-free projection onto non-analytic dual
    # divide out the divergent power (x_s - x)^(1-2k). dxm>0 real.
    power = dxm**(1.0-2.0*k)                             # (x_s-x)^(1-2k)
    cna = (proj/power) * math.exp(logscale)             # restore renorm magnitude
    return cna, proj, dxm, logabs

if __name__ == "__main__":
    print("=== Coarse kappa scan: boundedness minima + c_na sign (real k) ===")
    print(f"{'k':>7} {'log|hp|@sonic':>14} {'|c_na|':>13} {'Re c_na':>13} {'Im c_na':>13} {'(xs-x)':>10}")
    ks = np.arange(0.1, 3.61, 0.1)
    rows = []
    for k in ks:
        out = c_na(complex(k))
        if out is None:
            print(f"{k:7.2f}     (no shell reached)"); continue
        cna, proj, dxm, logabs = out
        rows.append((k, logabs, cna))
        print(f"{k:7.2f} {logabs:14.4f} {abs(cna):13.4e} {cna.real:13.4e} {cna.imag:13.4e} {dxm:10.4e}")
    # find sign changes in Re c_na
    print("\nSign changes in Re(c_na) (bracket eigenvalues):")
    for i in range(1,len(rows)):
        if rows[i-1][2].real * rows[i][2].real < 0:
            print(f"   between k={rows[i-1][0]:.2f} and k={rows[i][0]:.2f}")
    print("\nBoundedness minima in log|hp| (bracket eigenvalues):")
    for i in range(1,len(rows)-1):
        if rows[i][1] < rows[i-1][1] and rows[i][1] < rows[i+1][1]:
            print(f"   local min at k={rows[i][0]:.2f}  log|hp|={rows[i][1]:.3f}")
