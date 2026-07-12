# bg_build.py — build the kappa-INDEPENDENT background trajectory ONCE and pickle it.
# Stores, per RK4 step from center to just before the sonic point, the 4 stage-states + slopes,
# plus Dson at each node. Reused by all Stage-B kappa evaluations (huge speedup).
import numpy as np, math, pickle, sys, time
import css_core as C

S3 = math.sqrt(3.0)
A2STAR = 0.250067905275
NC = 1.0

def bg_slopes(Y):
    N,A,om,V = Y
    return (C.Nx(A,N,om,V), C.Ax(A,N,om,V), C.omx(A,N,om,V), C.Vx(A,N,om,V))

def build(z0, h, dmin):
    """Integrate center->sonic. Stop when Dson first rises above -dmin (i.e. within dmin of sonic,
    on the Dson<0 side). Store ALL steps' stage-states so downstream can pick any stop shell >= dmin."""
    Y = C.center_seed(NC, A2STAR, z0); x0 = math.log(z0); x = x0
    S1=[];S2=[];S3=[];S4=[]; Dnode=[]; xnode=[]
    prevD = C.Dson(Y[0], Y[3])
    n = int((3.0 - x0)/h)
    t0=time.time()
    for i in range(n):
        k1 = bg_slopes(Y)
        Y2 = tuple(Y[j]+0.5*h*k1[j] for j in range(4)); k2 = bg_slopes(Y2)
        Y3 = tuple(Y[j]+0.5*h*k2[j] for j in range(4)); k3 = bg_slopes(Y3)
        Y4 = tuple(Y[j]+h*k3[j] for j in range(4));      k4 = bg_slopes(Y4)
        S1.append(Y + k1); S2.append(Y2 + k2); S3.append(Y3 + k3); S4.append(Y4 + k4)
        Dnode.append(prevD); xnode.append(x)
        Ynext = tuple(Y[j]+(h/6)*(k1[j]+2*k2[j]+2*k3[j]+k4[j]) for j in range(4))
        Dnext = C.Dson(Ynext[0], Ynext[3])
        # stop once we are within dmin of the sonic point (Dnext<0 and -Dnext<=dmin)
        if Dnext < 0 and -Dnext <= dmin and abs(Ynext[3])>0.3:
            Dnode.append(Dnext); xnode.append(x+h)
            S1.append(Ynext + bg_slopes(Ynext))  # sentinel final node (slopes only)
            break
        # safety: if we overshoot to Dnext>0 (crossed sonic) stop too (shouldn't with dmin>0)
        if prevD < 0 and Dnext >= 0:
            Dnode.append(Dnext); xnode.append(x+h); break
        Y = Ynext; prevD = Dnext; x += h
    out = dict(S1=np.array(S1[:len(S2) if len(S1)>len(S2) else len(S1)]),
               S2=np.array(S2), S3=np.array(S3), S4=np.array(S4),
               Dnode=np.array(Dnode), xnode=np.array(xnode),
               h=h, z0=z0, dmin=dmin, nsteps=len(S2))
    # align S1 to nsteps (S1 may have one sentinel extra)
    out['S1_full'] = np.array(S1)
    print(f"  built nsteps={len(S2)} in {time.time()-t0:.1f}s ; final Dson={Dnode[-1]:.3e} at x={xnode[-1]:.4f}")
    print(f"  final state N,A,om,V = {S2[-1][:4] if len(S2) else '??'}")
    return out

if __name__ == "__main__":
    z0 = math.exp(float(sys.argv[1])) if len(sys.argv)>1 else math.exp(-12)
    h  = float(sys.argv[2]) if len(sys.argv)>2 else 5e-5
    dmin = float(sys.argv[3]) if len(sys.argv)>3 else 2e-3
    fn = sys.argv[4] if len(sys.argv)>4 else "bgcache.pkl"
    print(f"building background: z0=e^{math.log(z0):.1f} h={h} dmin={dmin} -> {fn}")
    bg = build(z0, h, dmin)
    pickle.dump(bg, open(fn,"wb"))
    print(f"saved {fn}  ({bg['nsteps']} steps)")
