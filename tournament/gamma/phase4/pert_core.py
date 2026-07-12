# pert_core.py — numeric core for the perturbation ODE hp'(x) = L(x;k) . hp.
# Loads Lmat.pkl, lambdifies L (complex k) and the sonic-regularity numerator.
import sympy as sp, pickle, numpy as np
_d = pickle.load(open("Lmat.pkl","rb"))
_A0,_N0,_o0,_V0 = sp.symbols('A0 N0 om0 V0')
_Nx,_Ax,_ox,_Vx = sp.symbols('N_x A_x om_x V_x')
_k = sp.symbols('k')
_Np,_Ap,_opp,_vpp = sp.symbols('Np Ap opp vpp')     # perturbation amplitudes hp=(Np,Ap,omp,Vp)
_L = sp.sympify(_d['Lmat'])
_args = (_N0,_A0,_o0,_V0, _Nx,_Ax,_ox,_Vx, _k)
_fL = sp.lambdify(_args, _L, 'numpy')

# sonic-regularity numerators: opP = numOmP/detf, vpP = numVP/detf, each LINEAR in hp amplitudes.
# detf ~ Dson. At sonic Dson=0, regularity requires numOmP=0 AND numVP=0 on the propagated hp.
_vpP = sp.sympify(_d['vpP']); _opP = sp.sympify(_d['opP']); _detf = sp.sympify(_d['detf'])
_numVP = sp.cancel(_vpP*_detf); _numOmP = sp.cancel(_opP*_detf)
_argsH = (_N0,_A0,_o0,_V0, _Nx,_Ax,_ox,_Vx, _k, _Np,_Ap,_opp,_vpp)
_f_numVP  = sp.lambdify(_argsH, _numVP, 'numpy')
_f_numOmP = sp.lambdify(_argsH, _numOmP, 'numpy')

def Lmat(N,A,om,V, Nb,Ab,ob,Vb, k):
    return np.array(_fL(N,A,om,V, Nb,Ab,ob,Vb, k), dtype=complex)

def numVP(N,A,om,V, Nb,Ab,ob,Vb, k, hp):
    return complex(_f_numVP(N,A,om,V, Nb,Ab,ob,Vb, k, hp[0],hp[1],hp[2],hp[3]))
def numOmP(N,A,om,V, Nb,Ab,ob,Vb, k, hp):
    return complex(_f_numOmP(N,A,om,V, Nb,Ab,ob,Vb, k, hp[0],hp[1],hp[2],hp[3]))

if __name__ == "__main__":
    import css_core as C
    S3=np.sqrt(3.0)
    # background slopes at sonic point (exact): N'=-2/S3, A'=3, om'=9/2, V'=2/S3 (branch B1)
    Nb,Ab,ob,Vb = -2/S3, 3.0, 4.5, 2/S3
    Npt = (2/S3, 1.5, 0.75, 1/S3)
    print("L at the sonic point is singular (detf->0); check indicial at CENTER instead.")
    # center: A=1,om=0,V=0,N large. background slopes there: N'=-N, A'->0, om'->0, V'->V(=0)
    # Use a near-center point
    nc=1.0; z=1e-6; N=nc/z; A=1.0; om=1.5*0.25*z*z; V=z/(2*nc)
    Nb=C.Nx(A,N,om,V); Ab=C.Ax(A,N,om,V); ob=C.omx(A,N,om,V); Vb=C.Vx(A,N,om,V)
    for k in [2.81, 0.357, 0.0]:
        L = Lmat(N,A,om,V, Nb,Ab,ob,Vb, k)
        ev = np.linalg.eigvals(L)
        print(f"  k={k}: eig(L) at x=ln(z)={np.log(z):.2f} -> {np.round(ev,4)}")
