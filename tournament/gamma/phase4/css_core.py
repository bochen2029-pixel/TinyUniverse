# css_core.py — fast scalar numeric core for the radiation-fluid CSS system.
# Loads css_syms.pkl once, lambdifies with 'math' (scalar), exposes:
#   slopes(Y) -> (N',A',om',V')   [the resolved autonomous ODE, correct full-4D system]
#   Dson(N,V), numV,numOm,det (r=1 normalized), and the metric slopes.
import sympy as sp, pickle, math

_d = pickle.load(open("css_syms.pkl", "rb"))
_r = sp.Symbol('r', positive=True)
_S = {k: sp.sympify(v).subs(_r, 1) for k, v in _d.items()}
_A0, _N0, _o0, _V0 = sp.symbols('A0 N0 om0 V0')
_sy = (_A0, _N0, _o0, _V0)

# scalar (math) lambdifications
_omx  = sp.lambdify(_sy, _S['omx'], 'math')
_Vxx  = sp.lambdify(_sy, _S['Vxx'], 'math')
_Ax   = sp.lambdify(_sy, _S['Axv'], 'math')
_Nx   = sp.lambdify(_sy, _S['Nxv'], 'math')
_numV = sp.lambdify(_sy, _S['numV'], 'math')
_numO = sp.lambdify(_sy, _S['numOm'], 'math')
_det  = sp.lambdify(_sy, _S['det'], 'math')

def Nx(A, N, om, V): return _Nx(A, N, om, V)
def Ax(A, N, om, V): return _Ax(A, N, om, V)
def omx(A, N, om, V): return _omx(A, N, om, V)
def Vx(A, N, om, V): return _Vxx(A, N, om, V)
def numV(A, N, om, V): return _numV(A, N, om, V)
def numOm(A, N, om, V): return _numO(A, N, om, V)
def det(A, N, om, V): return _det(A, N, om, V)

def Dson(N, V): return 3*N*N*V*V - N*N - 4*N*V - V*V + 3   # self-consistent singular locus

def slopes(Y):
    """Y=(N,A,om,V) -> dY/dx via the resolved autonomous system."""
    N, A, om, V = Y
    return (Nx(A, N, om, V), Ax(A, N, om, V), omx(A, N, om, V), Vx(A, N, om, V))

def rk4_step(Y, h):
    k1 = slopes(Y)
    k2 = slopes(tuple(Y[i] + 0.5*h*k1[i] for i in range(4)))
    k3 = slopes(tuple(Y[i] + 0.5*h*k2[i] for i in range(4)))
    k4 = slopes(tuple(Y[i] + h*k3[i] for i in range(4)))
    return tuple(Y[i] + (h/6.0)*(k1[i] + 2*k2[i] + 2*k3[i] + k4[i]) for i in range(4))

# center seed (z=e^x -> 0). DERIVED (center_series.py) and numerically verified (Vx/V->+1):
#   N = nc/z,   A = 1 + a2 z^2,   om = (3/2) a2 z^2,   V = +z/(2 nc)   (OUTGOING, v1=1/(2nc))
#   n1 = 0.  Free params: nc (x-translation gauge) and a2 (central density, the shooting parameter).
def center_seed(nc, a2, z0):
    """Return Y=(N,A,om,V) at x0=ln(z0), from the verified regular-center series."""
    N = nc/z0
    A = 1.0 + a2*z0*z0
    om = 1.5*a2*z0*z0
    V = z0/(2.0*nc)
    return (N, A, om, V)
