# check_branch.py — which separatrix branch does the integrated background follow at the sonic point?
# Compare the measured (om',V') approaching sonic against branch1 (om'=9/2,V'=2/sqrt3) and
# branch2 (om'=27/2,V'=4/sqrt3). Also check the center expansion self-consistency.
import numpy as np, math
import css_core as C
S3=math.sqrt(3.0)
A2STAR=0.250067905275

def bg_slopes(Y):
    N,A,om,V=Y; return (C.Nx(A,N,om,V),C.Ax(A,N,om,V),C.omx(A,N,om,V),C.Vx(A,N,om,V))

# integrate to near sonic, sample slopes just before
Y=C.center_seed(1.0,A2STAR,math.exp(-12)); x=math.log(math.exp(-12)); h=5e-5
prevY=Y
samples=[]
for i in range(int((1.0-x)/h)):
    Y2=C.rk4_step(Y,h)
    D=C.Dson(Y2[0],Y2[3])
    if D<0 and -D<0.02:
        sl=bg_slopes(Y2)
        samples.append((-D, Y2, sl))
    if D>=0: break
    Y=Y2; x+=h
print("Approaching sonic: (om',V') vs branch1(4.5,1.1547) branch2(13.5,2.3094)")
print(f"{'|D|':>10} {'N':>9} {'V':>9} {'om_slope':>10} {'V_slope':>10}")
for d,Y,sl in samples[::max(1,len(samples)//12)]:
    print(f"{d:10.2e} {Y[0]:9.5f} {Y[3]:9.5f} {sl[2]:10.4f} {sl[3]:10.4f}")
print(f"\nbranch1: om'=4.5000 V'=1.1547")
print(f"branch2: om'=13.5000 V'=2.3094")
# center check: does the seed satisfy the ODE? compare seed-implied slope to RHS slope
print("\nCenter seed consistency (ratio seed-slope / RHS-slope, should be ~1):")
for lz in [-12,-14,-16]:
    z=math.exp(lz); Y=C.center_seed(1.0,A2STAR,z)
    # analytic seed slopes d/dx: N=nc/z=nc e^-x -> dN/dx=-N; A=1+a2 z^2 -> dA/dx=2 a2 z^2; om=1.5a2z^2->dom/dx=3a2z^2; V=z/(2nc)->dV/dx=z/(2nc)=V
    seed_sl=(-Y[0], 2*A2STAR*z*z, 3*A2STAR*z*z, Y[3])
    rhs_sl=bg_slopes(Y)
    print(f"  z=e^{lz}: N {seed_sl[0]/rhs_sl[0]:.4f}  A {seed_sl[1]/rhs_sl[1]:.4f}  om {seed_sl[2]/rhs_sl[2]:.4f}  V {seed_sl[3]/rhs_sl[3]:.4f}")
