# Check h-convergence of the EC shoot near oi*~0.375. Print lhop(first crossing) vs oi for several h.
import numpy as np, math
import hka_shoot_oi as S

if __name__=="__main__":
    N_inf=1.0
    ois=np.arange(0.370,0.381,0.001)
    print(f"{'oi':>8}", end="")
    hs=[1e-4,5e-5,2.5e-5,1e-5]
    for h in hs: print(f"  lhop(h={h:.0e})", end="")
    print("   |  V@cross(h=1e-5)")
    for oi in ois:
        row=f"{oi:8.4f}"
        Vc=None
        for h in hs:
            val,r=S.F(round(oi,6),N_inf,h=h)
            row+=f"  {val:12.3e}" if val is not None else f"  {'break':>12}"
            if h==1e-5 and r.get('status')=='cross': Vc=r['V']
        row+=f"   |  {Vc if Vc is not None else float('nan'):.6f}"
        print(row)
