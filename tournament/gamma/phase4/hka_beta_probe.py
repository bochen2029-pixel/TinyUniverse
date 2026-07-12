# hka_beta_probe.py — targeted eigenvalue scan on the CORRECTED operator (hka_pert_derive L).
# Delta(kappa)=det[v_center | s1 s2 s3] (columns unit-normalized) + rejection ||c-QQ^H c||.
# Zero <=> center-regular mode lies in the sonic-analytic span <=> eigenvalue.
# Physical target kappa~2.81055255 -> beta=1/kappa=0.35580192. Gauge modes (0.357,1.0) discarded (fn15).
import numpy as np, hka_beta4 as B
B.bg(); B.bg_path()
print(f"sonic x_s={B.bg()['xs']:.5f}")
print(f"{'kappa':>10} {'|Delta|':>11} {'Re':>12} {'Im':>12} {'rej':>10}")
grid=[0.357,0.5,1.0,1.5,2.0,2.4,2.6,2.7,2.75,2.80,2.81055255,2.82,2.85,2.9,3.0,3.5,4.0,5.0]
for kap in grid:
    d=B.Delta(kap)
    if d is None: print(f"{kap:10.5f}   fail"); continue
    det,rej=d
    print(f"{kap:10.5f} {abs(det):11.3e} {det.real:+12.3e} {det.imag:+12.3e} {rej:10.2e}")
