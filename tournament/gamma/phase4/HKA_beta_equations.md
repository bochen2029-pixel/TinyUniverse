# HKA radiation-fluid CSS equations for β = 0.35580192 — primary-source retrieval

> **Provenance:** background research subagent (`abde345b`), 2026-07-12, TINY UNIVERSE autonomous overnight run.
> Retrieved & read in full: Evans–Coleman (gr-qc/9402041), Koike–Hara–Adachi PRL 1995 (gr-qc/9503007),
> Hara–Koike–Adachi long paper 1997/1999 (gr-qc/9607010 = PRD 59, 104008), Maison 1996 (gr-qc/9504008),
> Gundlach–Martín-García review (arXiv:0711.4620). Ori–Piran PRD 42, 1068 (1990) NOT retrieved (pre-arXiv,
> paywalled) — but **not needed**: its content (sonic-point analyticity classification) is reproduced in
> Maison Eq. 8 and HKA §IV C.
>
> **THE BUG DIAGNOSIS (why Stage-B gave β≈0.997):** HKA gr-qc/9607010 footnote 15 states the spurious
> GAUGE mode sits at **κ≈0.35699 in the sonic-point gauge N̄_p(0)=0** but **at κ=1 in the origin gauge
> N̄_p(−∞)=0** (Maison's gauge). Our Stage-B κ₀≈1.003→β≈0.997 IS this κ=1 origin-gauge gauge mode. The
> PHYSICAL relevant mode is **κ≈2.81055255 → β=1/κ=0.35580192** in EITHER gauge. Fix = (a) HKA Eq. 4.1
> background with a genuine regular center, (b) integrate perturbation Eq. 5.13 from the sonic point x=0
> outward in gauge N̄_p(0)=0, (c) shoot on κ to kill the growing mode at x→−∞, DISCARD the gauge mode.
>
> **Convention:** 0.3558 = the radiation-FLUID (Γ=4/3) mass exponent (β≡γ≡γ_BH=1/Re λ0=1/κ). The ≈0.374
> value is the massless SCALAR field — a DIFFERENT matter model, not a different convention. Do not mix them.
> Cleanest ready-to-code source: **HKA gr-qc/9607010** (one self-consistent gauge, all pieces present).

---

## (1) Variables, ansatz, gauge (HKA gr-qc/9607010 §III)

Metric (polar-radial), Eq. (3.1):  `ds² = −α²(t,r)dt² + a²(t,r)dr² + r²dΩ²`
4-velocity (3.2): `u^t=−α/√(1−V²), u^r=aV/√(1−V²)`; V = 3-velocity. EOS (3.4): **p=(γ−1)ρ**, radiation γ=4/3 ⇒ p=ρ/3.
Self-similar coords + rescaled fields (3.5):
```
s ≡ −ln(−t),   x ≡ ln(−r/t)
N ≡ α/(a·eˣ),   A ≡ a²,   ω ≡ 4π r² a² ρ
```
Phase-space = (A, N, ω, V), functions of x only for self-similar solutions.
**Gauge:** residual freedom t↦F⁻¹(t) ⇒ (s,x)↦(f⁻¹(s), x−s+f⁻¹(s)); A,ω,V invariant, N↦ḟN. **Fix N at ONE point per s-line.
HKA fix the SONIC POINT at x=0.** This sets β_scaling=1, so β_BH=1/κ with no extra factor.
Sound speed c_s=√(γ−1)=1/√3≈0.57735.

## (2) Background ODEs (self-similar; HKA Eq. 4.1)
```
A_x/A = 1 − A + 2ω[1+(γ−1)V²]/(1−V²)                                    (4.1a)
N_x/N = −2 + A − (2−γ)ω                                                 (4.1b)
(A_x/A)_x = −2γ N V ω/(1−V²)                                            (4.1c)
(1+NV)ω_x/ω + γ(N+V)V_x/(1−V²)
      = (3(2−γ)/2)NV − ((2+γ)/2)ANV + (2−γ)NVω                          (4.1d)
(γ−1)(N+V)ω_x/ω + γ(1+NV)V_x/(1−V²)
      = (2−γ)(γ−1)Nω + ((7γ−6)/2)N + ((2−3γ)/2)AN                       (4.1e)
```
Constraint (4.2, from 4.1a&4.1c):  `1 − A + 2ω[1+(γ−1)V²]/(1−V²) = −2γNVω/(1−V²)`.
Use {4.1a,4.1b,4.1d,4.1e} as the basic set; 4.1c/4.2 as checks.
Maison independent cross-check (gr-qc/9504008 Eq. 8):  `A² = 1 − 2ρ̃(1+(1+k)v/B+kv²)/(1−v²)`.

## (2C) Regular center, x→−∞ (HKA §IV C 1, Eq. 4.10–4.13)
M ≡ NV. Fixed point of (4.1):  **A=1, M=−2/(3γ), ω=V=0**  (4.11) — exactly ONE unstable direction (along M+2/(3γ)).
Asymptotics (4.12): `A∼1+A₋∞ e^{2x}, N∼N₋∞ e^{−x}, ω∼ω₋∞ e^{2x}, V∼V₋∞ eˣ`;
(4.13): **A₋∞=(2/3)ω₋∞,  N₋∞V₋∞=−2/(3γ)**. Regularity = analytic ∀x AND A=1,V=0 at x=−∞.
> A system that "lacks a regular center" will not have fixed point (4.11). Most likely typo culprit = a sign in
> the ω or A term of (4.1a/4.1b) or in constraint (4.2). Compare term-by-term.

## (2D) Sonic point — L'Hôpital / coordinate-singularity removal (HKA Eq. 4.4–4.9)
Write (4.1d,e) as `[[a,b],[c,d]]·[ω_x,V_x]ᵀ=[e,f]ᵀ` (4.4). Sonic point = det vanishes (4.5):
```
det ∝ γ{(1+NV)² − (γ−1)(N+V)²}/(1−V²) = 0        (flow speed = sound speed √(γ−1))
```
Finite derivatives ⇒ row-proportionality (L'Hôpital) (4.6): `a·f − e·c = 0`.
(4.5)+(4.6)+(4.2) fix sonic data in terms of single free param V₀≡V(0) (4.7–4.9):
```
N₀ = (1 − √(γ−1)V₀)/(√(γ−1) − V₀)                                        (4.7)
A₀ = [γ²+4γ−4 + 8(γ−1)^{3/2}V₀ − (3γ−2)(2−γ)V₀²]/[γ²(1−V₀²)]             (4.8)
ω₀ = 2√(γ−1)(√(γ−1)−V₀)(1+√(γ−1)V₀)/[γ²(1−V₀²)]                          (4.9)
```
Regular sonic solution = 1-param family in V₀. For γ=4/3: √(γ−1)=1/√3; A₀ encodes Misner–Sharp 2m/r=1−1/A₀=1/3.
Higher expansion coeffs: expand (4.1) about x=0 (first-order coeffs solve a quadratic; pick real-derivative branch).

## (3) Boundary conditions (two-sided shoot)
- **(a) Center x→−∞:** A=1, V=0; decays per (4.12)–(4.13). Only discrete V₀ reach it (one unstable dir, 4.11).
- **(b) Sonic x=0:** det=0 (4.5) AND a·f−e·c=0 (4.6); data A₀,N₀,ω₀ (4.7–4.9) param by V₀; gauge = sonic pt at x=0.
- **(c) Outer x→+∞:** stationary pt N=0, A=1+2ω[1+(γ−1)V²]/(1−V²) (4.14), attracting node if ω≤(1−V²)/[γ(1+V²)] (4.15).
  Physical (Evans–Coleman) solution = the one with EXACTLY ONE zero of V (first in the Bogoyavlenskii sequence). Not asymptotically flat (a→a∞≈1.07).
- **Shoot (HKA §IV E):** fix V₀; start from sonic series (4.7–4.9); integrate (4.1) x=0→−∞ (RK4). Reject V₀ if A<1
  or if a second sonic point appears. EC solution = bracketed root between the A<1 case and the second-sonic case.

## (4) Perturbation eigenvalue problem (HKA §V)
Perturb logs Ā=lnA, N̄=lnN, ω̄=lnω: `h(s,x)=H_ss(x)+ε h_p(x)e^{κs}`, s=−ln(−t) (growth as t→0⁻ ⇔ Re κ>0).
Linear system Ψ=(Ā_p,N̄_p,ω̄_p,V_p)ᵀ, HKA Eq. (5.13) with coefficients (5.5)–(5.10):
```
As=1, Bs=γV/(1−V²), Cs=(γ−1)V, Ds=γ/(1−V²)                               (5.5)
Ax=1+NV, Bx=γ(N+V)/(1−V²), Cx=(γ−1)(N+V), Dx=γ(1+NV)/(1−V²)             (5.6)
E1=−((γ+2)/2)ANV
E2=((6−3γ)/2)NV−((2+γ)/2)ANV+(2−γ)ωNV−NV ω̄_x−γN V_x/(1−V²)
E3=(2−γ)ωNV
E4=((6−3γ)/2)N−((2+γ)/2)AN+(2−γ)ωN−N ω̄_x−γ(1+2NV+V²)V_x/(1−V²)²          (5.7)
F1=((2−3γ)/2)AN
F2=(2−γ)(γ−1)ωN+((7γ−6)/2)N+((2−3γ)/2)AN−(γ−1)N ω̄_x−γNV V_x/(1−V²)
F3=(2−γ)(γ−1)ωN
F4=−(γ−1)ω̄_x−γ(N+2V+NV²)V_x/(1−V²)²                                     (5.8)
G1=−A, G3=2{1+(γ−1)V²}ω/(1−V²), G4=4γωV/(1−V²)²                          (5.9)
H1=A, H3=(γ−2)ω                                                          (5.10)
```
Matrix form (5.13):  `diag(1,1,Ax,Cx… )` block — d/dx of Ψ with the (Ā,N̄) rows first-order in (G,H) and the
(ω̄,V) rows given by [[Ax,Bx],[Cx,Dx]]·(ω̄_p',V_p') = [E·Ψ, F·Ψ] + κ[[As,Bs],[Cs,Ds]]·(ω̄_p,V_p). (Assemble from 5.5–5.10.)
Auxiliary (Ā) identity (5.14, check): `κĀ_p+∂_xĀ_p = −(2γNVω/(1−V²))(N̄_p+ω̄_p) − (2γNω(1+V²)/(1−V²)²)V_p`.

**BCs:** (i) analytic ∀x; (ii) regular center Ā_p=0 at x=−∞. Sonic x=0 is a REGULAR SINGULAR point
(AxDx−BxCx=0 = same as 4.5); row-proportionality + identity (5.15) + **gauge N̄_p(0)=0** fix the analytic
perturbation up to scale by the single complex κ. x→−∞ modes: {1,1,e^{−2x}, e^{−x}·(1,−1,(2−3γ)/(2(γ−1)),0)};
e^{−x} excluded by (5.15); **kill the growing e^{−2x} mode by choosing κ** (the eigenvalue condition). x→+∞
modes {1,e^{−κx},e^{−κx}} bounded iff Re κ>0.

**GAUGE MODE (5.20) — DISCARD:** `h̃_gauge=e^{κ̄s}·{h'_ss(x) for Ā,ω̄,V; h'_ss(x)+κ̄ for N̄}`, κ̄ arbitrary.
Sonic gauge N̄_p(0)=0 ⇒ κ̄=−dN̄_ss/dx|₀≈0.35699. **Origin gauge N̄_p(−∞)=0 ⇒ κ̄=1.** (← our β≈0.997 artifact.)

**Result (HKA §V G 1, Table I; KHA-1995 §5):** unique physical relevant mode
```
κ = 2.8105525488   ⇒   β_BH = 1/κ = 0.35580192          (verified: 1/2.81055255 = 0.35580192, 8 s.f.)
```
KHA-1995 Eq. (14–15): `M_BH ~ O(ε^{1/Re κ})  ⇒  β = 1/Re κ`. Maison (gr-qc/9504008 Eq. 9–11) independent: growth
e^{−ωτ}, τ=ln(−t); γ=1/ω, same 0.3558. Gundlach (0711.4620 Eq. 1,6–10): M∝(p−p*)^γ, γ=1/λ0, λ0=κ=ω≈2.811.

## (5) Numbers table
| source | printed | value (γ=4/3) | cite |
|---|---|---|---|
| Evans–Coleman 1994 | β (mass exp) | β≈0.36 (η_c=1.018828234) | gr-qc/9402041 |
| KHA 1995 PRL | κ; β=1/Reκ | **κ≈2.81055255 ⇒ β≈0.35580192** | gr-qc/9503007 |
| HKA 1997/1999 | κ; β_BH=1/κ | **κ=2.8105525488 ⇒ 0.35580192** | gr-qc/9607010 Tbl I/II |
| Maison 1996 | γ=1/ω | 0.3558 | gr-qc/9504008 Tbl 1 |
| Gundlach review | γ | 0.3558 (μ=(5/2)γ≈0.8895) | 0711.4620 §4.2.1 |
| EC similarity exp n | n≈1.1485 | (NOT β — background homothety) | gr-qc/9402041 |

Curvature scaling (Gundlach Eq. 31): `R_max ∝ (p−p*)^{2γ}`, 2γ≈0.7116. (ρ_max version is DERIVED via EC's
R_μνR^μν=(256/3)π²ρ², not literally quoted — treat as inference.)

## (6) Gaps
- Ori–Piran PRD 42,1068 (1990): not retrieved (pre-arXiv/paywalled); NOT needed (content in Maison Eq. 8, HKA §IV).
- HKA higher-order sonic coeffs: omitted in the paper ("we do not give explicit expressions") — generate by
  symbolic expansion of (4.1)/(5.13) about x=0 yourself.
- Bogoyavlenskii JETP 46,633 (1977): CSS-sequence source, not retrieved, not needed.
