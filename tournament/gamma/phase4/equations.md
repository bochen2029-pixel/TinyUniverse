# Self-Similar Gravitational Critical Collapse — Verified Equation Reconstruction

**Purpose.** Ground a single-file fp64 CPU ODE-eigenvalue solver for the critical exponents of
gravitational collapse. Every load-bearing equation and number below is transcribed from a primary
source (fetched, PDF-to-text converted, and read), not reconstructed from memory. Anything not
directly verifiable against a fetched source is tagged `[UNVERIFIED]`.

**Sources fetched and read (full text via PDF→text conversion):**

- **[KHA95]** T. Koike, T. Hara, S. Adachi, *Critical behaviour in gravitational collapse of radiation
  fluid — A renormalization group (linear perturbation) analysis*, **Phys. Rev. Lett. 74 (1995) 5170**,
  arXiv:**gr-qc/9503007**. — the CSS fluid oracle, read verbatim.
- **[G97]** C. Gundlach, *Understanding critical collapse of a scalar field*, **Phys. Rev. D 55 (1997) 695**,
  arXiv:**gr-qc/9604019**. — the DSS scalar construction, read verbatim.
- **[EC94]** C. R. Evans, J. S. Coleman, *Observation of critical phenomena and self-similarity in the
  gravitational collapse of radiation fluid*, **Phys. Rev. Lett. 72 (1994) 1782**, arXiv:**gr-qc/9402041**.
  — β≈0.36 from evolution (abstract/search-verified; full PDF not separately parsed).
- **[LRR]** C. Gundlach & J. M. Martín-García, *Critical Phenomena in Gravitational Collapse*,
  **Living Reviews in Relativity 10 (2007) 5**. Consolidated cross-check read via the freely-available
  arXiv preprint **gr-qc/0001046** (the 2000 version; identical on all points cited here). The published
  2007 update is arXiv:0711.4620.

> **Convention note (important for the coder).** The three primary sources use *different* coordinates
> and variable names. This doc keeps each source's own notation within its section (so equations are
> transcribed exactly and can be checked against the paper), and flags where they connect. Do **not**
> mix the KHA `(s,x,N,A,ω,V)` set with the Gundlach `(τ,ζ,a,g,X±)` set.

---

# PRIORITY 1 — CSS radiation perfect fluid, `p = ρ/3` (Evans–Coleman / Koike–Hara–Adachi)

This is the clean case: the homothetic (continuously self-similar, CSS) ansatz reduces the Einstein +
perfect-fluid equations to an **autonomous ODE system in one variable**. The critical solution is the
member analytic through the **sonic point** and regular at the center — a two-point boundary-value /
shooting eigenvalue problem. Its single unstable linear mode gives β. **This validates the whole method
on a known number: β ≈ 0.3558.** All equations in this section are transcribed from **[KHA95]**.

## 1.1 Similarity variables and metric ansatz  `[VERIFIED — KHA95 eqs (16),(18) and text §3]`

Equation of state: perfect fluid `p = (γ_ad − 1) ρ` with adiabatic index `γ_ad = 4/3`, i.e. **`p = ρ/3`**
(radiation fluid). Stress tensor `T_ab = ρ u_a u_b + p(g_ab + u_a u_b)`.

General spherically symmetric line element (**KHA eq. 16**):

$$
ds^2 = -\,\alpha^2(t,r)\,dt^2 + a^2(t,r)\,dr^2 + r^2\left(d\theta^2 + \sin^2\theta\,d\phi^2\right).
$$

The **only** residual coordinate freedom preserving this form is the time relabeling `t ↦ F^{-1}(t)`
(**KHA eq. 17**) — this is the gauge freedom you fix at the sonic point (see §1.3).

Fluid 4-velocity components, with `V` the fluid 3-velocity:

$$
u^t = \frac{a}{\sqrt{1-V^2}}\quad\text{(KHA writes }u_t=a(1-V^2)^{-1/2}\text{)},\qquad
u^r = \frac{aV}{\sqrt{1-V^2}}.
$$

**Similarity variables (KHA §2, eq. text):**

$$
\boxed{\,s \equiv -\ln(-t),\qquad x \equiv \ln\!\left(-\frac{r}{t}\right)\,}\qquad (t<0)
$$

`s` is the "scaling variable" (RG time; translation in `s` = the scaling symmetry), `x` the "spatial"
self-similar coordinate. Inverse relation used later: `r = e^{x-s}`.

**Reduced (scale-invariant) matter/metric variables (KHA §3):**

$$
N \equiv \alpha\,a^{-1}\,e^{-x},\qquad A \equiv a^2,\qquad \omega \equiv 4\pi r^2 a^2 \rho .
$$

These four `(N, A, ω, V)` are the dependent variables that make the system **autonomous**.

> **Relation to the Evans–Coleman ξ.** [EC94] and [LRR §3.2] use a self-similar coordinate
> `x_EC = r/f(t)` with `f(t)=t^n` (EC call it "self-similarity of the second kind", but the spacetime
> is genuinely homothetic / "first kind" — the difference is only that EC's `t` is proper time at
> infinity, not at the origin). KHA's `x = ln(−r/t)` is the log-form of the same self-similar slicing.
> Either is fine to code; KHA's is the one with the explicit autonomous ODEs below, so **use KHA**.

## 1.2 The autonomous ODE system  `[VERIFIED — KHA95 eq. (18), transcribed term-by-term]`

Full system of evolution equations (before imposing self-similarity). `_,x` and `_,s` are partial
derivatives. Transcribed **exactly** from KHA eq. (18):

**(a) Constraint / metric slope equations (already ODEs in `x`):**

$$
\frac{A_{,x}}{A} = 1 - A + \frac{2\omega}{1-V^2}\left(1 + \frac{V^2}{3}\right),
$$

$$
\frac{N_{,x}}{N} = -2 + A - \frac{2\omega}{3}.
$$

**(b) Fluid equations (the two dynamical PDEs):**

$$
\frac{\omega_{,s} + (1+NV)\,\omega_{,x}}{\omega}
+ \frac{4\{\,V V_{,s} + (N+V)V_{,x}\,\}}{3(1-V^2)}
- \frac{NV A_{,x}}{3A}
+ \frac{4V N_{,x}}{3}
+ 2N\left(1 + \frac{4\omega}{9(1-V^2)}\right) = 0,
$$

$$
\frac{4V\,\omega_{,s} + (4V + N + 3NV^2)\,\omega_{,x}}{\omega}
+ \frac{4\{(1+V^2)V_{,s} + (1+V^2+2NV)V_{,x}\}}{1-V^2}
+ \frac{N(1-V^2)A_{,x}}{A}
+ 4(1+V^2)N_{,x} + 2N(1+3V^2) = 0.
$$

> **Transcription caveat (coder must verify against the PDF).** The two fluid equations in KHA eq. (18)
> are typeset across several lines; in the second equation the outermost `{…}` and the overall bracket
> balance should be re-checked directly against **gr-qc/9503007 p.3** before trusting signs, because PDF
> line-wrapping makes one closing brace ambiguous. The *structure* (two first-order combinations of
> `ω_,s, ω_,x, V_,s, V_,x` plus algebraic source terms) is certain. This is flagged as the single most
> important line-level check. `[VERIFY-BRACES]`

**Self-similar reduction.** Assume `N = N_ss(x)`, `A = A_ss(x)` depend on `x` only. KHA state (and prove
in their ref. [8]) that any spherically symmetric self-similar spacetime can be put in this form using
the gauge freedom (17). It then **follows** that `ω = ω_ss(x)` and `V = V_ss(x)` are also functions of
`x` alone. Setting all `_,s → 0` collapses eq. (18) to a closed **autonomous ODE system in `x`** for
`(N_ss, A_ss, ω_ss, V_ss)`. The first two equations (a) are already in `dX/dx` form; the two fluid
equations (b), with `_,s = 0`, become two linear algebraic relations for `ω_,x` and `V_,x` that you
solve (2×2) to get `ω_ss'(x)` and `V_ss'(x)`. **This is the ODE system the solver integrates.**

## 1.3 Boundary / regularity conditions — the sonic-point eigenvalue condition  `[VERIFIED — KHA95 §4]`

Three conditions; the sonic-point one is the shooting/eigenvalue condition.

1. **Analyticity everywhere:** the self-similar solution must be analytic for all `x ∈ ℝ` (KHA cond. i).

2. **Regular center** `x → −∞` (KHA cond. ii): the spacetime and matter are regular there, i.e.

$$
\boxed{A = 1,\qquad V = 0 \quad\text{at the center } x=-\infty.}
$$

3. **Sonic point (the selection condition).** Gauge-fix by placing the sonic point at **`x = 0`**. The
   sonic point is where *"the velocity of the fluid particle seen from the observer on the constant-`x`
   line is equal to the speed of sound `1/√3`."* At this point the ODE system has a **singular point**
   (a coefficient — the `dω/dx, dV/dx` denominator from the 2×2 inversion — vanishes). Following
   **Ori & Piran, Phys. Rev. D 42 (1990) 1068**, analyticity requires the solution be analytic *through*
   this singular point. Power-series expansion about `x=0` shows the locally-analytic solutions form a
   **one-parameter family**, parameterized by (say) `V_ss(0)`. Imposing the regular-center condition (2)
   then restricts `V_ss(0)` to a **discrete set** — this is the eigenvalue/shooting condition that
   **selects the Evans–Coleman critical solution `H_ss`.**

   > **Coder's recipe.** (i) At `x=0` require the sound-speed condition to hold — for a fluid element on
   > a constant-`x` worldline the local flow speed equals `c_s = 1/√3`; combined with the vanishing of
   > the 2×2 denominator this fixes the sonic-point locus. (ii) Do a Frobenius/Taylor expansion at `x=0`
   > to get regular initial data (the L'Hôpital limit for `ω_,x, V_,x` at the removable singularity),
   > with `V_ss(0)` as the shooting parameter. (iii) Integrate outward (`x→+∞`) and inward toward the
   > center (`x→−∞`); adjust `V_ss(0)` until `A→1, V→0` at the center. `[VERIFIED-STRUCTURE; the explicit
   > sonic-point algebraic equation V_ss(0) satisfies is NOT written in closed form in KHA — you derive
   > it from the denominator=0 condition. This is the main derivation the builder must do — see §1.6.]`

## 1.4 Linear perturbations and the eigenvalue problem  `[VERIFIED — KHA95 §5, eqs (19),(7)]`

Linearize each field `h ∈ {N, A, ω, V}` about the critical solution `H_ss(x)` (**KHA eq. 19**):

$$
h(s,x) = H_{ss}(x) + \epsilon\, h_{\text{var}}(s,x).
$$

Fix the residual time-gauge (17) by requiring `N_var(s, 0) = 0`. Seek eigenmodes with **exponential
`s`-dependence** (this is the key ansatz — `s` is RG time, so growth is exponential, not periodic, in
the CSS case):

$$
\boxed{\,h_{\text{var}}(s,x) = h_p(x)\,e^{\kappa s},\qquad \kappa \in \mathbb{C}.\,}
$$

Substituting into eq. (18) and linearizing yields **linear, homogeneous, first-order ODEs** for the
mode profiles `(N_p, A_p, ω_p, V_p)`. Boundary/analyticity conditions mirror the background:

- perturbations analytic for all `x` (the sonic point is now a **regular singular point** of the
  perturbation ODEs);
- regular center: `A_p, V_p` finite at `x → −∞`.

Analyticity through the sonic point again leaves a one-parameter family labeled by `κ`; regularity at
the center then allows only **discrete `κ`**. This is the eigenvalue problem. The RG argument (KHA §2,
eqs 8–15) gives the master relation

$$
\boxed{\;\beta = \frac{1}{\operatorname{Re}\kappa_0}\;}\qquad\text{(KHA eq. 15),}
$$

where `κ_0` is the single **relevant** (`Re κ > 0`) mode. KHA verify numerically there is exactly **one**
relevant mode (checked over `0 ≤ Re κ ≤ 15, |Im κ| ≤ 14`).

> **Gauge-mode warning (do not mistake it for the physical mode).** KHA report an **unphysical gauge
> mode at `κ ≈ 0.35699`** in their gauge, arising from a coordinate transformation of the self-similar
> solution. It is numerically close-ish to nothing important but must be discarded; the physical relevant
> eigenvalue is `κ ≈ 2.81` (below). `[VERIFIED — KHA95 §5, last paragraph]`

## 1.5 Verified numbers  `[VERIFIED]`

| quantity | value | source |
|---|---|---|
| relevant eigenvalue `κ_0` (`= Re κ_0`, mode is real) | **`κ ≈ 2.81055255`** | KHA95 §5 (verbatim) |
| critical exponent `β = 1/Re κ_0` | **`β ≈ 0.35580192`** (abstract rounds to `0.3558019`) | KHA95 abstract + §5 |
| `β` from independent *evolution* (Evans–Coleman) | **`β ≈ 0.36`** | EC94 abstract; LRR §2.2 |
| spurious gauge mode (discard) | `κ ≈ 0.35699` | KHA95 §5 |

So the user's recollection is confirmed: **β ≈ 0.3558** (specifically `0.35580192`), from
**Re κ₀ ≈ 2.81055255**; and Evans–Coleman's evolution value is **β ≈ 0.36**. Cross-checked in [LRR],
which additionally credits Maison (Phys. Rev. D — via Comm. Math. Phys.) and KHA with the same
perturbative computation. `[Maison's own digit not independently re-verified here — UNVERIFIED to the
last decimal, but LRR states agreement.]`

## 1.6 The one gap a builder must close for Priority 1

**The explicit closed-form sonic-point regularity equation is not printed in KHA.** KHA give the *rule*
("analytic through the sonic point where flow speed = `1/√3`", singular point of the ODE, one free
parameter `V_ss(0)` fixed to discrete values by the center BC) but not the algebraic equation
`V_ss(0)` must satisfy, nor the L'Hôpital expressions for `(ω_ss', V_ss')` at `x=0`. The builder must
**derive these from the autonomous system in §1.2** by: (1) writing the 2×2 solve for `(ω_,x, V_,x)`
explicitly, (2) setting its determinant → 0 to locate the sonic point, (3) demanding the numerators
vanish simultaneously (removable singularity) to get the regular local expansion. **The canonical fully
worked-out version of exactly this is Ori & Piran, Phys. Rev. D 42 (1990) 1068–1090** (KHA's ref. [9]) —
fetch that paper if the derivation is error-prone. `[ACTION for builder — primary derivation task.]`

---

# PRIORITY 2 — DSS massless scalar field (Choptuik / Gundlach)

The ambitious target. The critical solution here is **discretely** self-similar (DSS): periodic in the
log-scale coordinate with period `Δ`. This turns the eigenvalue problem from a 1-D ODE (CSS) into a
**periodic hyperbolic boundary-value problem** (a nonlinear eigenvalue problem for `Δ` and a periodic
function `ξ_0(τ)`) — genuinely harder. All equations transcribed from **[G97]** (gr-qc/9604019).

## 2.1 Field equations (first-order form)  `[VERIFIED — G97 eqs (1)-(9)]`

Einstein + massless real scalar `φ`, minimally coupled:

$$
G_{ab} = 8\pi G\left(\phi_{,a}\phi_{,b} - \tfrac12 g_{ab}\phi_{,c}\phi^{,c}\right),\qquad \phi^{;c}{}_{,c}=0.
$$

Metric (same polar-radial form as the fluid case, **G97 eq. 2**):

$$
ds^2 = -\alpha(r,t)^2\,dt^2 + a(r,t)^2\,dr^2 + r^2(d\theta^2+\sin^2\theta\,d\phi^2).
$$

Auxiliary first-order matter fields (**G97 eq. 3**):

$$
X(r,t)\equiv \frac{\sqrt{2\pi G}\,r}{a}\,\phi_{,r},\qquad
Y(r,t)\equiv \frac{\sqrt{2\pi G}\,r}{\alpha}\,\phi_{,t},\qquad X_\pm = X \pm Y.
$$

`X_±` propagate along the (null) characteristics, controlled by `g ≡ a/α`. Complete first-order system
(**G97 eqs 5–8**), with `C_1, C_2, C_±` defined for convenience:

$$
r\,g_{,r} = (1-a^2)\,g, \tag{5}
$$
$$
r\,a_{,r} = \tfrac{1}{2a}\big[(1-a^2)+a^2(X_+^2+X_-^2)\big] \equiv C_1, \tag{6}
$$
$$
g\,r\,a_{,t} = \tfrac{1}{2}a^3(X_+^2 - X_-^2) \equiv C_2, \tag{7}
$$
$$
r(X_{\pm,r}\mp g\,X_{\pm,t}) = \Big[\tfrac12(1-a^2) - a^2 X_\mp^2\Big]X_\pm - X_\mp \equiv C_\pm. \tag{8}
$$

Ricci scalar (diagnostic): `R = 4 r^{-2} X_+ X_-` (G97 eq. 9). Basic variables: `X_+, X_-, g, a`.

> **Alternative variable set (Choptuik/LRR).** [LRR eqs 6–10] use `Φ ≡ φ_{,r}` and `Π ≡ (a/α)φ_{,t}`
> with `a_{,r}/a + (a^2-1)/(2r) - 2\pi r(\Pi^2+\Phi^2)=0` (Hamiltonian constraint) and
> `\alpha_{,r}/\alpha - a_{,r}/a - (a^2-1)/r = 0` (polar slicing). Equivalent to G97's `X_±` up to the
> definitions above; use whichever the solver prefers. `[VERIFIED — LRR §2.1.1]`

## 2.2 DSS definition and similarity coordinates  `[VERIFIED — G97 eqs (10)-(15)]`

**Continuous** self-similarity (CSS/homothety): a vector field `χ` with `£_χ g_{ab} = 2 g_{ab}` (eq. 10).

**Discrete** self-similarity (DSS, **G97 eq. 11**): a diffeomorphism `φ` and real constant `Δ` s.t. for
all integers `n`,

$$
(\varphi^*)^n g_{ab} = e^{2n\Delta} g_{ab}.
$$

In adapted coordinates `(σ, x^α)` where `φ` maps `(σ,x^α)↦(σ+Δ,x^α)`, DSS ⇔ (**eq. 12**)

$$
g_{\mu\nu}(\sigma,x^\alpha) = e^{2\sigma}\,\tilde g_{\mu\nu}(\sigma,x^\alpha),\qquad
\tilde g_{\mu\nu}(\sigma,x^\alpha)=\tilde g_{\mu\nu}(\sigma+\Delta,x^\alpha).
$$

`σ` = logarithm of spacetime scale. CSS is the `Δ→0` degenerate limit.

**Specific similarity coordinates used for the eigenvalue problem (G97 eq. 15):**

$$
\boxed{\;\tau \equiv \ln\!\left(\frac{t}{r_0}\right),\qquad
\zeta \equiv \ln\!\left(\frac{r}{t}\right) - \xi_0(\tau)\;}
$$

where `r_0` is an arbitrary fixed scale and **`ξ_0(τ)` is a periodic function of period `Δ`, to be
determined as part of the solution.** (`τ` here is the DSS log-time; the singularity is approached as
`τ → −∞` when `t→0^+` from the past light cone. Sign of `Y` flips if you reverse the collapse
convention.)

> **Map to the "textbook" scalar coordinates.** The commonly-quoted form is
> `τ = −ln(T_*−t)`, `x = r/(T_*−t)` (or `ln[r/(T_*−t)]`). G97's `(τ,ζ)` with the extra periodic shift
> `ξ_0(τ)` is a **gauge-refined** version of exactly that: `ξ_0` is chosen (eq. 24 below) to pin the
> singular surface to a coordinate line, which is what makes the numerics tractable. The `T_*−t` and
> `t/r_0` differ only by the choice of time origin/scale. `[Textbook form: standard, structurally
> matches; the precise `ξ_0` refinement is G97-specific — VERIFIED as G97's actual choice.]`

## 2.3 Reduced field equations in similarity coordinates  `[VERIFIED — G97 eqs (16)-(19)]`

$$
X_{\pm,\zeta} = \frac{C_\pm \pm e^{\zeta+\xi_0} g\, X_{\pm,\tau}}{1 \pm (1+\xi_0')\,e^{\zeta+\xi_0} g}, \tag{16}
$$
$$
g_{,\zeta} = (1-a^2)\,g, \tag{17}
$$
$$
a_{,\zeta} = C_1, \tag{18}
$$
$$
a_{,\tau} = e^{-(\zeta+\xi_0)} g^{-1} C_2 + (1+\xi_0')\,C_1, \tag{19}
$$

with `ξ_0' ≡ dξ_0/dτ`. Invariant under `τ`-translation (= change of `r_0`). The scalar field itself is
recovered from (**G97 eq. 22**)

$$
\phi_{,\tau} = (2\pi G)^{-1/2} a\big[(1+\xi_0')X + g^{-1}e^{-(\zeta+\xi_0)}Y\big],\qquad
\phi_{,\zeta} = (2\pi G)^{-1/2} a X,
$$

which forces `φ = (periodic in τ) + κτ` (eq. 23); **boundedness ⇒ `κ = 0`** — Gundlach imposes this
(the Choptuik attractor has bounded `φ`).

## 2.4 Boundary conditions + the singular point (sonic line / SSH)  `[VERIFIED — G97 eqs (20)-(25)]`

**Periodicity (DSS) + center regularity + gauge (eqs 20–21):**

$$
a(\zeta,\tau+n\Delta)=a(\zeta,\tau),\quad g(\zeta,\tau+n\Delta)=g(\zeta,\tau),
$$
$$
a(\zeta=-\infty,\tau)=1,\qquad g(\zeta=-\infty,\tau)=1 \quad(\text{regular center }r=0,\ \alpha=1).
$$

**Coordinate condition fixing the singular surface (G97 eq. 24).** The equations (16) become singular
where the denominator vanishes. Gundlach uses the freedom in `ξ_0(τ)` to put this on the line `ζ = 0`
"for all `τ` at once":

$$
\boxed{\;\big[\,1-(1+\xi_0')\,e^{\xi_0} g\,\big]_{\zeta=0}=0.\;}
$$

This line `ζ=0` is **null** — the **past self-similarity horizon (SSH)**, a.k.a. the **sonic line**
(the DSS analogue of the fluid sonic point). It is a Cauchy horizon for the periodic problem.

**Regularity (analyticity) condition at the SSH (G97 eq. 25) — the selection condition:**

$$
\boxed{\;\big[\,C_- - e^{\xi_0} g\, X_{-,\tau}\,\big]_{\zeta=0}=0.\;}
$$

Gundlach proves (via the local solution eqs 26–30, homogeneous mode `∝ F(...)`) that imposing this one
condition **automatically enforces analyticity** (the non-analytic homogeneous piece has `F ≡ 0`).

## 2.5 Construction of the DSS solution (the method)  `[VERIFIED — G97 §II.C]`

The DSS critical solution is obtained as a **nonlinear eigenvalue problem / periodic hyperbolic BVP**:

1. **Interchange space and time.** Treat `ζ` as the evolution ("time") coordinate and the periodic `τ`
   as the "space" coordinate. Eqs (16)–(18) are then first-order **evolution equations in `ζ`** for
   `(a, g, X_+, X_-)`, subject to the one constraint (19) (which propagates).
2. **Fourier in `τ`.** Because everything is periodic in `τ` and smooth, expand every field
   (`a, g, X_±, ξ_0`) in a **truncated Fourier series in `τ`** (rapid convergence). Reconstruct `a` at
   each `ζ`-step from the constraint (19) rather than evolving it.
3. **Impose a reflection symmetry** (κ=0 ⇒ specific even/odd frequency content): `a, g, ξ_0` carry only
   **even** frequencies in `τ`; `X_±` only **odd** frequencies.
4. **Two BCs at `ζ=−∞`** (regular center: `a=1, g=1`) and **two at `ζ=0`** (the coordinate condition
   eq. 24 + the analyticity condition eq. 25). Degree-of-freedom count (`4∞−∞−1+1+∞=4∞`) balances,
   so the solution set is discrete.
5. **Unknowns solved for:** the periodic Fourier fields **plus the period `Δ` itself plus the periodic
   function `ξ_0(τ)`** — solved by relaxation. Gundlach finds one locally-unique solution; it matches
   Choptuik's attractor to numerical precision.

(§III of G97 extends this to the future SSH via further coordinate changes `w, ρ` — **not needed for the
exponent**; skip for a `γ`-only solver.)

## 2.6 Perturbations → critical exponent γ  `[VERIFIED — G97 §IV, eqs (64)-(77)]`

Linearize about the DSS background `Z_* = (a,g,X_+,X_-)_*`. Because the background is `τ`-periodic but a
generic perturbation is **not**, expand (**G97 eq. 64**):

$$
\delta Z(\zeta,\tau) = \sum_{i} C_i\, e^{\lambda_i \tau}\,\delta_i Z(\zeta,\tau),\qquad
\delta_i Z\ \text{periodic in }\tau\ \text{with period }\Delta,\ \ \lambda_i\in\mathbb{C}.
$$

(Only `0 ≤ Im λ_i < π/Δ` need be considered; eigenvalues come in conjugate pairs.) If the linearized
background equation is `δZ_{,ζ} = A\,δZ + B\,δZ_{,τ}` (eq. 65), each Floquet mode obeys

$$
\delta_i Z_{,\zeta} = (A + \lambda_i B)\,\delta_i Z + B\,\delta_i Z_{,\tau}. \tag{66}
$$

(The explicit four-component linear system is G97 eqs 67–70.) Boundary conditions: regular center at
`ζ→−∞`, and at `ζ=0` the perturbed analyticity condition (vanishing numerator of the `δ_iX_-`
equation). `λ_i` are the eigenvalues of this **linear periodic hyperbolic BVP** (discrete set). Then the
black-hole-mass scaling (dimensional analysis, eqs 71–81) gives

$$
\boxed{\;\gamma = -\frac{1}{\lambda_1} = \frac{1}{\operatorname{Re}\lambda_0}\;}\qquad(\text{G97 eq. 77})
$$

where `λ_1` is the single mode in the **left half-plane** (growing toward the singularity `τ→−∞`).
Gundlach finds **exactly one** such growing mode (plus the `λ=0` gauge mode `δZ ∝ Z_{*,τ}`).

**DSS "wiggle" (bonus, G97 §IV.B / LRR eq. 40).** In the DSS case the clean power law acquires a
universal periodic modulation:

$$
\ln M = \gamma\,\ln(p-p_*) + c + f\big[\gamma\ln(p-p_*)+c\big],\qquad f(z)=f(z+\Delta),
$$

period in `ln(p−p_*)` is `Δ/(2γ) ≈ 4.61`. Not needed for `γ` itself.

## 2.7 Verified numbers  `[VERIFIED]`

| quantity | value | source |
|---|---|---|
| echoing period `Δ` | **`Δ = 3.4453 ± 0.0005`** | G97 abstract + §II.C (verbatim) |
| growing-mode Lyapunov exponent `λ_1` | **`λ_1 = −2.674 ± 0.009`** (real) | G97 §IV.A + §V.C |
| critical exponent `γ = −1/λ_1` | **`γ = 0.374 ± 0.001`** | G97 abstract + §V.C |
| `γ` from best collapse simulation (agreement) | `≈ 0.374` | G97 §V.C (cites Choptuik) |
| CSS/DSS cross-check (review) | `γ ≈ 0.37`, `Δ ≈ 3.44` | LRR §2.1, §3.3 |

Confirmed: **`γ = 0.374 ± 0.001`** and **`Δ = 3.4453 ± 0.0005`** (Gundlach). (G97 supersedes his own
earlier `Δ = 3.4439 ± 0.0004` from the Letter [16]; use **3.4453**.)

## 2.8 The one gap a builder must close for Priority 2

**DSS is a periodic BVP with the period `Δ` and the gauge function `ξ_0(τ)` themselves unknowns** — a
Fourier-spectral relaxation solver, materially harder than the CSS shoot. G97 gives the equations and
the *scheme* (Fourier-in-τ, evolve in ζ, reconstruct `a` from the constraint, relax on `Δ, ξ_0` +
Fourier coefficients under BCs 21/24/25) but the **appendices (A–G) with the actual Fourier machinery,
the local expansions at the singular line, and the numerical relaxation details were not read here** (only
the main text pp.1–15 was parsed). Before coding the DSS solver, **fetch and read G97 appendices C, E, F**
(Fourier method, perturbation numerics, error estimate). Recommendation: **prove the method on the CSS
fluid (Priority 1) first** — same conceptual pipeline, one dimension simpler, known answer β≈0.3558 —
then step up to DSS.

---

# Single biggest gap/uncertainty overall

For **Priority 1 (the warm-up oracle)**: the autonomous ODE system, the three boundary conditions, the
sonic-point selection principle, the perturbation ansatz `h_p(x)e^{κs}`, and the numbers
(`κ₀≈2.81055255`, `β≈0.35580192`) are **fully verified and codeable**. The one derivation the builder
must perform (not printed in KHA) is the **explicit sonic-point regularity equation + the L'Hôpital
local expansion at `x=0`** — obtained by writing the 2×2 solve for `(ω_,x, V_,x)`, setting its
determinant to zero, and demanding numerators vanish there; **Ori & Piran PRD 42 (1990) 1068** is the
worked reference if needed. Secondary check: **`[VERIFY-BRACES]`** in §1.2 — confirm the brace balance
of the second fluid equation directly against gr-qc/9503007 p.3.

---

## Provenance / reproducibility

- Papers fetched as PDFs, converted to text with a local PDF→text tool, and read line-by-line:
  gr-qc/9503007 (KHA, 4 pp., all equations transcribed above); gr-qc/9604019 (G97, pp.1–15 read;
  appendices A–G **not** read); gr-qc/0001046 (LRR/Gundlach 2000 preprint, read for cross-check of
  numbers, coordinate conventions, and the `M ∝ (p−p_*)^{1/λ_0}` / DSS-wiggle framing).
- Not independently digit-verified: Maison's perturbative β (LRR asserts agreement with KHA);
  Evans–Coleman's full ODE system (only the β≈0.36 evolution result and the ξ=r/t^n convention were
  confirmed, via abstract + LRR §3.2, not the EC PDF body).
- Every boxed equation and every table number traces to a fetched primary source as cited inline.
