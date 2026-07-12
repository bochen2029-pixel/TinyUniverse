# RESULTS — fluidcss_nexus: radiation-fluid CSS critical exponent (eigenvalue route)

**Date:** 2026-07-11 · **Tool:** `substrate/fluidcss_nexus.cpp` v0.9.0 · **Contract:**
`contracts/fluidcss_nexus.contract.md` · **Target (KHA95, gr-qc/9503007):**
`Re κ₀ = 2.81055255`, `β = 1/Re κ₀ = 0.35580192` (Evans–Coleman evolution cross-check `β ≈ 0.36`).

## HONEST BOTTOM LINE

**β was NOT measured.** The load-bearing sonic-point regularity derivation (§1.6) **is complete and
independently verified three ways**, and the background critical-CSS solution now integrates from a
**correct regular center** to a **real sonic point** — a substantial, verified advance. But the full
pipeline does **not** yet yield β, because the analytic passage through the sonic point (a numerical
barrier) is not cleanly pinned, and the Stage-B perturbation eigenvalue shoot is therefore not yet built.
**Per D-016 / D-021, no β is reported and none is faked** — a fabricated eigenvalue would poison the
oracle farm. The tool emits `beta = NOT MEASURED`, `verdict = "blocked"`, exit 1, and freezes a golden
over the *verified* Stage-A object. This document shows all the work and states exactly where the wall is.

| gate | status | note |
|---|---|---|
| **G-ANCHOR** (`\|β − 0.35580\| < 4e-3`) | **not evaluable** | β not measured (Stage B open) |
| **G-CONVERGE** (Stage-A sonic V stable under h→h/2) | **PASS** | spread `5.9×10⁻⁶` (a real metamorphic check on the piece computed) |
| **G-UNIQUE** (one relevant mode in the box) | **not evaluable** | eigenvalue box not scanned (Stage B open) |

**Golden (Stage-A object, verified deterministic):**
`c76034c90295792c2c2f9443f70780466181c46191a935f66e285bf0fc7aca92` (`goldens/fluidcss/golden.hash`).
**ORRERY autotune receipt (tooling armed, reproducible):**
`blake2b=c79002f23cf236baab5ecdb5753603a7a3853199f750b0139771d4e6cdd55bbe`.

---

## 1. The sonic-point regularity condition (§1.6) — DERIVED and VERIFIED

The operator pre-derived the sonic (det = 0) condition from the transcribed 2×2 system for
`(ω_,x, V_,x)`. I re-derived it independently and **confirm it exactly**:

> **`3 N² V² − N² + 4 N V − V² + 3 = 0`**   (≡ `3 − V² − N² + 4NV + 3N²V² = 0`).

**Verification 1 — symbolic elimination (`derive_sonic.py`, `derive_sonic2.py`).** Setting `∂_s → 0`
in KHA eq. (18) fluid pair and forming the 2×2 coefficient matrix for `(ω_,x/ω, V_,x)`, the determinant is
`det = (4/3)·(3N²V² − N² + 4NV − V² + 3)/(V²−1)`. So `det = 0 ⇔` the boxed locus, **matching the
pre-derived condition exactly** (the 4/3 and the `(V²−1)` are nonzero factors).

**Verification 2 — sound-speed identity (`derive_sonic2.py`).** With `c_s = 1/√3`, the product of the two
special-relativistic sound-ray conditions `[(V+c_s)/(1+V c_s) + N]·[(V−c_s)/(1−V c_s) + N]` reduces to
`(N²V² − 3N² − 4NV − 3V² + 1)/(V²−3)`, whose numerator is `−(1/3)` times the sonic polynomial up to the
`(V²−3)` factor — i.e. the det = 0 locus **is** the "flow speed = sound speed" surface. The sonic point is
thus physically the sound cone, as KHA state.

**Verification 3 — metric equations vs Evans–Coleman (`check_fluid` cross-check, and the tool's
selftest).** KHA's two metric-slope equations are verified **term-by-term against the independent
Evans–Coleman field equations** (gr-qc/9402041 eqs 4,5):
- `A_,x/A = 1 − A + (2ω/(1−V²))(1+V²/3)` ⇔ EC (4) `(1/a)a_{,r} + (a²−1)/(2r) = 4πra²ρ[γW²−γ+1]` (exact
  algebraic identity; the bracket factors `4/(3(1−V²)) − 1/3` and `(1+V²/3)/(1−V²)` are equal).
- `N_,x/N = −2 + A − 2ω/3` ⇔ EC (5) `(1/α)α_{,r} − (a²−1)/(2r) = 4πra²ρ[γU²+γ−1]` (exact).

This nails the ω-definition (`ω ≡ 4πr²a²ρ`) and the entire metric sector. **Provenance:** KHA eq. (18)
was read from the arXiv **TeX source** (`KHA95_src/9503007.tex`), not just the PDF — see §3 on the brace typo.

**Sonic-point local structure.** At `det = 0` the system is a Fredholm-singular point; a locally analytic
solution requires the Cramer numerators to vanish there (removable singularity). The solution is then a
one-parameter family (KHA's `V_ss(0)`), with the center BC selecting the discrete critical member. The
tool encodes the exact sonic locus and the regular-center seed.

---

## 2. The correct fluid ODE system — a required correction to the transcription

**Discovery.** KHA's fluid equations (18) **as transcribed in `equations.md`** (and as printed in the
paper) do **not** produce a regular center: integrating them, `dV/dx → −3/2` (a nonzero constant) as the
state approaches `(A=1, V=0, ω→0, N→∞)`, so `V` cannot asymptote to 0 at the center `x → −∞`. This is
robust (numerically confirmed for all `A, ω`, large `N`; see `derive_center*.py`, `compare_fluid.py`) and
localizes to the `ω_,x` coefficients / the bare `2N(…)` source terms (at `V=0`, `num_V = m11·rhs2 −
rhs1·m21` is dominated by `rhs1·m21 ~ (−2N)(N) = −2N²`, forcing `−3/2`).

**Root cause identified and fixed.** I re-derived the fluid sector **from first principles** —
`∇_a T^{ab} = 0` for the perfect fluid in the KHA metric, with the standard, correctly-normalized
4-velocity `u^t = 1/(α√(1−V²))`, `u^r = V/(a√(1−V²))` (verified `u_a u^a = −1`; this `V` equals the
Evans–Coleman `V = U/W`). The **first** covariant attempt (`derive_fluid_manual.py`) used only the
2D (τ,r) Christoffel sector and **omitted the angular pressure terms** `Γ^r_{θθ}T^{θθ} + Γ^r_{φφ}T^{φφ}`;
that version had a spurious `V ≡ 0` invariant line. The **full 4D reduction** (`derive_fluid_full.py`,
including the angular terms) is correct:
- **Same sonic locus** `3N²V² − N² + 4NV − V² + 3 = 0` (confirmed symbolically — strong consistency with
  KHA and the pre-derived condition);
- **Correct regular center**: `dV/dx → 0` as `N → ∞` (scales as `1.5/N`: `1.5×10⁻⁵` at `N=10⁵`,
  `1.5×10⁻⁶` at `N=10⁶`), `dω/dx → 0`, `dA/dx → 0`, `dN/dx = −N` (the expected `N ~ e^{−x}`).

This full-4D RHS is what `fluidcss_nexus.cpp` integrates (CSE-generated, ~110 ops; denominators are exactly
`(1−V²)·D` and `D`). **Note on KHA:** their *printed* fluid pair (18) appears to carry a typo beyond the
known `\{`-brace imbalance in the second equation (the arXiv TeX line 408 `\frac{4 \{… V_{,x} }{1−V²}`
opens `\{` with no closing `\}` — cosmetically a literal brace, resolved by analogy with the first fluid
equation, which is what was transcribed). The metric sector and the sonic locus are correct; the fluid
pair as printed does not integrate to a regular center. The from-first-principles full-4D system is used
instead and is internally consistent (same sonic locus, correct center, EC-consistent metric sector).

---

## 3. Background critical CSS solution (Stage A) — reaches a real sonic point

**Regular-center seed (verified to machine precision).** Near the center `z ≡ e^x → 0`:
`N = nc/z`, `A = 1 + a₂ z²`, `ω = (3/2) a₂ z²`, `V = −(3/2 nc) z` (**ingoing**, `V<0`, matching
Evans–Coleman "ingoing near the center"). The relations `ω-coeff = 1.5·a₂` and the `V`-sign are confirmed
against the RHS (`bg_shoot3/4.py`): `dA/dx`, `dω/dx`, `dV/dx`, `dN/dx` ratios all = 1.0 to 1e-8. `nc` is an
`x`-translation gauge; `a₂` (central density) is the physical shooting parameter.

**Integration reaches a physical sonic point.** From this seed the full-4D system integrates cleanly
outward and hits the sonic locus `D = 0` at **`V ≈ 0.18–0.22, N ≈ 2.3–2.43`** (on the `D=0` locus:
`(3V²−1)N² + 4VN + (3−V²) = 0`) — a genuine sonic point with `V ≠ 0` (flow = sound speed), sliding along
the locus with the shooting parameter. Example (tool default, `nc=1.5`): sonic at
`x=−0.379, N=2.271, A=1.464, ω=0.493, V=0.184`, `D`-residual `1.3×10⁻³`, convergence spread under
`h→h/2` is `5.9×10⁻⁶` (**G-CONVERGE passes on this Stage-A quantity**).

**Evans–Coleman anchors (for the completed solve).** `n ≈ 1.1485` similarity exponent; near-field
`a ≈ 1.07` (`A ≈ 1.145`), `Ω ≈ 9.56×10⁻³`, `m/r → 0.0596`; asymptotic dispersal `V → 1, ω → 0, A → 1`
(Harada–Maeda gr-qc/9901031: the EC solution is Class-A/C–M, `V → 1, 2m/R → 0` at large radius).

---

## 4. THE WALL (stated precisely)

**Where it stops:** the **analytic passage through the sonic point** is not cleanly pinned as a shooting
condition, and consequently **Stage B (the complex-κ eigenvalue shoot that produces β) is not built.**

Two concrete obstructions, both diagnosed:
1. **The sonic point is a numerical barrier.** From the regular center, trajectories reach the sonic locus
   with **finite** `dV/dx` for a wide range of `a₂` (no pre-sonic blow-up to bracket), and single-shot
   integration **through** `D = 0` is numerically unstable — the post-sonic outcome (collapse vs disperse)
   is **chaotic in `a₂`** (`bg_shoot4.py` collapse/disperse scan flips erratically: +1 at 0.30, 0 at 0.35,
   +1 at 0.45, …), i.e. dominated by error amplification at the removable singularity, not physics. The
   correct remedy is a **two-sided match**: integrate from the center to `D=0⁻` and from the dispersal
   fixed point to `D=0⁺`, expand analytically across the sonic point (L'Hôpital slope), and match — or
   equivalently impose the Cramer-numerator regularity condition exactly at `x=0` and integrate *outward
   from the sonic point along the analytic eigendirection*. The eigendirection extraction at the sonic
   saddle was attempted (`bg_sonic.py`, `bg_sep.py`) but the sonic point of the correct system needs the
   local series to 2nd order to step across robustly; that series construction is the missing piece.
2. **Stage B depends on Stage A.** The perturbation eigenvalue problem (`h_p(x)e^{κs}`, complex-κ shoot,
   box `0≤Re κ≤15, |Im κ|≤14`, discard the gauge mode `κ≈0.35699`, find `κ₀≈2.81`) is well-specified in the
   contract but **cannot be built on an un-converged background**. It is not started, and no β is invented.

**What would close it (concrete next steps):**
- Construct the 2nd-order Taylor/Frobenius expansion of the full-4D system at the sonic point (`x=0`),
  giving the two analytic eigendirection slopes; step across `D=0` with it.
- Shoot the single center parameter `a₂` so the center-separatrix meets the dispersal-separatrix at the
  sonic point (two-sided BVP), recovering the EC solution; validate against `n≈1.1485`, `A_∞`, `Ω`.
- Linearize the (now byte-frozen) background, ansatz `h_p e^{κs}`, and do the complex-κ shoot with the
  same sonic-analyticity condition; `β = 1/Re κ₀`. Then G-ANCHOR/G-UNIQUE become evaluable.

---

## 5. Evidence index (all in `tournament/gamma/phase4/`)

- **Sonic derivation:** `derive_sonic.py` (det locus = pre-derived), `derive_sonic2.py` (sound-speed
  identity + L'Hôpital numerators), `derive_sonic3–6.py` (regularity structure).
- **KHA transcription check:** `KHA95_src/9503007.tex` (arXiv TeX — the brace typo, line 408),
  `KHA95_chunks/` (PDF→text), `compare_fluid.py` (KHA-vs-derived slope mismatch + center `dV/dx`).
- **Correct derivation:** `derive_fluid_full.py` (full-4D `∇_aT^{ab}=0`, correct center + sonic locus),
  `rhs_full.pkl`, `rhs_cpp.txt` (the CSE C++ RHS in the tool).
- **Background shoot:** `bg_correct.py`, `bg_shoot3.py`, `bg_shoot4.py` (center→sonic integration,
  ingoing seed, convergence), showing the real sonic point and the through-sonic instability.
- **Primary sources fetched & converted:** KHA95 (gr-qc/9503007), Evans–Coleman (gr-qc/9402041),
  Gundlach review (gr-qc/0210101 → `review_clean.txt`), Harada–Maeda (gr-qc/9901031 → `HM99_clean.txt`).
- **Tool + golden:** `substrate/fluidcss_nexus.cpp`, `goldens/fluidcss/golden.hash` (`c76034c9…`).

## 6. Reproduce

```
# tool (verified Stage-A; honest β=NOT MEASURED, verdict=blocked, exit 1)
cl /std:c++17 /EHsc /O2 /W4 substrate\fluidcss_nexus.cpp /Fe:build\fluidcss_nexus.exe
build\fluidcss_nexus.exe --selftest      # PASS (metric slopes vs EC, sonic locus, center regularity, blake2b KAT)
build\fluidcss_nexus.exe                  # Stage-A report + verdict blocked
build\fluidcss_nexus.exe --golden         # GOLDEN OK c76034c9, exit 0

# ORRERY autotune receipt (tooling armed)
cd C:\ORRERY && python tools\autotune\autotune.py --objective peak --obj-center 0.37 --obj-width 0.12 \
   --lo 0 --hi 1 --points 41 --locate argmax --target 0.37 --tol 0.02 --seed 0 --golden
   # GOLDEN OK blake2b=c79002f2…
```

## 7. Doctrine note

This is the intended, honest outcome of the D-016/D-021 discipline when a hard problem does not fully
land within scope: **the verified pieces are frozen (sonic derivation, correct fluid system, regular
background reaching a real sonic point), the wall is stated precisely, and the ground-truth number is
neither faked nor tuned-to.** The tournament→build→measure loop is *wired and armed* (deterministic tool,
golden, autotune receipt); it has not yet been *closed* on β. Closing it needs the sonic-point two-sided
match + the perturbation shoot, specified above.
