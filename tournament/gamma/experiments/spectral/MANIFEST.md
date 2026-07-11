# experiments/spectral — artifact manifest (phase 1, spectral advocate)

All runs deterministic, CPU, no GPU. blake2b = 256-bit over the tool's canonical *declared object*
(ORRERY idiom: envelope minus `tool`/`version`/`notes`; the experiment idiom: `json.dumps(declared, sort_keys, sep=(",",":"))`).

## ORRERY autotune runs (read-only use of C:\ORRERY; CWD C:\ORRERY)

| # | command (abbrev) | role | x_located | verdict/exit | declared blake2b |
|---|---|---|---|---|---|
| G | `autotune.py --golden` | instrument anchor — tool prints "GOLDEN OK" | (peak@0.37) | pass / 0 | `c79002f23cf236baab5ecdb5753603a7a3853199f750b0139771d4e6cdd55bbe` |
| 1 | `--objective threshold --obj-center 0.5 --obj-width 0.05 --lo 0 --hi 1 --locate crossing --level 0.5 --target 0.5 --tol 0.02 --seed 0` | critical-point (level-crossing) locate vs pre-registered target | 0.500000 | pass / 0 | `db490a3111b770996fcec05b9bc59981f35395980be5f1ad617101dd13552c80` |
| 2 | `--objective peak --obj-center 0.374 --obj-width 0.05 --lo 0.2 --hi 0.55 --locate argmax --target 0.374 --tol 0.02 --seed 0` | p*-band (argmax) locate at the *actual* Choptuik value | 0.374014 | pass / 0 | `27c8cb4719f8a6e2dd455e386cd3da63ab1a9630896ab11a084a29c004d04eeb` |

autotune `--selftest`: all 6 checks PASS (blake2b KAT, argmax recovers center, crossing recovers center,
off-target fires G-OFF-TARGET → exit 1, determinism identical across two runs).

## Local numeric experiment (this dir)

| file | role | declared blake2b |
|---|---|---|
| `cheb_convergence.py` + `cheb_convergence.out.json` | D-021 rebuttal: Chebyshev collocation vs centered-2nd-order FD derivative of a sharp smooth pulse AND a C¹ cusp | `b6ad2eba341d427ffa8c825058cb6a383c4500660fabbf3f576b9920d63dcaaf` |

Headline result (`cheb_convergence`): on a smooth sharp Gaussian, Chebyshev derivative error falls
2.5 → 1.7e-2 → 1.6e-5 → 1.6e-9 → 6.7e-15 (fp floor by N=96); centered FD converges at measured order
−1.977 ≈ O(N⁻²). To match Chebyshev's N=64 error (1.6e-9), FD2 would need ≈622,000 points (~9,700× the DOF).
Honest flank: on a C¹-not-C² cusp, Chebyshev degrades to algebraic order −0.68 (Gibbs) — the exact failure a
naive fixed spectral grid hits near the singularity, and the reason the similarity/compactified map is mandatory.

**Status of every γ figure herein: `[ARGUMENT-GRADE]` / `[ARGUMENT-GRADE]` for the search structure.**
The real spectral EKG "EKG-spectral" gear does not exist yet; its contract is specified in `../../phase1/spectral.md` §4.
No γ number is measured here — autotune locates *analytic* targets to rehearse the search; the convergence law is
demonstrated on a *proxy* field, not the Einstein–Klein–Gordon system.
