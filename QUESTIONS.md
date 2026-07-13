# QUESTIONS.md — open design questions (parked, owned, non-blocking unless marked)

**Q-001 · Player embodiment. RESOLVED 2026-07-11 (default taken):** god-hand v1, frozen in frame.contract.md v1.0.0 (input-trace tool enum). Probe-mode ("your camera is a seat with a τ") remains a candidate M7+ mode — reopening it is an input-trace MINOR addition (new tool ids), not a breaking change.

**Q-002 · Particle species v1.** Current plan: one massive neutral species + photons. Charge/EM adds gameplay (plasma, lightning) at real cost (short-range solver). *Default: defer charge; keep the mass dial clean.*

**Q-003 · Product framing.** Standalone sandbox, or branded as the playable face of finaltheoryofeverything.org (the Ratchet/inscription framing is already the mechanic either way)? Affects naming/HUD copy only — zero physics impact. *No default; operator call, anytime before public release.*

**Q-004 · Quantum bubble budget on 16 GB.** 64³ fp32 complex ψ ≈ 2 MB/bubble (+FFT scratch); hundreds are plausible but unmeasured. If 3D is tight, sanctioned fallbacks: 2D bubbles (visually honest for thin scenes) or fewer/larger. *Resolved by measurement at M6 — no speculation.*

**Q-006 · The combined SR-inertia + 1PN-field precession coefficient. RESOLVED 2026-07-12 (D-026, `nexus/precession_nexus.cpp`).** The fp64 isolation experiments Q-006 owed, run: measured with the orbit's **actual** (force-distorted) semi-latus rectum, sr = 1.00π (exact Sommerfeld), pn1 = 6.03π (the correct GR 6π), **combined = 6.95π ≈ 7π — the linear superposition (π + 6π = 7π) is CONFIRMED.** The app's original 6.41π was a **normalization artifact**: it divided by the nominal (undistorted) semi-latus rectum p=6.40 instead of the orbit's actual p=6.90 (the forces expand the orbit ~8%); 6.95 × 6.40/6.90 = 6.45 ≈ 6.41π. So it was the *"implementation deviates from the declared model"* branch — the measurement normalization, not the superposition reasoning. (Separately, the combined model's 7π over-counts the correct full-GR 6π = pn1 = N3 `curve`; so retiring the 1PN-field-with-SR-inertia hybrid, D-016, stands — the app ships SR-only, exact, and the full 6π comes from N3.) *Closed.*

**Q-005 · Vulkan vs DX12.** Default Vulkan (cuda-samples parity, Mímir reference). Revisit only if Streamline/DLSS or tooling friction is measured to be materially worse than DX12. *One-way door closes at M1 freeze.*
