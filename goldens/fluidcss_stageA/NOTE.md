# goldens/fluidcss_stageA — supersession note (D-032, 2026-07-19)

The original Stage-A golden `b4f4e46372d71aa3e768dc7a0e9de93eb163e27c4f13b42a66d8b93c92343503`
(frozen 2026-07-12, D-027, tool v0.9.1) banked the **collapsing flat FRIEDMANN solution
mislabeled as Evans-Coleman**: the v0.9.x construction used "V = -c_s" as the sonic
criterion, which is precisely the Friedmann point of the V0-parametrized sonic line
(rflanl.tex 4.7-4.9). Measured proof: V = -sqrt(1-1/A) holds to 1.9e-10 along that
background (the FRW free-fall identity); its "exact" invariants (A=1+2om/3, oi*=3/8,
2m/r=1/3) are all FRW radiation identities. See D-032 and
tournament/gamma/phase4/RESULTS_hka_beta.md (2026-07-19 addendum).

Superseded 2026-07-19 (operator-directed autonomous session) by the TRUE Evans-Coleman
background golden `27af792062bd9b8c...` (tool v1.0.0: the paper's own V0-shoot
construction; V0* = 0.1124394014, Nbar'(sonic) = -0.355699, one V-zero, two-sided
closure 4.2e-6). The old Friedmann construction is retained in the tool as the
`--friedmann` control face (closed-form exact solution; correctly labeled).
