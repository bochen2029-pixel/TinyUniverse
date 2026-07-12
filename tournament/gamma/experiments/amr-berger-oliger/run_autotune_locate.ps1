# run_autotune_locate.ps1 — AMR/Berger-Oliger advocate, ORRERY-armed search receipt.
#
# Drives ORRERY `autotune` (Python, no GPU) to demonstrate a DETERMINISTIC level-crossing
# locate of a critical point against a PRE-REGISTERED --target, with a real gate. This stands
# in — as an evidence-grade proof-of-mechanism — for the p* level-crossing an AMR EKG tool
# would perform (does the pulse collapse? crosses 0->1 at p*). It does NOT compute gamma; the
# real "AMR EKG gear" does not exist yet (see amr-berger-oliger.md sec.4 for its contract).
#
# Emits: three declared JSON envelopes + their blake2b to autotune_provenance.json.
# Read-only use of ORRERY. Reproduce: pwsh -File run_autotune_locate.ps1
$ErrorActionPreference = "Stop"
Set-Location C:\ORRERY

Write-Host "== ORRERY autotune selftest ==" -ForegroundColor Cyan
python tools\autotune\autotune.py --selftest 2>&1 | Select-Object -Last 1

Write-Host "`n== Pre-registered locate: level-crossing 0.5 -> target 0.5 (stands in for p*) ==" -ForegroundColor Cyan
python tools\autotune\autotune.py --objective threshold --obj-center 0.5 --obj-width 0.05 `
    --lo 0 --hi 1 --locate crossing --level 0.5 --target 0.5 --tol 0.02 --seed 0 --json

Write-Host "`n== Provenance (declared blake2b for positive / neg-control / tight) ==" -ForegroundColor Cyan
$py = @'
import sys, json, hashlib
sys.path.insert(0, r"tools\autotune")
import autotune as at
def rec(A):
    r, gs, v, code = at.run_autotune(A, None)
    obj = at.declared_object(0, A, r, gs, v)
    h = hashlib.blake2b(obj.encode("utf-8"), digest_size=32).hexdigest()
    return {"x_located": r["x_located"], "located_error": r["located_error"],
            "on_target": r["on_target"], "verdict": v, "exit": code, "blake2b": h}
base = {"objective":"threshold","obj_center":0.5,"obj_width":0.05,"tool":"","sweep":"","metric":"",
        "fixed":"","lo":0.0,"hi":1.0,"points":41,"locate":"crossing","level":0.5,"target":0.5,"tol":0.02}
out = {
  "tool": "autotune", "orrery_version": "1.0.0",
  "positive":     rec(base),
  "neg_control":  rec({**base, "target":0.60}),
  "tight_tol":    rec({**base, "tol":0.005}),
}
print(json.dumps(out, indent=2))
'@
python -c $py | Tee-Object -FilePath "C:\TinyUniverse-tournament\tournament\gamma\experiments\amr-berger-oliger\autotune_provenance.json"
