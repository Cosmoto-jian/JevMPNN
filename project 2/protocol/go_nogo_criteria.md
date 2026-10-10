# Go / No-Go Criteria (Frozen)

Primary quantity: Delta AP_apo = AP(Static+Mechanics) - AP(Static), per protein-pocket pair,
evaluated within-pair over positive+negative residues, grouped 5-fold CV, cluster bootstrap CI.

## Gates (all must be evaluated on the identical fold assignment)

| Gate | Criterion | Pass condition |
|------|-----------|----------------|
| G1 | Statistical direction | cluster-bootstrap 95% CI lower bound of median Delta AP_apo > 0 |
| G2 | Minimum practical effect (SESOI) | median Delta AP_apo >= 0.01 |
| G3 | Protein-level consistency | win rate > 50% excluding ties (epsilon = 0.005) |
| G4 | Model-form confirmation | LightGBM median Delta AP > 0 AND CI lower bound > 0 |
| G5 | Pipeline integrity | real median Delta AP > profile-swap null 95th percentile; swap-null center ~ 0 |
| G6 | Spatial-context control | Static+Mechanics > Static-smoothed required for strong mechanics claim |
| G0a | Leakage audit | CryptoBench sealed-test contamination = 0; fold overlap = 0 |

## Decision rules

- All of G1-G6 pass -> **Go**: proceed to full AHoJ development study (Steps 9-12).
- G1-G3 pass but G4 CI crosses 0, or G6 weak -> **Conditional Go**: expand development sample for
  confirmation without changing feature definitions or the sealed test.
- Any of: median Delta AP < 0.01; CI lower <= 0 after reasonable sample expansion; LightGBM
  direction opposite to GAM; swapped mechanics produces similar gain; leakage not removable ->
  **No-Go**: stop the "mechanics improves apo->holo prediction" main story; fall back to
  descriptive / redundancy / subgroup questions. No additional model classes may be added.

## Post-pilot power analysis

Using pilot SD(Delta AP) and cluster-level variance, estimate the required number of
protein-pocket pairs / sequence clusters for the full study; write to
analysis/power_analysis_after_pilot.md with the git commit hash at that time.
