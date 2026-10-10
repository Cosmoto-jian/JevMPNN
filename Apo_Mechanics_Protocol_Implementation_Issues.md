# Apo Mechanics Protocol Implementation Issues

Date: 2026-10-04

Scope: audit of the current pipeline against `Apo_Mechanics_Binding_Site_Research_Execution_Protocol.md` and `project 2/protocol/analysis_protocol_v1.0.md`.

## Executive Summary

The current implementation is not yet a strict execution of the frozen protocol. It has run part of the pilot workflow, and the cached main CV already shows a No-Go trend, but the result cannot be treated as a complete frozen-protocol result because several hard requirements are missing or implemented differently from the protocol.

Most important blockers:

1. CryptoBench homology exclusion is not implemented inside the reproducible pipeline.
2. MMseqs2 clustering does not use the frozen 30% identity / 80% coverage parameters.
3. The ANM response matrix does not match the protocol definition of `R_ij`.
4. Pilot did not complete profile-swap, collinearity audit, summary outputs, or Go/No-Go reporting.
5. `step5_gonogo.py` has an import-order bug and will fail before generating the report.

Until these are fixed and the pipeline is rerun from the affected stages, the current pilot should be treated only as an incomplete smoke/pilot attempt.

## Issue 1: CryptoBench Leakage Control Is Not Reproducible In The Pipeline

Severity: Critical

Protocol requirement:

- Development sequences must be screened against the sealed CryptoBench test sequences.
- Any development pair with sequence identity >= 30% and coverage >= 80% must be excluded as `E11_CRYPTOBENCH_HOMOLOG`.
- The result must be auditable through the manifest and leakage report.

Observed implementation:

- `project 2/scripts/step2c_build_pairs.py` documents a CryptoBench homology screen in the script header.
- The actual code does not run MMseqs2 `easy-search`.
- Included rows are written with empty `cryptobench_overlap_flag` and empty `sequence_cluster_id`.

Evidence:

- `project 2/scripts/step2c_build_pairs.py:14`
- `project 2/scripts/step2c_build_pairs.py:260`
- `project/manifests/pair_manifest.tsv` has empty `cryptobench_overlap_flag` and `sequence_cluster_id` fields for all included rows.

Additional observation:

- `project/audit/step2c_run.log` reports `212 included pairs, 188 exclusions`.
- Current manifests contain 210 included pairs and `exclusion_log.tsv` contains two `E11_CRYPTOBENCH_HOMOLOG` rows.
- This suggests the E11 exclusions were likely added by a post-processing step outside the recorded `step2c` program.

Impact:

- The strict no-leakage claim is not reproducible from the current pipeline.
- Any downstream CV result is not fully auditable under the frozen protocol.

Recommended fix:

- Add the MMseqs2 development-vs-sealed-test search to `step2c_build_pairs.py` or a dedicated frozen step.
- Use the exact protocol parameters: `--min-seq-id 0.3 -c 0.8 --cov-mode 0`.
- Populate `cryptobench_overlap_flag`.
- Write excluded pairs as `E11_CRYPTOBENCH_HOMOLOG` directly from the pipeline.
- Regenerate `pair_manifest.tsv`, `residue_mapping.tsv`, `exclusion_log.tsv`, and `dev_sequences.fasta`.

## Issue 2: MMseqs2 Clustering Does Not Use Frozen Parameters

Severity: Critical

Protocol requirement:

- Development clustering must use 30% sequence identity and 80% coverage.
- Frozen command parameters are documented in `project 2/protocol/sequence_cluster_rules.md`.

Observed implementation:

- `step4_pilot.py` calls `mmseqs easy-cluster` without passing `--min-seq-id 0.3`, `-c 0.8`, or `--cov-mode 0`.

Evidence:

- `project 2/protocol/sequence_cluster_rules.md:8`
- `project 2/protocol/sequence_cluster_rules.md:15`
- `project 2/scripts/step4_pilot.py:181`

Impact:

- The 5-fold grouped CV may not correspond to the frozen leakage-control definition.
- The reported cluster count and fold assignment are not guaranteed to satisfy the intended homology separation.

Recommended fix:

- Rebuild clusters with explicit frozen parameters.
- Write cluster IDs into `pair_manifest.tsv` or a dedicated split manifest.
- Regenerate the fold assignment and leakage audit.

## Issue 3: ANM Response Matrix Does Not Match Protocol Definition

Severity: Critical

Protocol requirement:

- `R_ij` must be the Frobenius norm of the 3x3 block `(i,j)` of the full pseudo-inverse `H+`.
- %DFI and response dispersion must be computed from that response matrix.

Observed implementation:

- `features_mechanics.py` constructs `W = V / sqrt(lambda)`, reshapes it into residue blocks, then computes `R2 = einsum("iab,jab->ij", A, A)`.
- This is not the same as taking the Frobenius norm of each 3x3 block of `H+`.
- The script does compute `Hplus` for QC, but it does not use `Hplus` to compute `R_ij`.

Evidence:

- `project 2/scripts/features_mechanics.py:49`
- `project 2/scripts/features_mechanics.py:53`
- `project 2/scripts/features_mechanics.py:56`

Impact:

- The central mechanical features may not be the frozen %DFI/dispersion features.
- Any AP gain or loss from mechanics may be about a different mechanical descriptor than the one preregistered.

Recommended fix:

- Compute `Hplus = (V[:, nz] / lam[nz]) @ V[:, nz].T`.
- For each residue pair `(i,j)`, compute:

```python
block = Hplus[3*i:3*i+3, 3*j:3*j+3]
R[i, j] = np.linalg.norm(block, ord="fro")
```

- Recompute all mechanics features and rerun Step 4.

## Issue 4: Pilot Did Not Complete All Frozen Gates

Severity: Critical

Protocol requirement:

- Pilot must evaluate:
  - GAM Static vs Static+Mechanics.
  - LightGBM Static vs Static+Mechanics.
  - Static spatial smoothing control.
  - Mechanical profile-swap null with 100 seeds.
  - Mechanics/WCN collinearity audit.
  - Go/No-Go report.

Observed outputs:

- `project 2/results/pilot/run3.log` stops after static-smoothing.
- No `pilot_summary.json`, `pilot_per_pair.tsv`, `cv_state.pkl`, or `go_nogo_report.md` was found in the local snapshots.

Observed pilot numbers before stopping:

- GAM median Delta AP = -0.0001, 95% CI [-0.0004, 0.0002].
- LightGBM median Delta AP = +0.0005, 95% CI [-0.0010, 0.0020].
- Static-smoothing median gain = -0.0012.

Impact:

- The main CV already indicates No-Go territory, but the formal gate set was not fully evaluated.
- The result cannot be described as a completed pilot under the frozen protocol.

Recommended fix:

- Fix upstream issues first.
- Rerun Step 4 from a clean cache.
- Confirm that profile-swap, collinearity, summaries, and per-pair outputs are produced.

## Issue 5: Go/No-Go Script Has A Runtime Bug

Severity: High

Observed implementation:

- `step5_gonogo.py` uses `csv.DictReader` before importing `csv` inside the function.
- Because `import csv` appears later inside the same function, Python treats `csv` as a local variable, causing `UnboundLocalError`.

Evidence:

- `project 2/scripts/step5_gonogo.py:19`
- `project 2/scripts/step5_gonogo.py:20`

Impact:

- Even if `pilot_summary.json` and `pilot_per_pair.tsv` existed, Step 5 would fail before report generation.

Recommended fix:

- Move `import csv` to the module-level imports.
- Remove the duplicated `per_pair` read.

## Issue 6: LightGBM Does Not Do Protocol-Specified Fold-Internal Hyperparameter Selection

Severity: Medium

Protocol requirement:

- LightGBM should use the same folds and feature sets.
- Limited hyperparameter selection must occur only inside the development fold.

Observed implementation:

- LightGBM uses fixed parameters:
  - `n_estimators=400`
  - `num_leaves=31`
  - `learning_rate=0.05`
  - `min_child_samples=20`

Evidence:

- `project 2/scripts/step4_pilot.py:9`
- `project 2/scripts/step4_pilot.py:98`

Impact:

- This is conservative in one sense, but it is not what the frozen protocol says.
- The confirmatory LightGBM gate is not implemented as specified.

Recommended fix:

- Either implement the preregistered limited inner-fold grid search, or formally record a protocol deviation before using the fixed-parameter result.

## Issue 7: Leakage Audit And Split Audit Are Incomplete

Severity: Medium

Protocol requirement:

- Each run should generate split/audit artifacts with:
  - train/test clusters,
  - fold leakage checks,
  - UniProt overlap = 0,
  - cluster overlap = 0,
  - CryptoBench sealed-test overlap = 0.

Observed implementation:

- `project/audit/leakage_report.md` contains only fold counts and cluster counts.
- It reports cluster overlap by construction but does not document sealed-test overlap or UniProt overlap.
- No `splits/run_x` folder or train/test cluster tables were found.

Evidence:

- `project/audit/leakage_report.md`
- `project 2/scripts/step4_pilot.py:200`

Impact:

- The analysis is not traceable enough for the protocol's audit standard.

Recommended fix:

- Write explicit split manifests for each fold.
- Include UniProt overlap, cluster overlap, and sealed-test overlap checks.
- Save the exact commands and hashes used for clustering.

## Issue 8: Manifest Fields Are Not Fully Populated

Severity: Medium

Observed data:

- `pair_manifest.tsv` includes `cryptobench_overlap_flag` and `sequence_cluster_id`, but both are empty for all included pairs.

Impact:

- Important audit fields exist structurally but do not carry the information required by the protocol.

Recommended fix:

- Populate these fields during Step 2c / Step 4.
- Avoid relying on external or manual post-processing.

## Issue 9: Duplicate Apo PDBs Exist Across Multiple Pockets

Severity: Low to Medium

Observed data:

- Several apo PDB IDs appear in multiple pair IDs, for example:
  - `1tcm`: `P00042_1tcm`, `P00044_1tcm`, `P00041_1tcm`
  - `5i7e`: `P00002_5i7e`, `P00001_5i7e`
  - `3amk`: `P00125_3amk`, `P00126_3amk`

Interpretation:

- This is not automatically wrong, because the research unit is protein-pocket pair and one protein may have multiple physical pockets.
- However, it requires strong masking and cluster grouping to avoid false negatives and leakage-like dependence.

Recommended fix:

- Add an audit table listing all repeated apo structures and confirming:
  - whether pockets are physically distinct,
  - whether other-pocket residues are masked,
  - whether all pairs from the same UniProt/cluster stay in the same CV fold.

## Current Pilot Status

The current cached pilot main CV produced:

| Model | Median Delta AP | 95% CI | Win Rate |
|---|---:|---:|---:|
| GAM | -0.0001 | [-0.0004, 0.0002] | 0.518 |
| LightGBM | +0.0005 | [-0.0010, 0.0020] | 0.531 |

Under the frozen Go/No-Go criteria, these values would fail:

- G1: CI lower bound must be > 0.
- G2: median Delta AP must be >= 0.01.
- G4: LightGBM CI lower bound must be > 0.

However, because the implementation is not yet protocol-faithful and the pilot did not complete all gates, this should be reported as:

> Incomplete pilot result showing No-Go tendency, not a completed frozen-protocol No-Go decision.

## Recommended Fix Order

1. Fix `features_mechanics.py` so `R_ij` is computed from the 3x3 blocks of `H+`.
2. Implement reproducible CryptoBench homology screening and E11 exclusion.
3. Rebuild MMseqs2 clusters with explicit frozen parameters.
4. Populate `cryptobench_overlap_flag` and `sequence_cluster_id`.
5. Regenerate manifests and features from the corrected pipeline.
6. Delete stale pilot cache and rerun Step 4.
7. Fix and rerun Step 5.
8. Generate the complete audit outputs required by the protocol.

## Bottom Line

The current pipeline does not yet strictly execute the frozen research protocol. The most serious problems are leakage-screen reproducibility, sequence-cluster parameter mismatch, and the central ANM response calculation. These should be corrected before interpreting the pilot result scientifically or making any project-level Go/No-Go decision.
