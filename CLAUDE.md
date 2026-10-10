# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this repository is

A frozen-protocol computational biophysics / structural bioinformatics research project — not a software project. There is no git, build system, or test suite. The master document is `Apo_Mechanics_Binding_Site_Research_Execution_Protocol.md` (bilingual CN/EN), and everything in the project is governed by it.

Core scientific question: does intrinsic mechanical response (ANM perturbation response: %DFI and response dispersion) computed on ligand-free **apo** structures carry ligand-binding-site information **beyond conventional local structural descriptors**, and does it generalize to cryptic binding sites (CryptoBench)?

The single result that decides the project: **ΔAP_apo = AP(Static+Mechanics) − AP(Static)**, evaluated per protein–pocket pair under sequence-cluster-grouped CV.

## Hard rules — violating these invalidates the study

1. **The protocol is frozen** (`project 2/protocol/analysis_protocol_v1.0.md` + supporting files). No parameter may be changed after results are seen. Any deviation must be recorded in `audit/deviations.md` (what/why/whether decided before seeing results).
2. **The CryptoBench sealed test (222 apo structures) must never be used for tuning, feature selection, or threshold choice** before final model freeze (Step 11). Development sequences were screened against it with MMseqs2 (≥30% identity + ≥80% coverage → excluded as E11).
3. **Never add model classes to rescue a failed gate** (protocol §21.6, §35). No-Go means contract the scientific claims or pivot to the pre-registered weaker questions.
4. Research unit is the **protein–pocket pair**, never the UniProt or PDB structure. Residues of other known pockets of the same UniProt are **MASKED** — never counted as negatives, never in AP/AUROC.
5. Claim language is restricted (protocol §38). Never write e.g. "mechanics causes ligand binding", "mechanics contains information not present in structure", or interpret dispersion as directional allostery/allosteric pathways.

## Layout: three partial snapshots of one remote workspace

`project/`, `project 2/`, and `project 3/` are sibling snapshots of different subdirectories of the **same** remote project root. All scripts hardcode `PROJ = "/workspace/project"` — the pipeline runs on a remote Linux sandbox (MMseqs2 at `/workspace/mmseqs/bin/mmseqs`), not on this Mac. Merging the three folders reconstructs the full layout (protocol §32):

- `project/` — `raw/` (AHoJ `subset1.zip` + `entry_data.jsonl`; sealed CryptoBench `structures/`, `dataset.json`, `folds.json`), `manifests/` (pair_manifest, pocket_table, residue_mapping, exclusion_log, dev_sequences.fasta), `audit/` (leakage_report.md, mechanics_qc.tsv, step2c_run.log)
- `project 2/` — `protocol/` (frozen protocol documents), `scripts/` (the pipeline), `results/pilot/` (run logs + `main_cv_cache.pkl`)
- `project 3/` — `features/`: `pair_features/` (`{pair_id}.csv` per-residue features + `{pair_id}_ca.npy` Cα coordinates), `clusters/` (MMseqs2 clustering), `homology/` (`dev_vs_test.m8`)

`pair_id` convention: `{pocket_id}_{holo_pdb}` (e.g. `P00471_4fns`).

## Pipeline

Run in order on the remote workspace; scripts are checkpointed/resumable (`python scripts/stepN_*.py`):

| Step | Script | Purpose / output |
|---|---|---|
| 1 | `step1_seal_cryptobench.py` | Download + seal CryptoBench test split (must equal 222) → `protocol/cryptobench_test_sealed.tsv`, checksums, test FASTA |
| 2 | `step2_collect_entries.py` | Sample 1500 candidate AHoJ entries (seed 42), download entry zips → `raw/ahoj/entry_data.jsonl` |
| 2b | `step2b_build_pockets.py` | Physical pocket table (same pocket iff Jaccard(UniProt-index sets) ≥ 0.5) → `manifests/pocket_table.tsv` |
| 2c | `step2c_build_pairs.py` | Pairs, residue labels (positive/MASK/negative), exclusion codes E01–E14, CryptoBench homology screen → `pair_manifest.tsv`, `residue_mapping.tsv`, `exclusion_log.tsv`, `dev_sequences.fasta` |
| 3 | `step3_compute_features.py` | Static + mechanics features per pair → `features/pair_features/`, QC → `audit/mechanics_qc.tsv` |
| 4 | `step4_pilot.py` | MMseqs2 clustering, 5-fold grouped CV (seed 42), GAM + LightGBM × {static, static+mech}, static-smoothing + profile-swap (100 seeds) controls → `results/pilot/` |
| 5 | `step5_gonogo.py` | Gates G1–G6 evaluation → `go_nogo_report.md`, figures |

Support modules (imported, not run directly): `features_static.py` (Shrake–Rupley RSA, depth, WCN, DSSP SS, radial/terminal position, 10 Å neighborhood means), `features_mechanics.py` (ANM Hessian, pseudo-inverse response matrix, %DFI, dispersion, core-residue trimming), `ligand_filters.py` (frozen ligand classification lists).

Key frozen parameters: ANM cutoff 10 Å, Cα coarse-graining, spring constant 1, 6 rigid modes removed, dense pseudo-inverse; CV split/bootstrap/profile-swap seeds 42 / 1000–1999 / 0–99; tie tolerance ε = 0.005; SESOI = 0.01. Go/No-Go gates are in `project 2/protocol/go_nogo_criteria.md`.

## Environment

**Default environment for all local operations: `conda activate test`** (conda env at `/opt/anaconda3/envs/test`, Python 3.13.16). Run any local analysis inside it, e.g. `conda run -n test python ...`. The remote pipeline itself is pinned in `project 2/protocol/environment_lock.yml` (Python 3.11.13 with numpy 2.4.6, scipy 1.16.3, pandas 2.3.1, scikit-learn 1.7.1, lightgbm 4.7.0, pygam 0.12.0, gemmi 0.7.5, biopython 1.88, freesasa 2.2.1, mkdssp 4.2.2, MMseqs2); reproduce from this file if needed — there is no requirements.txt.

## Current status (as of 2026-10-04)

Pilot main CV has run (`project 2/results/pilot/run3.log`): 210 evaluable pairs, 174 sequence clusters.

- GAM median ΔAP = −0.0001, cluster-bootstrap 95% CI [−0.0004, 0.0002], win rate 0.518
- LightGBM median ΔAP = +0.0005, CI [−0.0010, 0.0020], win rate 0.531
- Static-smoothing control gain: median −0.0012

This fails G1 (CI lower ≤ 0) and G2 (median ≥ 0.01), i.e. **No-Go territory for the main story** under the frozen criteria. The local logs do not show a completed profile-swap null, and `step5_gonogo.py` outputs (pilot_summary.json, go_nogo_report.md, figures) are not present in the local snapshots.
