# Sequence Cluster Rules (Frozen)

## Development clustering

- Tool: MMseqs2 v18-8cc5c (static binary).
- Input: one representative sequence per development protein chain set (full observed construct,
  FASTA), one sequence per UniProt entry per distinct construct.
- Clustering: `mmseqs cluster` with `--min-seq-id 0.3 --coverage 0.8 -c 0.8`
  (coverage defined as fraction of the shorter sequence aligned, both sequence- and
  target-coverage set to 0.8 via `-c 0.8 --cov-mode 0` on preclustered identical sequences).
- Sensitivity default: `easy-cluster` default sensitivity is sufficient; no target-effort tuning.

## Sealed-test homology screen

- `mmseqs easy-search` development sequences vs CryptoBench sealed-test sequences,
  `--min-seq-id 0.3 -c 0.8 --cov-mode 0`, report alignment coverage of the shorter sequence.
- Exclusion rule: identity >= 30% AND coverage >= 80% -> exclude development pair (E11).

## Fold assignment

- 5-fold grouped CV; grouping key = MMseqs2 cluster id; all pairs of one cluster in one fold.
- Split seed 42; identical fold assignment reused by GAM, LightGBM, and all controls.
- Per-run audit: splits/run_x/{train,test}_clusters.tsv + leakage_check.tsv
  (UniProt overlap = 0, cluster overlap = 0, sealed-test overlap = 0).
