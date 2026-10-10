# Analysis Protocol v1.0 (Frozen)

Project: Does intrinsic mechanical response in ligand-free apo structures provide
ligand-binding-site information beyond conventional local structural descriptors?
Source protocol: Apo_Mechanics_Binding_Site_Research_Execution_Protocol.md (v1.0 Frozen Protocol).
This file restates every frozen parameter that governs computation. No parameter below may be
changed after Protocol Lock except through audit/deviations.md.

## 1. Research unit and labels

- Unit: protein-pocket pair (one physical pocket of one UniProt). Multiple holos of the same
  physical pocket may pool binding-residue labels (union within pocket). Never union across pockets.
- Positive (y=1): target-pocket lining residues defined by the holo structure and mapped to apo.
- Masked (y=MASK): residues of OTHER known physical pockets of the same UniProt. Excluded from
  AP/AUROC evaluation for this pocket; never negative.
- Negative (y=0): all other evaluable residues outside the target pocket and outside any other
  known pocket.

## 2. Data sources

- Development: AHoJ-DB (apoholo.cz), Subset 1 entry list (publication subset, v1 archive listing,
  downloaded 2026-10-03) + AHoJ-DB REST API for per-entry apo/holo mapping and pocket metadata.
- Sealed test: CryptoBench (OSF pz4a9), official folds.json test split = 222 apo structures,
  dataset.json v2024-12 (downloaded 2026-10-03). Sealed before any training; used only at Step 12.

## 3. Ligand scope (frozen)

- Biologically relevant small organic ligands; AHoJ/BioLiP2 relevance flag as first filter.
- Excluded: water; single-atom ions; PEG, glycerol and other crystallization additives; nucleic
  acids; peptide/protein ligands; obvious non-physiological buffer molecules.
- Minimum heavy-atom count: >= 6.

## 4. Structure quality and apo selection (frozen)

Per included pair: holo pocket mappable to UniProt residues; unique traceable apo-holo residue
mapping; no severe missing residues in pocket; >= 5 positive residues after mapping; identical or
provably equivalent biological assembly; membrane proteins excluded from primary analysis
(UniProt subcellular location + PDB keyword double criterion); interface pockets keep all chains
forming the pocket.

Primary apo per pocket - fixed priority order:
1. biological assembly consistent with holo; 2. complete pocket mapping; 3. higher resolution /
better structure quality; 4. fewer missing residues; 5. PDB ID lexicographic order.
One primary apo per pocket in primary analysis; additional apos reserved for reproducibility
analysis only.

Strict-apo sensitivity subset: target pocket ligand-free AND no other bound molecule capable of
major conformational change (whitelist/blacklist of benign ligands: water only permitted beyond
the target ligand; anything else flags the structure non-strict).

## 5. Static features (frozen)

Per residue: AA identity (20 classes); RSA (Shrake-Rupley, FreeSASA or Biopython fallback,
normalized by theoretical maximum accessibility); residue depth (mean over residue atoms of the
distance to the nearest molecular solvent-accessible surface point, Shrake-Rupley surface);
WCN (sum of 1/d^2 over C-alpha atoms within 10 A, excluding self); secondary structure (DSSP,
3-state: H/E/C); radial position (distance of C-alpha to structure centroid divided by maximum
residue-centroid distance); terminal distance (min(i-1, L-i)/L over observed residues).
Neighborhood context at 10 A C-alpha radius: mean RSA, mean Depth, mean WCN of neighbors
(excluding self). GAM formula per source protocol section 13.1 including te(RSA, Depth).

## 6. Mechanics features (frozen)

ANM: C-alpha coarse-graining; cutoff 10 A; spring constant 1 for contacts; 3N x 3N Hessian;
remove 6 rigid-body zero modes; full pseudo-inverse H+ of the non-rigid block (dense pinv).
Response R_ij = Frobenius norm of the 3x3 block (i,j) of H+ (isotropic unit-force response).
%DFI_i = percentile rank of DFI_i = sum_j R_ij (j != i) within core residues.
Response dispersion: p_ij = R_ij / sum_{k != i} R_k; D_i = -sum_{j != i} p_ij log p_ij / log(N-1),
converted to core-residue percentile.
core_residue_mask (fixed algorithm, uniform over all data): from each sequence terminus, trim
consecutive unresolved or low-contact residues until the first residue with >= 3 contacts within
10 A; bidirectional mapping between original and trimmed indices is stored.
Collinearity pre-check in Pilot: Spearman(%DFI, WCN), Spearman(Dispersion, WCN),
Spearman(%DFI, Dispersion); if |rho| > 0.8 for Dispersion vs WCN, Dispersion is demoted to
secondary feature and the model reduces to Static + %DFI (pre-registered downgrade).

## 7. Models and cross-validation (frozen)

- GAM (primary): pyGAM LogisticGAM, splines s() on RSA, Depth, WCN, Radial, Terminal, %DFI,
  Dispersion; categorical AA and SS; te(RSA, Depth); pair-level weights w_pi = 1/N_p.
- LightGBM (confirmatory): identical folds and feature sets; small hyperparameter grid selected
  inside development folds only.
- CV: MMseqs2 clustering of development sequences at 30% sequence identity, 80% coverage;
  5-fold grouped CV; same cluster/UniProt never spans folds; identical fold assignment for all
  models. Fixed split seed 42.

## 8. Metrics and statistics (frozen)

Primary: per-pair average precision AP_p; Delta AP_p = AP_p^{Static+Mech} - AP_p^{Static};
median Delta AP; cluster bootstrap 95% CI (bootstrap unit = sequence cluster, 1000 resamples,
seeds 1000-1999); win/tie/loss with epsilon = 0.005 (win rate excludes ties); AUROC secondary.
Evaluation within each protein-pocket pair over positive + negative residues (masked excluded).

## 9. Negative controls (frozen)

1. Mechanical profile swap: proteins binned by core length; within-bin random swap of %DFI and
   dispersion profiles (donor interpolated by normalized sequence coordinate if needed); 100
   seeds (0-99); full Static+swapped-Mechanics rerun; real median Delta AP must exceed the null
   95th percentile.
2. Static spatial smoothing: compare Static vs Static-smoothed (10 A neighborhood mean of
   p_static) vs Static+Mechanics.
3. B-factor comparator (X-ray subset, protein-internal percentile): Static / Static+B / Static+M
   / Static+B+M (main-paper analysis, not a Pilot gate).

## 10. Go / No-Go (Pilot gates, frozen)

G1 cluster-bootstrap 95% CI lower bound of median Delta AP_apo > 0;
G2 median Delta AP_apo >= 0.01 (SESOI);
G3 win rate > 50% after excluding ties;
G4 LightGBM median Delta AP > 0 and CI lower bound > 0;
G5 real-mechanics gain > profile-swap null 95th percentile;
G6 Static+Mechanics must also beat Static-smoothed for any strong mechanics claim.
All gates pass -> Go. GAM passes but LightGBM CI crosses 0 or smoothing control weak ->
Conditional Go. Any No-Go condition -> stop expansion, contract claims (source protocol 21.6).

## 11. Sealed test rules (frozen)

CryptoBench test (222 apo structures) sealed before Pilot. Development pool excludes any protein
with sequence identity >= 30% AND alignment coverage >= 80% to any sealed-test sequence
(MMseqs2 prefilter + easy-search). No tuning, feature selection, or threshold choice may use the
sealed test before final model freeze (Step 11).

## 12. Random seed table (frozen)

master seed 42; CV split seed 42; cluster bootstrap seeds 1000-1999; profile-swap seeds 0-99;
LightGBM seed 42; any tie-breaking random draws seed 42.
