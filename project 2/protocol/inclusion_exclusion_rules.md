# Inclusion / Exclusion Rules (Frozen at Protocol Lock)

Every excluded pair receives one reason code below in exclusion_log.tsv. Rules are applied in
this exact order.

## Reason codes

| Code | Rule |
|------|------|
| E01_LIGAND_CLASS | ligand is water, single-atom ion, nucleic acid, peptide/protein, or a crystallization additive / buffer molecule (PEG, GOL, EDO, MES, ACT, SO4, PO4, CIT, TRS, DMS, MPD, and similar; fixed additive list maintained in scripts/ligand_filters.py) |
| E02_HEAVY_ATOMS | ligand heavy-atom count < 6 |
| E03_NOT_BIOREL | AHoJ/BioLiP2 biological-relevance flag absent or negative |
| E04_MAPPING | holo pocket cannot be uniquely mapped to UniProt / apo residues |
| E05_MISSING_POCKET | pocket residues severely unobserved in apo (any positive residue unobserved) |
| E06_FEW_POSITIVES | < 5 positive residues after mapping |
| E07_ASSEMBLY | apo/holo biological assembly state mismatch that cannot be resolved to an equivalent assembly |
| E08_MEMBRANE | membrane protein (UniProt subcellular location annotation "Cell membrane"/"Membrane" OR PDB structure keywords/COMPOUND indicate membrane) - primary analysis only |
| E09_INTERFACE_INCOMPLETE | interface pocket without all chains forming the pocket present in the apo structure |
| E10_NO_APO | no apo structure available for the pocket |
| E11_CRYPTOBENCH_HOMOLOG | sequence identity >= 30% AND coverage >= 80% to a CryptoBench sealed-test sequence |
| E12_DUPLICATE | duplicate PDB/chain/pocket identity (keeps first by fixed priority order) |
| E13_PARSE_FAIL | structure file or AHoJ entry unparseable / API failure after retries |
| E14_SIZE | observed core residues < 30 or > 2000 (numerical limits of dense ANM pipeline, fixed a priori) |

## Primary-apo selection (per pocket, frozen priority)

1. biological assembly consistent with holo; 2. complete pocket mapping; 3. higher resolution;
4. fewer missing residues; 5. PDB ID lexicographic order. Non-selected apos are flagged
`reproducibility_only` and excluded from primary per-pocket evaluation.

## Pocket identity (frozen)

Pockets are expressed in UniProt residue indices (AHoJ residue_mappings provide the PDB ->
UniProt index mapping). Two AHoJ entries belong to the SAME physical pocket iff
Jaccard(uniport_index_set_A, uniprot_index_set_B) >= 0.5. Entries of the same UniProt with
Jaccard < 0.5 against every existing pocket define a NEW pocket (these become MASK sources for
the target pocket). All holos of one physical pocket pool their positive labels (union);
different pockets never pool.

## Strict-apo flag

strict_apo_flag = TRUE iff the apo structure contains no ligand at the target pocket and no other
bound molecule outside the benign whitelist (water only). Used for the sensitivity subset only.

## Audit questions (source protocol section 7)

All seven audit questions are answered in audit/leakage_report.md before any model training.
