#!/usr/bin/env python
"""Step 2c - Structure-level pair construction, labels, and audits.

Pipeline (frozen):
1. Candidate pockets from pocket_table.tsv with inclusion OK and n_pos_unp >= 5.
2. Random subsample of POOL_N pockets (seed 42) for structure-level processing.
3. Download apo + holo mmCIF (RCSB); extract observed protein residues.
4. Map observed residues to UniProt indices by pairwise alignment (BLOSUM62, gap -11/-1);
   anchors validated against AHoJ pocket mappings.
5. Ligand heavy-atom count from the holo structure (>= 6 required).
6. Membrane exclusion via UniProt subcellular location annotation.
7. Assembly proxy: observed chain count of the target UniProt must match apo vs holo.
8. Labels: positive = target pocket; MASK = other pockets (same UniProt); negative = rest.
9. CryptoBench homology screen (MMseqs2, 30% id / 80% cov) -> E11 exclusions.
Outputs: pair_manifest.tsv, residue_mapping.tsv, exclusion_log.tsv, dev_sequences.fasta
"""
import csv, json, os, subprocess, sys, time, urllib.request
from collections import defaultdict
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
PROJ = "/workspace/project"
MAN = f"{PROJ}/manifests"
STRUCT = f"{PROJ}/structures"
os.makedirs(STRUCT, exist_ok=True)

POOL_N = 400          # structure-level processing pool (pre-registered subsample)
PILOT_N = 250         # pilot target
SEED = 42
MIN_POS = 5
MIN_HEAVY = 6
CORE_MIN, CORE_MAX = 30, 2000

# ---------------------------------------------------------------- utilities
def http_json(url, retries=3):
    for a in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=60) as r:
                return json.load(r)
        except Exception:
            time.sleep(2 * (a + 1))
    return None

def download_cif(pdb, dst):
    if os.path.exists(dst) and os.path.getsize(dst) > 1000:
        return True
    for a in range(3):
        subprocess.run(["curl", "-sL", "--max-time", "180", "-o", dst,
                        f"https://files.rcsb.org/download/{pdb.upper()}.cif"])
        if os.path.exists(dst) and os.path.getsize(dst) > 1000:
            return True
        time.sleep(2 * (a + 1))
    return False

_unp_cache = {}
def uniprot_seq(acc):
    if acc in _unp_cache:
        return _unp_cache[acc]
    d = http_json(f"https://rest.uniprot.org/uniprotkb/{acc}.json")
    if not d or "sequence" not in d:
        # obsolete/invalid accession
        _unp_cache[acc] = (None, False)
        return None, False
    seq = d["sequence"]["value"]
    membrane = False
    if d:
        for kw in d.get("keywords", []):
            if kw.get("name", "").lower() in {"membrane", "transmembrane"}:
                membrane = True
        for sl in d.get("uniProtKBCrossReferences", []):
            pass
        locs = d.get("comments", [])
        for c in locs:
            if c.get("commentType") == "SUBCELLULAR LOCATION":
                for loc in c.get("subcellularLocations", []):
                    val = loc.get("location", {}).get("value", "").lower()
                    if "membrane" in val:
                        membrane = True
    _unp_cache[acc] = (seq, membrane)
    return seq, membrane

def observed_residues(cif_path, target_unp):
    """Observed protein residues of ALL chains; returns list of dicts with chain, resnum, resname."""
    import gemmi
    st = gemmi.read_structure(cif_path)
    st.setup_entities()
    st.remove_alternative_conformations()
    st.remove_hydrogens()
    out = []
    for ch in st[0]:
        for r in ch:
            info = gemmi.find_tabulated_residue(r.name)
            if info is not None and info.is_amino_acid() and info.is_standard():
                out.append(dict(chain=ch.name, resnum=r.seqid.num, resname=r.name))
    return out

def ligand_heavy_atoms(cif_path, lig, chain, resnum):
    import gemmi
    st = gemmi.read_structure(cif_path)
    st.setup_entities()
    st.remove_hydrogens()
    st.remove_waters()
    n = 0
    for ch in st[0]:
        for r in ch:
            if r.name == lig and (chain is None or ch.name == chain):
                if resnum is not None and r.seqid.num != resnum:
                    continue
                n += len(r)
    return n

_aligner = None
def map_observed_to_uniprot(obs, unp_seq):
    """Map observed residues (list of dicts with resname) to UniProt indices via
    global alignment of the concatenated observed sequence. Returns list of uniprot
    indices (or None) aligned to obs order."""
    from Bio import Align
    from Bio.Align import substitution_matrices
    global _aligner
    if _aligner is None:
        _aligner = Align.PairwiseAligner()
        _aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
        _aligner.open_gap_score = -11
        _aligner.extend_gap_score = -1
        _aligner.mode = "global"
    three2one = {
        "ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E","GLY":"G",
        "HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P","SER":"S",
        "THR":"T","TRP":"W","TYR":"Y","VAL":"V","MSE":"X","SEC":"X","PYL":"X","ASX":"X",
        "GLX":"X","UNK":"X"}
    obs_seq = "".join(three2one.get(r["resname"], "X") for r in obs)
    aln = _aligner.align(unp_seq, obs_seq)[0]
    # map obs index -> uniprot index
    mapping = [None] * len(obs_seq)
    ui, oi = 0, 0
    A = aln.aligned  # blocks of aligned segments
    for (us, ue), (os_, oe) in zip(A[0], A[1]):
        while ui < us:
            ui += 1  # gap in obs
        while oi < os_:
            oi += 1
        while oi < oe and ui < ue:
            mapping[oi] = ui
            ui += 1; oi += 1
    return [mapping[k] if mapping[k] is not None else None for k in range(len(obs_seq))], aln.score

def main():
    table = list(csv.DictReader(open(f"{MAN}/pocket_table.tsv"), delimiter="\t"))
    cands = [r for r in table if r["inclusion"] == "OK" and int(r["n_pos_unp"]) >= MIN_POS]
    print(f"candidates (OK, >= {MIN_POS} pos): {len(cands)}")
    import random
    rng = random.Random(SEED)
    rng.shuffle(cands)
    pool = cands[:POOL_N]
    print(f"processing pool: {len(pool)}")

    pair_rows, resmap_rows, excl_rows = [], [], []
    dev_seqs = {}
    n_done = 0
    for r in pool:
        unp = r["uniprot"]
        pair_id = f"{r['pocket_id']}_{r['apo_pdb']}"
        def excl(code):
            excl_rows.append(dict(pair_id=pair_id, uniprot=unp, pocket_id=r["pocket_id"],
                                  reason=code))
        # UniProt sequence + membrane
        seq, membrane = uniprot_seq(unp)
        if seq is None:
            excl("E13_PARSE_FAIL"); continue
        if membrane:
            excl("E08_MEMBRANE"); continue
        # structures
        apo_pdb, holo_pdb = r["apo_pdb"], r["holo_pdb"]
        apo_cif, holo_cif = f"{STRUCT}/{apo_pdb}.cif", f"{STRUCT}/{holo_pdb}.cif"
        if not (download_cif(apo_pdb, apo_cif) and download_cif(holo_pdb, holo_cif)):
            excl("E13_PARSE_FAIL"); continue
        try:
            apo_obs = observed_residues(apo_cif, unp)
            holo_obs = observed_residues(holo_cif, unp)
        except Exception:
            excl("E13_PARSE_FAIL"); continue
        if not apo_obs or not holo_obs:
            excl("E13_PARSE_FAIL"); continue
        # ligand heavy atoms (from holo)
        holo_chain = r["holo_chain"]
        try:
            heavy = ligand_heavy_atoms(holo_cif, r["ligand"], None, None)
        except Exception:
            heavy = 0
        if heavy < MIN_HEAVY:
            excl("E02_HEAVY_ATOMS"); continue
        # assembly proxy: observed chain count of this UniProt's polymer chains
        # (chains containing the target uniprot residues) - compare apo vs holo
        try:
            apo_map, apo_score = map_observed_to_uniprot(apo_obs, seq)
            holo_map, holo_score = map_observed_to_uniprot(holo_obs, seq)
        except Exception:
            excl("E04_MAPPING"); continue
        # anchor validation: AHoJ holo pocket uniprot indices must be present in holo mapping
        holo_unp_idx = set(map(int, r["unp_idx"].split(",")))
        holo_mapped_idx = {m for m in holo_map if m is not None}
        anchor_ok = len(holo_unp_idx & holo_mapped_idx) / max(len(holo_unp_idx), 1) >= 0.9
        if not anchor_ok:
            excl("E04_MAPPING"); continue
        # observed chains of the target uniprot (chains containing >= 30 mapped residues)
        def chains_of(obs, mp):
            cc = defaultdict(int)
            for o, m in zip(obs, mp):
                if m is not None:
                    cc[o["chain"]] += 1
            return {c for c, n in cc.items() if n >= 30}
        apo_chains, holo_chains = chains_of(apo_obs, apo_map), chains_of(holo_obs, holo_map)
        if len(apo_chains) != len(holo_chains):
            excl("E07_ASSEMBLY"); continue
        # positive residues in apo (observed + mapped)
        pos_apo = [(o, m) for o, m in zip(apo_obs, apo_map)
                   if m is not None and m in holo_unp_idx]
        if len(pos_apo) < MIN_POS:
            excl("E06_FEW_POSITIVES"); continue
        # missing pocket residues in apo?
        if len(pos_apo) < len(holo_unp_idx):
            excl("E05_MISSING_POCKET"); continue
        # other pockets -> MASK
        other_idx = set()
        if r["other_pocket_unp_idx"]:
            other_idx = set(map(int, r["other_pocket_unp_idx"].split(",")))
        mask_apo = [(o, m) for o, m in zip(apo_obs, apo_map)
                    if m is not None and m in other_idx and m not in holo_unp_idx]
        # core size check (observed mapped residues)
        n_core = sum(1 for m in apo_map if m is not None)
        if n_core < CORE_MIN or n_core > CORE_MAX:
            excl("E14_SIZE"); continue
        # strict apo flag: no non-query ligands in apo pocket + no other ligands beyond water
        strict = (r["apo_nonquery_ligs"].strip() == "")
        # sequence for clustering / homology screen
        seq_key = f"{pair_id}|{unp}"
        obs_seq = "".join({o["chain"]: "" for o in apo_obs} and [""])  # placeholder
        chain_seqs = defaultdict(list)
        for o, m in zip(apo_obs, apo_map):
            if m is not None:
                chain_seqs[o["chain"]].append((m, o["resnum"]))
        # longest chain sequence as representative
        best_ch = max(chain_seqs, key=lambda c: len(chain_seqs[c]))
        three2one = {"ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E",
                     "GLY":"G","HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F",
                     "PRO":"P","SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V"}
        rep_seq = "".join(three2one.get(o["resname"], "X") for o in apo_obs
                          if o["chain"] == best_ch)
        dev_seqs[pair_id] = (unp, rep_seq)
        pair_rows.append(dict(
            pair_id=pair_id, uniprot_id=unp, pocket_id=r["pocket_id"],
            apo_pdb=apo_pdb, holo_pdb=holo_pdb,
            apo_chains=",".join(sorted(apo_chains)), holo_chains=",".join(sorted(holo_chains)),
            assembly_id="1", ligand_id=r["ligand"], ligand_heavy_atoms=heavy,
            positive_residue_count=len(pos_apo),
            masked_other_pocket_residue_count=len(mask_apo),
            core_residue_count=n_core,
            membrane_flag="FALSE", missing_pocket_residue_flag="FALSE",
            strict_apo_flag="TRUE" if strict else "FALSE",
            cryptobench_overlap_flag="", sequence_cluster_id="",
            inclusion_status="INCLUDED", exclusion_reason="",
            apo_resolution=r["apo_resolution"], holo_resolution=r.get("holo_resolution", ""),
        ))
        # residue mapping rows
        pos_idx = set(m for _, m in pos_apo)
        mask_idx = set(m for _, m in mask_apo)
        for o, m in zip(apo_obs, apo_map):
            if m is None:
                continue
            label = "positive" if m in pos_idx else ("masked" if m in mask_idx else "negative")
            resmap_rows.append(dict(pair_id=pair_id, uniprot=unp, chain=o["chain"],
                                    resnum=o["resnum"], resname=o["resname"],
                                    uniprot_index=m, label=label))
        n_done += 1
        if n_done % 25 == 0:
            print(f"processed {n_done} included pairs...", flush=True)

    # write manifests
    with open(f"{MAN}/pair_manifest.tsv", "w") as f:
        cols = list(pair_rows[0].keys())
        f.write("\t".join(cols) + "\n")
        for r in pair_rows:
            f.write("\t".join(str(r[c]) for c in cols) + "\n")
    with open(f"{MAN}/residue_mapping.tsv", "w") as f:
        cols = list(resmap_rows[0].keys())
        f.write("\t".join(cols) + "\n")
        for r in resmap_rows:
            f.write("\t".join(str(r[c]) for c in cols) + "\n")
    with open(f"{MAN}/exclusion_log.tsv", "w") as f:
        f.write("pair_id\tuniprot\tpocket_id\treason\n")
        for r in excl_rows:
            f.write(f"{r['pair_id']}\t{r['uniprot']}\t{r['pocket_id']}\t{r['reason']}\n")
    with open(f"{MAN}/dev_sequences.fasta", "w") as f:
        for pid, (unp, s) in dev_seqs.items():
            f.write(f">{pid} {unp}\n{s}\n")
    print(f"DONE: {len(pair_rows)} included pairs, {len(excl_rows)} exclusions")

if __name__ == "__main__":
    main()
