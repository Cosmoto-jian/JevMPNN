#!/usr/bin/env python
"""Step 2b - Build the physical pocket table from collected AHoJ entry zips.

- Pocket identity: within a UniProt, two entries are the SAME physical pocket iff
  Jaccard(UniProt-index sets) >= 0.5 (frozen rule in inclusion_exclusion_rules.md).
- Primary apo per pocket: frozen priority (mapping completeness -> resolution ->
  fewer unobserved -> PDB ID lexicographic).
- Output: manifests/pocket_table.tsv (pre structure-download), used by step2c.
"""
import csv, io, json, os, sys
from collections import defaultdict

PROJ = "/workspace/project"
RAW = f"{PROJ}/raw/ahoj"
CKPT = f"{RAW}/entry_data.jsonl"
OUT = f"{PROJ}/manifests/pocket_table.tsv"
JACCARD_T = 0.5

def parse_mappings(files):
    """Return dict: (struct, poi) -> list of (pdb_chain, pdb_index, uniprot_id, uniprot_index)."""
    out = defaultdict(list)
    inner = files["__residue_mappings__"]
    for name, content in inner.items():
        if not name.endswith(".csv"):
            continue
        rdr = csv.DictReader(io.StringIO(content))
        for r in rdr:
            out[(r["query_struct"], r["POI"])].append(
                (r["PDB_chain"], int(r["PDB_index"]), r["UniProt_ID"], int(r["UniProt_index"])))
    return out

def jaccard(a, b):
    a, b = set(a), set(b)
    if not a or not b:
        return 0.0
    return len(a & b) / len(a | b)

def main():
    entries = []
    with open(CKPT) as f:
        for line in f:
            try:
                rec = json.loads(line)
            except Exception:
                continue
            if rec.get("failed"):
                continue
            entries.append(rec)
    print(f"collected entries: {len(entries)}")

    # per-entry: uniprot, holo pocket uniprot-index set, apo candidates
    pocket_recs = []
    for rec in entries:
        files = rec["files"]
        key = rec["entry_key"]
        pdb, chain, lig, idx = key.split("-")
        try:
            mappings = parse_mappings(files)
        except Exception:
            continue
        poi = f"{chain}_{lig}_{idx}"
        holo_rows = mappings.get((pdb, poi), [])
        if not holo_rows:
            continue
        uniprots = sorted({u for _, _, u, _ in holo_rows})
        if len(uniprots) != 1:
            continue  # multi-UniProt pocket (interface across proteins) - excluded E04
        unp = uniprots[0]
        holo_unp_idx = sorted({ui for _, _, _, ui in holo_rows})
        # apo candidates from apo_filtered_sorted_results.csv
        apo_cands = []
        rdr = csv.DictReader(io.StringIO(files["apo_filtered_sorted_results.csv"]))
        for r in rdr:
            if r.get("apoholo_assignment") != "A":
                continue
            try:
                reso = float(r["resolution"]) if r["resolution"] else 99.0
            except ValueError:
                reso = 99.0
            apo_cands.append(dict(
                apo_pdb=r["structure"], apo_pocket=r["pocket"],
                mapped_pct=float(r["%mapped_obs"] or 0),
                unobserved=int(r["unobserved_rsd"] or 99),
                resolution=reso,
                exp_method=r["exp_method"],
                nonquery_ligs=r["nonquery_ligs"],
                query_ligs=r["query_ligs"],
                n_pocket=int(r["#resis3"] or 0),
            ))
        # holo quality
        qinfo = list(csv.DictReader(io.StringIO(files["query_pocket_info.csv"])))[0]
        pocket_recs.append(dict(
            entry_key=key, uniprot=unp, ligand=lig,
            holo_pdb=pdb, holo_chain=chain,
            holo_unp_idx=holo_unp_idx,
            holo_resolution=float(qinfo["resolution"] or 99),
            holo_exp=qinfo["exp_method"],
            holo_nonquery_ligs=qinfo["nonquery_ligs"],
            apo_cands=apo_cands,
        ))

    print(f"entries with usable mapping: {len(pocket_recs)}")

    # group into physical pockets per UniProt
    by_unp = defaultdict(list)
    for pr in pocket_recs:
        by_unp[pr["uniprot"]].append(pr)
    pockets = []
    for unp, prs in by_unp.items():
        groups = []
        for pr in prs:
            placed = False
            for g in groups:
                if jaccard(pr["holo_unp_idx"], g["unp_idx"]) >= JACCARD_T:
                    g["members"].append(pr)
                    placed = True
                    break
            if not placed:
                groups.append(dict(unp_idx=pr["holo_unp_idx"], members=[pr]))
        for g in groups:
            pockets.append(dict(uniprot=unp, members=g["members"],
                                unp_idx=sorted(set(g["unp_idx"]))))
    print(f"physical pockets: {len(pockets)} across {len(by_unp)} UniProts")

    # assemble pocket table with primary apo selection
    rows = []
    for pi, p in enumerate(pockets):
        unp = p["uniprot"]
        # primary holo: best mapping then resolution
        members = sorted(p["members"], key=lambda m: (m["holo_resolution"], m["entry_key"]))
        holo = members[0]
        # all other-pocket uniprot indices (Jaccard < 0.5 against this pocket by construction)
        other_idx = set()
        for q in pockets:
            if q is p or q["uniprot"] != unp:
                continue
            other_idx |= set(q["unp_idx"])
        other_idx -= set(p["unp_idx"])
        # primary apo: frozen priority
        cands = [c for m in p["members"] for c in m["apo_cands"]]
        if not cands:
            rows.append(dict(uniprot=unp, pocket_id=f"P{pi:05d}", ligand=holo["ligand"],
                             holo_pdb=holo["holo_pdb"], holo_chain=holo["holo_chain"],
                             n_pos_unp=len(p["unp_idx"]), n_apo_cands=0,
                             apo_pdb="", inclusion="E10_NO_APO", other_pocket_unp_idx=""))
            continue
        cands.sort(key=lambda c: (-c["mapped_pct"], c["unobserved"], c["resolution"],
                                  -c["n_pocket"], c["apo_pdb"]))
        best = cands[0]
        rows.append(dict(uniprot=unp, pocket_id=f"P{pi:05d}", ligand=holo["ligand"],
                         holo_pdb=holo["holo_pdb"], holo_chain=holo["holo_chain"],
                         holo_entry=holo["entry_key"],
                         n_pos_unp=len(p["unp_idx"]), n_apo_cands=len(cands),
                         apo_pdb=best["apo_pdb"], apo_pocket=best["apo_pocket"],
                         apo_mapped_pct=best["mapped_pct"], apo_unobserved=best["unobserved"],
                         apo_resolution=best["resolution"], apo_exp=best["exp_method"],
                         apo_nonquery_ligs=best["nonquery_ligs"],
                         other_pocket_unp_idx=",".join(map(str, sorted(other_idx))),
                         unp_idx=",".join(map(str, p["unp_idx"])),
                         inclusion="OK"))
    os.makedirs(f"{PROJ}/manifests", exist_ok=True)
    cols = ["uniprot","pocket_id","ligand","holo_pdb","holo_chain","holo_entry","n_pos_unp",
            "n_apo_cands","apo_pdb","apo_pocket","apo_mapped_pct","apo_unobserved",
            "apo_resolution","apo_exp","apo_nonquery_ligs","unp_idx","other_pocket_unp_idx","inclusion"]
    with open(OUT, "w") as f:
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")
    n_ok = sum(1 for r in rows if r["inclusion"] == "OK")
    print(f"pocket table written: {len(rows)} pockets, {n_ok} with apo candidate")

if __name__ == "__main__":
    main()
