#!/usr/bin/env python
"""Step 3 - Compute Static + Mechanics features for every included pair.

Reads manifests/pair_manifest.tsv + residue_mapping.tsv; writes per-pair feature tables
features/pair_features/{pair_id}.csv with one row per evaluable residue:
  static: aa, rsa, depth, wcn, ss, radial, terminal, nb_rsa, nb_depth, nb_wcn
  mechanics: pct_dfi, dispersion (percentile within core residues)
  b_factor_percentile (X-ray only)
QC per pair is appended to audit/mechanics_qc.tsv.
"""
import csv, os, sys
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from features_static import compute_static_features
from features_mechanics import mechanics_features, percentile_rank

PROJ = "/workspace/project"
STRUCT = f"{PROJ}/structures"
FEAT = f"{PROJ}/features/pair_features"
os.makedirs(FEAT, exist_ok=True)

def process_pair(p):
    """Compute features for one pair; returns qc row. Skips if output exists."""
    pid = p["pair_id"]
    out_path = f"{FEAT}/{pid}.csv"
    if os.path.exists(out_path):
        return dict(pair_id=pid, status="SKIP", error="")
    chains = [c for c in p["apo_chains"].split(",") if c]
    cif = f"{STRUCT}/{p['apo_pdb']}.cif"
    try:
        rows, ca = compute_static_features(cif, chains)
    except Exception as e:
        return dict(pair_id=pid, status="FAIL_STATIC", error=str(e)[:100])
    if not rows or ca is None:
        return dict(pair_id=pid, status="FAIL_EMPTY", error="")
    try:
        mech = mechanics_features(ca, [r["n_contacts"] for r in rows])
    except Exception as e:
        return dict(pair_id=pid, status="FAIL_ANM", error=str(e)[:100])
    qc = mech["qc"]
    core = mech["core_mask"]
    pct_dfi = np.full(len(rows), np.nan)
    pdisp = np.full(len(rows), np.nan)
    pct_dfi[core] = percentile_rank(mech["dfi"][core])
    pdisp[core] = percentile_rank(mech["dispersion"][core])
    # B-factor percentile (protein-internal, X-ray or cryo-EM B-factors if present)
    import gemmi
    st = gemmi.read_structure(cif)
    st.setup_entities()
    bmap = {}
    for ch in st[0]:
        if ch.name in chains:
            for r in ch:
                bs = [a.b_iso for a in r]
                bmap[(ch.name, r.seqid.num)] = float(np.mean(bs))
    bvals = np.array([bmap.get((r["chain"], r["resnum"]), np.nan) for r in rows])
    if np.isfinite(bvals).sum() > 10:
        from scipy.stats import rankdata
        ok = np.isfinite(bvals)
        bv = np.full(len(rows), np.nan)
        bv[ok] = 100.0 * (rankdata(bvals[ok], method="average") - 0.5) / ok.sum()
        bfac = bv
    else:
        bfac = bvals
    # labels
    lab = {}
    for rr in RESMAP:
        if rr["pair_id"] == pid:
            lab[(rr["chain"], int(rr["resnum"]))] = rr["label"]
    with open(out_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["pair_id","chain","resnum","resname","label","aa","rsa","depth","wcn",
                    "ss","radial","terminal","nb_rsa","nb_depth","nb_wcn",
                    "pct_dfi","dispersion","b_factor_percentile","in_core"])
        for k, r in enumerate(rows):
            key = (r["chain"], r["resnum"])
            label = lab.get(key, "unmapped")
            w.writerow([pid, r["chain"], r["resnum"], r["resname"], label, r["aa"],
                        f"{r['rsa']:.4f}", f"{r['depth']:.3f}", f"{r['wcn']:.4f}",
                        r["ss"], f"{r['radial']:.4f}", f"{r['terminal']:.4f}",
                        f"{r['nb_rsa']:.4f}", f"{r['nb_depth']:.3f}", f"{r['nb_wcn']:.4f}",
                        f"{pct_dfi[k]:.3f}" if np.isfinite(pct_dfi[k]) else "",
                        f"{pdisp[k]:.3f}" if np.isfinite(pdisp[k]) else "",
                        f"{bfac[k]:.2f}" if np.isfinite(bfac[k]) else "",
                        int(core[k])])
    np.save(f"{FEAT}/{pid}_ca.npy", ca)
    return dict(pair_id=pid, status="OK", n_res=len(rows),
                n_labeled=sum(1 for r in rows if (r["chain"], r["resnum"]) in lab),
                n_zero_modes=qc["n_zero_modes"],
                sym_err=f"{qc['sym_err']:.2e}",
                proj_err=f"{qc['proj_err']:.2e}",
                min_eig=f"{qc['min_nonzero_eig']:.3e}",
                trim_lo=mech["trim_lo"], trim_hi=mech["trim_hi"])

RESMAP = None

def main():
    global RESMAP
    pairs = list(csv.DictReader(open(f"{PROJ}/manifests/pair_manifest.tsv"), delimiter="\t"))
    RESMAP = list(csv.DictReader(open(f"{PROJ}/manifests/residue_mapping.tsv"), delimiter="\t"))
    print(f"pairs: {len(pairs)}", flush=True)
    from multiprocessing import Pool
    qc_rows = []
    with Pool(6) as pool:
        for i, res in enumerate(pool.imap_unordered(process_pair, pairs, chunksize=1)):
            qc_rows.append(res)
            if (i + 1) % 20 == 0:
                n_fail = sum(1 for r in qc_rows if r["status"].startswith("FAIL"))
                print(f"{i+1}/{len(pairs)} done (fail={n_fail})", flush=True)
    with open(f"{PROJ}/audit/mechanics_qc.tsv", "w") as f:
        cols = ["pair_id","status","n_res","n_labeled","n_zero_modes","sym_err","proj_err",
                "min_eig","trim_lo","trim_hi","error"]
        f.write("\t".join(cols) + "\n")
        for r in qc_rows:
            f.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")
    n_fail = sum(1 for r in qc_rows if r["status"].startswith("FAIL"))
    n_ok = sum(1 for r in qc_rows if r["status"] in ("OK", "SKIP"))
    print(f"DONE: ok={n_ok} fail={n_fail}")

if __name__ == "__main__":
    main()
