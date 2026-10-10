#!/usr/bin/env python
"""Step 1 - Seal the CryptoBench official test split (protocol section 4.2).

Outputs (all under project/protocol/ and project/raw/cryptobench_sealed/):
  cryptobench_test_sealed.tsv   one row per test apo structure
  test_sequences.fasta          sequences of the 222 test apo structures
  cryptobench_sealed_checksums.sha256
Read-only afterwards. Download date and versions recorded in environment_lock.yml.
"""
import json, hashlib, os, subprocess, sys, time

PROJ = "/workspace/project"
RAW = f"{PROJ}/raw/cryptobench_sealed"
PROTO = f"{PROJ}/protocol"
os.makedirs(RAW, exist_ok=True)

# 1. Fetch official files from OSF if not already present
for fname, url in [("dataset.json", "https://osf.io/download/ta2ju/"),
                   ("folds.json", "https://osf.io/download/5s93p/"),
                   ("README.md", "https://osf.io/download/x9dce/")]:
    dst = f"{RAW}/{fname}"
    if not (os.path.exists(dst) and os.path.getsize(dst) > 0):
        subprocess.run(["curl", "-sL", "-o", dst, url], check=True)
        print(f"downloaded {fname}: {os.path.getsize(dst)} bytes")

# auxiliary noncryptic pockets (used only at CryptoBench stage for MASK labels)
aux_dir = f"{RAW}/auxiliary-data/non_cryptic_pockets"
os.makedirs(aux_dir, exist_ok=True)
aux_dst = f"{aux_dir}/noncryptic_pockets.json"
if not os.path.exists(aux_dst):
    # direct OSF download link (verified 2026-10-03: auxiliary-data/non-cryptic-pockets/)
    aux_url = "https://osf.io/download/37ztm/"
    subprocess.run(["curl", "-sL", "-o", aux_dst, aux_url], check=True)
    print(f"downloaded noncryptic_pockets.json: {os.path.getsize(aux_dst)} bytes")

dataset = json.load(open(f"{RAW}/dataset.json"))
folds = json.load(open(f"{RAW}/folds.json"))
test_ids = sorted(folds["test"])
print(f"dataset apo structures: {len(dataset)}; official test: {len(test_ids)}")
assert len(test_ids) == 222, "official test split must be 222"

# 2. Download the 222 test apo structures (mmCIF) and extract sequences
struct_dir = f"{RAW}/structures"
os.makedirs(struct_dir, exist_ok=True)
def cif_path(pdb): return f"{struct_dir}/{pdb}.cif"

missing = [p for p in test_ids if not os.path.exists(cif_path(p))]
print(f"structures to download: {len(missing)}")
for i, pdb in enumerate(missing):
    url = f"https://files.rcsb.org/download/{pdb.upper()}.cif"
    for attempt in range(3):
        r = subprocess.run(["curl", "-sL", "-o", cif_path(pdb), url], capture_output=True)
        if os.path.exists(cif_path(pdb)) and os.path.getsize(cif_path(pdb)) > 1000:
            break
        time.sleep(2 * (attempt + 1))
    if (i + 1) % 50 == 0:
        print(f"  downloaded {i+1}/{len(missing)}")

# 3. Extract sequences + UniProt mapping with gemmi; build sealed TSV
import gemmi
rows, fasta = [], []
failures = []
for pdb in test_ids:
    try:
        st = gemmi.read_structure(cif_path(pdb))
        st.setup_entities()
        # polymer sequences per chain (first seqres entity model)
        seqs = {}
        for ch in st[0]:
            if ch.get_polymer():
                seq = gemmi.one_letter_code([r.name for r in ch])
                if seq.strip("X"):
                    seqs[ch.name] = seq
        entries = dataset[pdb]
        uniprots = sorted({e["uniprot_id"] for e in entries})
        ligands = sorted({e["ligand"] for e in entries})
        n_pockets = len(entries)
        rows.append(dict(apo_pdb=pdb, n_pocket_pairs=n_pockets,
                         uniprot_ids=",".join(uniprots), ligands=",".join(ligands),
                         n_chains_with_sequence=len(seqs)))
        for ch, seq in seqs.items():
            fasta.append(f">{pdb}_{ch}\n{seq}\n")
    except Exception as e:
        failures.append((pdb, str(e)))

print(f"sealed rows: {len(rows)}; failures: {len(failures)} {failures[:5]}")

with open(f"{PROTO}/cryptobench_test_sealed.tsv", "w") as f:
    cols = ["apo_pdb", "n_pocket_pairs", "uniprot_ids", "ligands", "n_chains_with_sequence"]
    f.write("\t".join(cols) + "\n")
    for r in rows:
        f.write("\t".join(str(r[c]) for c in cols) + "\n")

with open(f"{RAW}/test_sequences.fasta", "w") as f:
    f.writelines(fasta)

# 4. Checksums of everything sealing-related
to_hash = [f"{PROTO}/cryptobench_test_sealed.tsv", f"{RAW}/dataset.json",
           f"{RAW}/folds.json", f"{RAW}/test_sequences.fasta", aux_dst]
with open(f"{PROTO}/cryptobench_sealed_checksums.sha256", "w") as f:
    for p in to_hash:
        h = hashlib.sha256(open(p, "rb").read()).hexdigest()
        f.write(f"{h}  {p}\n")
print(open(f"{PROTO}/cryptobench_sealed_checksums.sha256").read())
print("SEALED. This directory is read-only from now on.")
