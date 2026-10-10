#!/usr/bin/env python
"""Step 2a - Sample candidate AHoJ-DB entries and download entry zips (frozen seed 42).

Sampling design (pre-registered in protocol lock):
  - universe: Subset 1 Success_jobs.txt entries whose ligand passes the frozen class filter
    and whose query has >= 1 apo pocket (query_summaries.csv);
  - restrict num_apo_pockets <= 100 to avoid mega-family domination of the pilot;
  - uniform random sample of 1500 candidate entries, seed 42;
  - download each entry zip from the AHoJ-DB API with retries; checkpoint to JSONL.
"""
import csv, io, json, os, random, subprocess, sys, time, zipfile
from concurrent.futures import ThreadPoolExecutor, as_completed

sys.path.insert(0, os.path.dirname(__file__))
from ligand_filters import ligand_class as lig_class

PROJ = "/workspace/project"
RAW = f"{PROJ}/raw/ahoj"
CKPT = f"{RAW}/entry_data.jsonl"
SAMPLE_N = 1500
SEED = 42
MAX_APO_POCKETS = 100

z = zipfile.ZipFile(f"{RAW}/subset1.zip")
keys = [l.strip() for l in z.read("Success_jobs.txt").decode().splitlines() if l.strip()]
uniprot_of = {}
apo_pockets_of = {}
with z.open("query_summaries.csv") as f:
    for r in csv.DictReader(io.StringIO(f.read().decode())):
        uniprot_of[r["ahoj_query"]] = r["qUNPs"]
        apo_pockets_of[r["ahoj_query"]] = int(r["num_apo_pockets"])

def query_key_of(entry_key):
    # entry key '101m-A-HEM-155' -> query '101m A HEM 155'
    pdb, chain, lig, idx = entry_key.split("-")
    return f"{pdb} {chain} {lig} {idx}"

cands = []
seen = set()
for k in keys:
    lig = k.split("-")[2]
    if lig_class(lig) != "OK":
        continue
    q = query_key_of(k)
    n_apo = apo_pockets_of.get(q, 0)
    if n_apo <= 0 or n_apo > MAX_APO_POCKETS:
        continue
    if k in seen:
        continue
    seen.add(k)
    cands.append(k)
print(f"candidate entries after ligand + apo filters: {len(cands)}")

rng = random.Random(SEED)
rng.shuffle(cands)
sample = cands[:SAMPLE_N]
print(f"sampling {len(sample)} entries (seed {SEED})")

done = set()
if os.path.exists(CKPT):
    with open(CKPT) as f:
        for line in f:
            try:
                done.add(json.loads(line)["entry_key"])
            except Exception:
                pass
print(f"already collected: {len(done)}")

def fetch(entry_key):
    url = f"https://apoholo.cz/api/db/entry/{entry_key}/download/{entry_key}.zip"
    dst = f"/tmp/{entry_key}.zip"
    for attempt in range(3):
        r = subprocess.run(["curl", "-sL", "--max-time", "120", "-o", dst, url],
                           capture_output=True)
        if os.path.exists(dst) and os.path.getsize(dst) > 500:
            try:
                zf = zipfile.ZipFile(dst)
                data = {n: zf.read(n).decode("utf-8", errors="replace")
                        for n in zf.namelist() if not n.endswith(".gz")}
                inner = {}
                if "residue_mappings.zip" in data:
                    iz = zipfile.ZipFile(io.BytesIO(zf.read("residue_mappings.zip")))
                    for n in iz.namelist():
                        inner[n] = iz.read(n).decode("utf-8", errors="replace")
                data["__residue_mappings__"] = inner
                os.remove(dst)
                return {"entry_key": entry_key, "files": data}
            except Exception as e:
                time.sleep(2 * (attempt + 1))
        else:
            time.sleep(2 * (attempt + 1))
    return {"entry_key": entry_key, "failed": True}

t0 = time.time()
with open(CKPT, "a") as out:
    with ThreadPoolExecutor(max_workers=6) as ex:
        futs = {ex.submit(fetch, k): k for k in sample if k not in done}
        n_ok = n_fail = 0
        for i, fut in enumerate(as_completed(futs)):
            res = fut.result()
            if res.get("failed"):
                n_fail += 1
            else:
                n_ok += 1
                out.write(json.dumps(res) + "\n")
                out.flush()
            if (i + 1) % 100 == 0:
                el = time.time() - t0
                print(f"{i+1}/{len(futs)} ok={n_ok} fail={n_fail} elapsed={el:.0f}s", flush=True)
print(f"DONE ok={n_ok} fail={n_fail}")
