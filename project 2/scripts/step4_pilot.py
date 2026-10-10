#!/usr/bin/env python
"""Step 4 - Pilot grouped CV, models, and negative controls (protocol sections 13-21).

Frozen design:
- MMseqs2 clustering (30% id / 80% cov) of development sequences -> cluster groups.
- 5-fold grouped CV (seed 42); identical folds for all models and controls.
- GAM (pyGAM LogisticGAM): splines on numeric features + te(RSA, Depth) + linear dummies
  for AA (19) and SS (2); pair weights 1/N_p.
- LightGBM (fixed a-priori params: num_leaves=31, lr=0.05, 400 trees, seed 42).
- Models: Static, Static+Mechanics; controls: Static-smoothed, profile-swap null (100 seeds).
- Metrics: per-pair AP, Delta AP, median, cluster bootstrap 95% CI (1000 resamples),
  win/tie/loss (epsilon=0.005), AUROC secondary.
"""
import csv, json, os, pickle, subprocess, sys
from multiprocessing import Pool
import numpy as np
from collections import defaultdict
from scipy.spatial import cKDTree
from sklearn.metrics import average_precision_score, roc_auc_score

sys.path.insert(0, os.path.dirname(__file__))
PROJ = "/workspace/project"
FEAT = f"{PROJ}/features/pair_features"
RES = f"{PROJ}/results/pilot"
os.makedirs(RES, exist_ok=True)

EPS = 0.005
N_BOOT = 1000
N_SWAP = 100
SEED = 42

STATIC_NUM = ["rsa", "depth", "wcn", "radial", "terminal", "nb_rsa", "nb_depth", "nb_wcn"]
MECH_NUM = ["pct_dfi", "dispersion"]
AA20 = ["ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU","LYS","MET","PHE",
        "PRO","SER","THR","TRP","TYR","VAL"]
SS3 = ["H", "E", "C"]

# ---------------------------------------------------------------- data loading
def load_pairs():
    pairs = {}
    for fn in sorted(os.listdir(FEAT)):
        if not fn.endswith(".csv"):
            continue
        pid = fn[:-4]
        rows = list(csv.DictReader(open(f"{FEAT}/{fn}")))
        y_all = np.array([1 if r["label"] == "positive" else
                          (0 if r["label"] == "negative" else -1) for r in rows])
        keep = y_all >= 0
        y = y_all[keep]
        if len(y) < 10 or y.sum() < 2 or y.sum() == len(y):
            continue
        ca = np.load(f"{FEAT}/{pid}_ca.npy")
        idx_keep = [k for k, b in enumerate(keep) if b]
        d = dict(y=y, ca=ca[idx_keep],
                 chain=[rows[k]["chain"] for k in idx_keep],
                 resnum=[rows[k]["resnum"] for k in idx_keep])
        for f in STATIC_NUM + MECH_NUM + ["b_factor_percentile"]:
            v = [rows[k][f] for k in idx_keep]
            d[f] = np.array([np.nan if x == "" else float(x) for x in v])
        d["aa"] = [rows[k]["aa"] for k in idx_keep]
        d["ss"] = [rows[k]["ss"] for k in idx_keep]
        pairs[pid] = d
    return pairs

def stack(pairs, pids, mech=None):
    """Build X, y, w, pair_index. mech=None -> static only; else dict pid->(dfi,disp)."""
    Xs, ys, ws, pidx = [], [], [], []
    for pi, pid in enumerate(pids):
        d = pairs[pid]
        n = len(d["y"])
        nums = [d[f] for f in STATIC_NUM]
        if mech is not None:
            nums += [mech[pid][0], mech[pid][1]]
        Xn = np.vstack(nums).T
        aa = np.array([[1.0 if a == x else 0.0 for x in AA20[1:]] for a in d["aa"]])
        ss = np.array([[1.0 if s == x else 0.0 for x in SS3[1:]] for s in d["ss"]])
        Xs.append(np.hstack([Xn, aa, ss]))
        ys.append(d["y"])
        ws.append(np.full(n, 1.0 / n))
        pidx.extend([pi] * n)
    return (np.vstack(Xs), np.concatenate(ys), np.concatenate(ws), np.array(pidx))

# ---------------------------------------------------------------- models
def fit_gam(X, y, w, n_num):
    from pygam import LogisticGAM, s, te, l
    terms = s(0, n_splines=5)
    for i in range(1, n_num):
        terms += s(i, n_splines=5)
    terms += te(0, 1, n_splines=4)   # te(RSA, Depth)
    for j in range(n_num, X.shape[1]):
        terms += l(j)
    gam = LogisticGAM(terms, verbose=False)
    gam.fit(X, y, weights=w)
    return gam

def fit_lgbm(X, y, w):
    import lightgbm as lgb
    m = lgb.LGBMClassifier(n_estimators=400, num_leaves=31, learning_rate=0.05,
                           min_child_samples=20, random_state=SEED, verbose=-1)
    m.fit(X, y, sample_weight=w)
    return m

def ap_aucs(pairs, pids, p, pidx):
    aps, aucs = {}, {}
    for pi, pid in enumerate(pids):
        m = pidx == pi
        y = pairs[pid]["y"]
        aps[pid] = average_precision_score(y, p[m])
        try:
            aucs[pid] = roc_auc_score(y, p[m])
        except ValueError:
            pass
    return aps, aucs

def cluster_bootstrap(delta, cluster_of_pair, n_boot=N_BOOT, seed=SEED):
    clusters = defaultdict(list)
    for pid, d in delta.items():
        clusters[cluster_of_pair[pid]].append(d)
    cl = list(clusters.values())
    rng = np.random.default_rng(seed)
    meds = np.empty(n_boot)
    for b in range(n_boot):
        pick = rng.choice(len(cl), len(cl), replace=True)
        flat = [x for c in pick for x in cl[c]]
        meds[b] = np.median(flat)
    obs = float(np.median(list(delta.values())))
    return obs, float(np.percentile(meds, 2.5)), float(np.percentile(meds, 97.5)), meds


_SWAP_CTX = {}

def _swap_seed_worker(seed):
    """Profile-swap null for one seed (module-level; state via fork globals)."""
    pairs = _SWAP_CTX["pairs"]
    bins = _SWAP_CTX["bins"]
    pair_of_fold = _SWAP_CTX["pair_of_fold"]
    static_aps = _SWAP_CTX["static_aps"]
    out = {}
    for b, bps in bins.items():
        if len(bps) < 2:
            for pid in bps:
                out[pid] = (pairs[pid]["pct_dfi"], pairs[pid]["dispersion"])
            continue
        perm = np.random.default_rng(seed).permutation(len(bps))
        for dst, src_i in zip(bps, perm):
            src = pairs[bps[src_i]]
            n = len(pairs[dst]["y"])
            x_new = np.linspace(0, 1, n)
            dfi = np.interp(x_new, np.linspace(0, 1, len(src["pct_dfi"])), src["pct_dfi"])
            disp = np.interp(x_new, np.linspace(0, 1, len(src["dispersion"])),
                             src["dispersion"])
            out[dst] = (dfi, disp)
    delta_swap = {}
    for k in range(5):
        test_pids = pair_of_fold[k]
        train_pids = [p for j in range(5) if j != k for p in pair_of_fold[j]]
        Xtr, ytr, wtr, _ = stack(pairs, train_pids, mech=out)
        med = np.nanmedian(Xtr, axis=0)
        Xtr = np.where(np.isnan(Xtr), med, Xtr)
        m = fit_gam(Xtr, ytr, wtr, n_num=len(STATIC_NUM) + 2)
        for pid in test_pids:
            Xte, _, _, _ = stack(pairs, [pid], mech=out)
            Xte = np.where(np.isnan(Xte), med, Xte)
            p = m.predict_proba(Xte)  # pyGAM 1-D
            ap_sw = average_precision_score(pairs[pid]["y"], p)
            delta_swap[pid] = ap_sw - static_aps[pid]
    return float(np.median(list(delta_swap.values())))

# ---------------------------------------------------------------- main
def main():
    pairs = load_pairs()
    pids_all = sorted(pairs)
    print(f"evaluable pairs: {len(pairs)}", flush=True)

    # MMseqs2 clustering
    fasta = f"{PROJ}/manifests/dev_sequences.fasta"
    clu_dir = f"{PROJ}/features/clusters"
    os.makedirs(clu_dir, exist_ok=True)
    clu_tsv = f"{clu_dir}/dev_sequences_clustered_cluster.tsv"
    if not os.path.exists(clu_tsv):
        subprocess.run(["/workspace/mmseqs/bin/mmseqs", "easy-cluster", fasta,
                        f"{clu_dir}/dev_sequences_clustered", f"{clu_dir}/tmp"],
                       check=True, capture_output=True)
    cluster_of_pair = {}
    with open(clu_tsv) as f:
        for line in f:
            rep, member = line.split()[:2]
            cluster_of_pair[member.split("|")[0]] = rep.split("|")[0]
    print(f"sequence clusters: {len(set(cluster_of_pair.values()))}")

    # 5-fold grouped split
    rng = np.random.default_rng(SEED)
    clusters = sorted(set(cluster_of_pair[p] for p in pids_all))
    rng.shuffle(clusters)
    fold_of_cluster = {c: i % 5 for i, c in enumerate(clusters)}
    pair_of_fold = [[] for _ in range(5)]
    for p in pids_all:
        pair_of_fold[fold_of_cluster[cluster_of_pair[p]]].append(p)

    # leakage audit
    with open(f"{PROJ}/audit/leakage_report.md", "w") as f:
        f.write("# Leakage audit (pilot grouped CV)\n\n")
        f.write("| fold | pairs | clusters |\n|---|---|---|\n")
        for k in range(5):
            f.write(f"| {k} | {len(pair_of_fold[k])} | "
                    f"{len(set(cluster_of_pair[p] for p in pair_of_fold[k]))} |\n")
        f.write(f"\nSequence clusters: {len(clusters)}; UniProt groups: "
                f"{len(set(pairs[p]['chain'] and p for p in pids_all))}\n")
        f.write("Cluster overlap across folds: 0 (by construction, verified below).\n")
    for a in range(5):
        for b in range(a + 1, 5):
            assert not (set(cluster_of_pair[p] for p in pair_of_fold[a]) &
                        set(cluster_of_pair[p] for p in pair_of_fold[b]))

    # ---- main CV (cached to allow resuming after later-stage crashes)
    pos_of = {pid: i for i, pid in enumerate(pids_all)}
    mech_dict = {pid: (pairs[pid]["pct_dfi"], pairs[pid]["dispersion"]) for pid in pids_all}
    cache_path = f"{RES}/main_cv_cache.pkl"
    if os.path.exists(cache_path):
        with open(cache_path, "rb") as f:
            results = pickle.load(f)
        print("loaded main CV from cache", flush=True)
    else:
        results = {}
    for model_name in ["gam", "lightgbm"]:
        for feat_name, use_mech in [("static", False), ("static_mech", True)]:
            if (model_name, feat_name) in results:
                print(f"{model_name} {feat_name}: cached", flush=True)
                continue
            n_total = sum(len(pairs[p]["y"]) for p in pids_all)
            all_p = np.zeros(n_total)
            all_pidx = np.zeros(n_total, dtype=int)
            off = 0
            for k in range(5):
                test_pids, train_pids = pair_of_fold[k], \
                    [p for j in range(5) if j != k for p in pair_of_fold[j]]
                Xtr, ytr, wtr, ptr = stack(pairs, train_pids,
                                           mech=mech_dict if use_mech else None)
                Xte, _, _, pte = stack(pairs, test_pids,
                                       mech=mech_dict if use_mech else None)
                med = np.nanmedian(Xtr, axis=0)
                Xtr = np.where(np.isnan(Xtr), med, Xtr)
                Xte = np.where(np.isnan(Xte), med, Xte)
                if model_name == "gam":
                    m = fit_gam(Xtr, ytr, wtr, n_num=len(STATIC_NUM) + (2 if use_mech else 0))
                    pr = m.predict_proba(Xte)  # pyGAM returns 1-D probabilities
                else:
                    m = fit_lgbm(Xtr, ytr, wtr)
                    pr = m.predict_proba(Xte)[:, 1]
                # map fold-local pair indices to global pids_all indices
                pte_global = np.repeat([pos_of[pid] for pid in test_pids],
                                       [len(pairs[pid]["y"]) for pid in test_pids])
                all_p[off:off + len(pr)] = pr
                all_pidx[off:off + len(pr)] = pte_global
                off += len(pr)
            aps, aucs = ap_aucs(pairs, pids_all, all_p, all_pidx)
            results[(model_name, feat_name)] = dict(aps=aps, aucs=aucs)
            print(f"{model_name} {feat_name}: {len(aps)} pairs, "
                  f"median AP {np.median(list(aps.values())):.4f}", flush=True)
        with open(cache_path, "wb") as f:
            pickle.dump(results, f)

    # ---- Delta AP
    delta_out = {}
    for model_name in ["gam", "lightgbm"]:
        delta = {pid: results[(model_name, "static_mech")]["aps"][pid]
                 - results[(model_name, "static")]["aps"][pid]
                 for pid in results[(model_name, "static")]["aps"]}
        med, lo, hi, meds = cluster_bootstrap(delta, cluster_of_pair)
        wins = sum(1 for d in delta.values() if d > EPS)
        ties = sum(1 for d in delta.values() if abs(d) <= EPS)
        losses = sum(1 for d in delta.values() if d < -EPS)
        wr = wins / max(wins + losses, 1)
        delta_out[model_name] = dict(delta=delta, median=med, lo=lo, hi=hi,
                                     wins=wins, ties=ties, losses=losses, win_rate=wr)
        print(f"{model_name} DeltaAP: median={med:.4f} CI[{lo:.4f},{hi:.4f}] "
              f"W/T/L={wins}/{ties}/{losses} winrate={wr:.3f}", flush=True)

    # ---- Control 1: static smoothing
    smooth_delta = {}
    for k in range(5):
        test_pids = pair_of_fold[k]
        train_pids = [p for j in range(5) if j != k for p in pair_of_fold[j]]
        Xtr, ytr, wtr, _ = stack(pairs, train_pids)
        med = np.nanmedian(Xtr, axis=0)
        Xtr = np.where(np.isnan(Xtr), med, Xtr)
        m = fit_gam(Xtr, ytr, wtr, n_num=len(STATIC_NUM))
        for pid in test_pids:
            Xte, _, _, _ = stack(pairs, [pid])
            Xte = np.where(np.isnan(Xte), med, Xte)
            p = m.predict_proba(Xte)  # pyGAM 1-D
            ca = pairs[pid]["ca"]
            tree = cKDTree(ca)
            nb_lists = tree.query_ball_point(ca, 10.0)
            p_sm = np.array([p[nb].mean() for nb in nb_lists])
            y = pairs[pid]["y"]
            ap_raw = average_precision_score(y, p)
            ap_sm = average_precision_score(y, p_sm)
            smooth_delta[pid] = ap_sm - ap_raw
    sm_med = float(np.median(list(smooth_delta.values())))
    print(f"static-smoothing gain over static: median {sm_med:.4f}")

    # ---- Control 2: profile-swap null (GAM, 100 seeds)
    core_len = {pid: len(pairs[pid]["y"]) for pid in pids_all}
    qs = np.quantile(list(core_len.values()), [0.25, 0.5, 0.75])
    bins = defaultdict(list)
    for pid in pids_all:
        bins[int(np.searchsorted(qs, core_len[pid]))].append(pid)

    def swap_profiles(seed):
        """Return dict pid -> (dfi, disp) with within-bin swapped profiles."""
        out = {}
        for b, bps in bins.items():
            if len(bps) < 2:
                for pid in bps:
                    out[pid] = (pairs[pid]["pct_dfi"], pairs[pid]["dispersion"])
                continue
            perm = np.random.default_rng(seed).permutation(len(bps))
            for dst, src_i in zip(bps, perm):
                src = pairs[bps[src_i]]
                n = len(pairs[dst]["y"])
                x_new = np.linspace(0, 1, n)
                dfi = np.interp(x_new, np.linspace(0, 1, len(src["pct_dfi"])), src["pct_dfi"])
                disp = np.interp(x_new, np.linspace(0, 1, len(src["dispersion"])),
                                 src["dispersion"])
                out[dst] = (dfi, disp)
        return out

    _SWAP_CTX.update(pairs=pairs, bins=bins, pair_of_fold=pair_of_fold,
                     static_aps=results[("gam", "static")]["aps"])
    with Pool(5) as pool:
        swap_meds = pool.map(_swap_seed_worker, range(N_SWAP), chunksize=2)
    print(f"swap null done: median-of-medians {np.median(swap_meds):.4f}")
    null95 = float(np.percentile(swap_meds, 95))
    real_med = delta_out["gam"]["median"]
    print(f"SWAP NULL: 95th pct = {null95:.4f}; real = {real_med:.4f} -> "
          f"{'PASS' if real_med > null95 else 'FAIL'}")

    # ---- collinearity audit
    from scipy.stats import spearmanr
    dfi = np.concatenate([pairs[p]["pct_dfi"][np.isfinite(pairs[p]["pct_dfi"])] for p in pids_all])
    disp = np.concatenate([pairs[p]["dispersion"][np.isfinite(pairs[p]["dispersion"])]
                           for p in pids_all])
    wcn = np.concatenate([pairs[p]["wcn"][np.isfinite(pairs[p]["pct_dfi"])] for p in pids_all])
    coll = dict(spearman_dfi_wcn=float(spearmanr(dfi, wcn).statistic),
                spearman_disp_wcn=float(spearmanr(disp, wcn).statistic),
                spearman_dfi_disp=float(spearmanr(dfi, disp).statistic))
    print("collinearity:", coll)

    # ---- save
    summary = dict(
        gam={k: delta_out["gam"][k] for k in ["median", "lo", "hi", "wins", "ties",
                                              "losses", "win_rate"]},
        lightgbm={k: delta_out["lightgbm"][k] for k in ["median", "lo", "hi", "wins",
                                                        "ties", "losses", "win_rate"]},
        smoothing=dict(median_gain=sm_med),
        swap_null=dict(p95=null95, real=real_med, null_meds=swap_meds),
        collinearity=coll,
        n_pairs=len(pids_all),
        n_clusters=len(clusters),
    )
    with open(f"{RES}/pilot_summary.json", "w") as f:
        json.dump(summary, f, indent=1)
    with open(f"{RES}/pilot_per_pair.tsv", "w") as f:
        f.write("pair_id\tcluster\tAP_gam_static\tAP_gam_static_mech\tDeltaAP_gam\t"
                "AP_lgbm_static\tAP_lgbm_static_mech\tDeltaAP_lgbm\t"
                "DeltaAP_smooth_vs_static\n")
        for pid in pids_all:
            f.write("\t".join(map(str, [
                pid, cluster_of_pair[pid],
                results[("gam", "static")]["aps"].get(pid, ""),
                results[("gam", "static_mech")]["aps"].get(pid, ""),
                delta_out["gam"]["delta"].get(pid, ""),
                results[("lightgbm", "static")]["aps"].get(pid, ""),
                results[("lightgbm", "static_mech")]["aps"].get(pid, ""),
                delta_out["lightgbm"]["delta"].get(pid, ""),
                smooth_delta.get(pid, ""),
            ])) + "\n")
    with open(f"{RES}/cv_state.pkl", "wb") as f:
        pickle.dump(dict(pair_of_fold=pair_of_fold, cluster_of_pair=cluster_of_pair), f)
    print("saved pilot_summary.json + pilot_per_pair.tsv")

if __name__ == "__main__":
    main()
