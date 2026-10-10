#!/usr/bin/env python
"""Step 5 - Go/No-Go evaluation (protocol section 21) and pilot figures."""
import json, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

matplotlib.rcParams["font.family"] = ["Liberation Sans", "DejaVu Sans"]
matplotlib.rcParams["svg.fonttype"] = "none"

PROJ = "/workspace/project"
RES = f"{PROJ}/results/pilot"
FIG = f"{PROJ}/figures"
os.makedirs(FIG, exist_ok=True)

def main():
    s = json.load(open(f"{RES}/pilot_summary.json"))
    per_pair = list(csv.DictReader(open(f"{RES}/pilot_per_pair.tsv"), delimiter="\t"))
    import csv
    per_pair = list(csv.DictReader(open(f"{RES}/pilot_per_pair.tsv"), delimiter="\t"))

    gam, lgb, swap, sm = s["gam"], s["lightgbm"], s["swap_null"], s["smoothing"]
    coll = s["collinearity"]

    # ---- gates
    gates = {
        "G1_CI_lower_gt0": gam["lo"] > 0,
        "G2_SESOI_ge_0.01": gam["median"] >= 0.01,
        "G3_winrate_gt_50": gam["win_rate"] > 0.5,
        "G4_LGBM_directional": lgb["median"] > 0 and lgb["lo"] > 0,
        "G5_swap_null": swap["real"] > swap["p95"],
        "G6_beats_smoothing": None,  # decided after comparing numbers
    }
    # G6: mechanics model must beat static-smoothed comparison logic:
    # strong claim requires Static+Mech > Static AND the smoothing gain must not fully
    # explain mechanics gain. We compare median DeltaAP(mech) vs median gain(smoothing).
    g6_strong = gam["median"] > sm["median_gain"]

    decision = "GO" if all(v is True for v in gates.values() if v is not None) and g6_strong \
        else ("CONDITIONAL GO" if gates["G1_CI_lower_gt0"] and gates["G2_SESOI_ge_0.01"]
              and gates["G3_winrate_gt_50"] else "NO-GO")

    report = f"""# Pilot Go / No-Go Report

## Data
- Evaluable protein-pocket pairs: {s['n_pairs']}
- Sequence clusters (MMseqs2 30% id / 80% cov): {s['n_clusters']}
- CV: 5-fold grouped by sequence cluster (seed 42), identical folds for all models.

## Primary result (GAM, apo features -> holo-defined labels)
- median Delta AP = {gam['median']:.4f}, cluster-bootstrap 95% CI [{gam['lo']:.4f}, {gam['hi']:.4f}]
- win/tie/loss (eps=0.005): {gam['wins']}/{gam['ties']}/{gam['losses']}, win rate (excl. ties) = {gam['win_rate']:.3f}

## Confirmatory model (LightGBM, identical folds)
- median Delta AP = {lgb['median']:.4f}, 95% CI [{lgb['lo']:.4f}, {lgb['hi']:.4f}]
- win/tie/loss: {lgb['wins']}/{lgb['ties']}/{lgb['losses']}, win rate = {lgb['win_rate']:.3f}

## Negative controls
- Profile-swap null (100 seeds, GAM): 95th percentile of median Delta AP = {swap['p95']:.4f};
  real = {swap['real']:.4f} -> {'PASS' if gates['G5_swap_null'] else 'FAIL'}
- Static spatial smoothing (10 A neighborhood of static predictions): median AP gain = {sm['median_gain']:.4f}
  vs mechanics median Delta AP = {gam['median']:.4f} -> mechanics {'EXCEEDS' if g6_strong else 'DOES NOT EXCEED'} smoothing gain

## Collinearity audit
- Spearman(%DFI, WCN) = {coll['spearman_dfi_wcn']:.3f}
- Spearman(Dispersion, WCN) = {coll['spearman_disp_wcn']:.3f}
- Spearman(%DFI, Dispersion) = {coll['spearman_dfi_disp']:.3f}
- Dispersion demotion rule (|rho|>0.8 vs WCN): {'TRIGGERED' if abs(coll['spearman_disp_wcn'])>0.8 else 'not triggered'}

## Gate summary
| Gate | Criterion | Value | Verdict |
|------|-----------|-------|---------|
| G1 | CI lower > 0 | {gam['lo']:.4f} | {'PASS' if gates['G1_CI_lower_gt0'] else 'FAIL'} |
| G2 | median Delta AP >= 0.01 | {gam['median']:.4f} | {'PASS' if gates['G2_SESOI_ge_0.01'] else 'FAIL'} |
| G3 | win rate > 50% | {gam['win_rate']:.3f} | {'PASS' if gates['G3_winrate_gt_50'] else 'FAIL'} |
| G4 | LGBM median > 0 and CI lower > 0 | {lgb['median']:.4f} / {lgb['lo']:.4f} | {'PASS' if gates['G4_LGBM_directional'] else 'FAIL'} |
| G5 | real > swap-null 95th pct | {swap['real']:.4f} vs {swap['p95']:.4f} | {'PASS' if gates['G5_swap_null'] else 'FAIL'} |
| G6 | mechanics > smoothing gain | {gam['median']:.4f} vs {sm['median_gain']:.4f} | {'PASS' if g6_strong else 'WEAK (contract claim)'} |

## Decision: **{decision}**
"""
    with open(f"{RES}/go_nogo_report.md", "w") as f:
        f.write(report)
    print(report)

    # ---- figures
    d_gam = [float(r["DeltaAP_gam"]) for r in per_pair if r["DeltaAP_gam"]]
    d_lgb = [float(r["DeltaAP_lgbm"]) for r in per_pair if r["DeltaAP_lgbm"]]
    fig, axes = plt.subplots(1, 2, figsize=(8, 3.2))
    axes[0].hist(d_gam, bins=30, color="#0279EE", alpha=0.75)
    axes[0].axvline(0, color="k", lw=1)
    axes[0].axvline(np.median(d_gam), color="#FF9400", lw=1.5,
                    label=f"median {np.median(d_gam):.3f}")
    axes[0].set_xlabel("$\\Delta$AP (Static+Mechanics $-$ Static)")
    axes[0].set_ylabel("# protein-pocket pairs")
    axes[0].legend(frameon=False)
    axes[0].set_title("GAM")
    axes[1].hist(d_lgb, bins=30, color="#75A025", alpha=0.75)
    axes[1].axvline(0, color="k", lw=1)
    axes[1].axvline(np.median(d_lgb), color="#FF9400", lw=1.5,
                    label=f"median {np.median(d_lgb):.3f}")
    axes[1].set_xlabel("$\\Delta$AP (Static+Mechanics $-$ Static)")
    axes[1].legend(frameon=False)
    axes[1].set_title("LightGBM")
    fig.tight_layout()
    fig.savefig(f"{FIG}/fig3_delta_ap_distribution.png", dpi=200)
    fig.savefig(f"{FIG}/fig3_delta_ap_distribution.svg")
    # swap null
    fig, ax = plt.subplots(figsize=(4, 3.2))
    ax.hist(swap["null_meds"], bins=20, color="#bbbbbb", alpha=0.9,
            label="profile-swap null")
    ax.axvline(swap["real"], color="#0279EE", lw=2, label="real mechanics")
    ax.axvline(swap["p95"], color="#FF9400", ls="--", lw=1.5, label="null 95th pct")
    ax.set_xlabel("median $\\Delta$AP")
    ax.set_ylabel("# null replicates")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(f"{FIG}/fig3_swap_null.png", dpi=200)
    fig.savefig(f"{FIG}/fig3_swap_null.svg")
    print("figures saved")

if __name__ == "__main__":
    main()
