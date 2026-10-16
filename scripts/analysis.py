import numpy as np, pandas as pd, matplotlib.pyplot as plt
from scipy import stats
from sklearn.metrics import roc_auc_score, roc_curve
from pathlib import Path

CSV_PATH = "data/processed/final_dataset.csv"    
COL_PRED_NM = "pred_ic50_nM"
COL_EXP_LOG = "exp_log10_ic50"
COL_EVIDENCE = "evidence"

N_BOOT = 2000 
RNG = np.random.default_rng(0)

df = pd.read_csv(CSV_PATH)
df["pred_log10_ic50"] = np.log10(df[COL_PRED_NM])
print(f"Peptides with predictions: {len(df):,}")
print(df[COL_EVIDENCE].value_counts().to_string(), "\n")

acc = df.dropna(subset=[COL_EXP_LOG]).copy()
x = acc[COL_EXP_LOG].to_numpy() 
y = acc["pred_log10_ic50"].to_numpy() 
n = len(acc)
print(f"Peptides compared (n): {n:,}")

rho, _ = stats.spearmanr(x, y)
boot = []
for _ in range(N_BOOT):
    i = RNG.integers(0, n, n)
    boot.append(stats.spearmanr(x[i], y[i])[0])
lo, hi = np.percentile(boot, [2.5, 97.5])
print(f"Spearman rho = {rho:.3f}  (95% CI {lo:.3f} to {hi:.3f})")

mae = np.mean(np.abs(y - x))
print(f"MAE = {mae:.3f} log10 units  ->  about {10**mae:.1f}-fold")

slope, intercept = stats.linregress(x, y)[:2]
print(f"Fit: y = {slope:.2f}x + {intercept:.2f}; crosses y = x at "
      f"{10 ** (intercept / (1 - slope)):.0f} nM\n")

for thr_nm in (500, 50):
    is_binder = (x < np.log10(thr_nm)).astype(int)
    auc = roc_auc_score(is_binder, -y)
    print(f"AUC at {thr_nm} nM = {auc:.3f}  "
          f"({is_binder.sum():,} binders vs {(1 - is_binder).sum():,} weaker)")

t = df.loc[df[COL_EVIDENCE] == "T-cell+", "pred_log10_ic50"].to_numpy()
b = df.loc[df[COL_EVIDENCE] == "binding+", "pred_log10_ic50"].to_numpy()
print("\nMedian predicted IC50 (nM):")
for name, grp in df.groupby(COL_EVIDENCE):
    print(f"  {name:9s} n = {len(grp):,}  median = {grp[COL_PRED_NM].median():.0f} nM")

res = stats.mannwhitneyu(t, b, alternative="two-sided")
u_weaker = res.statistic                  
u_stronger = len(t) * len(b) - u_weaker   
ps = u_stronger / (len(t) * len(b)) 
print(f"\nMann-Whitney U = {u_weaker:,.0f}  (complement {u_stronger:,.0f}; "
      f"n1*n2 = {len(t) * len(b):,})")
print(f"p = {res.pvalue:.2e}")
print(f"Probability that a T-cell+ peptide is the stronger predicted binder = {ps:.3f}")

OUT_DIR = Path("results")
OUT_DIR.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({
    "figure.facecolor": "white",
    "axes.facecolor": "white",
    "savefig.facecolor": "white",
    "text.color": "#1a1a1a",
    "axes.labelcolor": "#1a1a1a",
    "xtick.color": "#1a1a1a",
    "ytick.color": "#1a1a1a",
    "axes.edgecolor": "#888888",
    "font.size": 12,
})
grid = "#e3e3e3"

# Figure 1: Predicted vs experimental binding
fig, a = plt.subplots(figsize=(7.5, 6))
hb = a.hexbin(
    x, y, gridsize=45, bins="log", cmap="viridis",
    mincnt=1, extent=(-2, 4.2, 0, 4.8)
)
a.plot([-2, 4.8], [-2, 4.8], "--", color="#d62828",
       linewidth=1.8, label="y = x")
xs = np.array([-1.8, 4.0])
a.plot(xs, slope * xs + intercept, color="#e07a00", linewidth=2,
       label=f"regression fit")
a.text(
    0.03, 0.97,
    f"n = {len(x):,}\nSpearman ρ = {rho:.2f} (95% CI {lo:.2f}–{hi:.2f})\n"
    f"MAE = {mae:.2f} log10 units",
    transform=a.transAxes, va="top", fontsize=11,
    bbox=dict(facecolor="white", edgecolor="#cccccc",
              boxstyle="round,pad=0.4")
)
a.set_xlim(-2, 4.2)
a.set_ylim(0, 4.8)
a.set_xlabel("Experimental IC50, log10 (nM)")
a.set_ylabel("NetMHCpan 4.1 BA predicted IC50, log10 (nM)")
a.set_title("Predicted vs experimental binding (HLA-A*02:01, 9-mers)",
            fontsize=12.5)
a.legend(loc="lower right", frameon=True, facecolor="white",
         edgecolor="#cccccc", framealpha=0.95)
a.grid(color=grid, linewidth=0.6)
fig.colorbar(hb, ax=a, pad=0.02).set_label("peptides per bin (log)")
fig.tight_layout()
fig.savefig(OUT_DIR / "plot1_scatter.png", dpi=200)
plt.close(fig)

# Figure 2: ROC curves for 500 nM and 50 nM thresholds
fig, b = plt.subplots(figsize=(6.2, 5.6))
for threshold, color in ((500, "#1f77b4"), (50, "#d6336c")):
    is_binder = (x < np.log10(threshold)).astype(int)
    fpr, tpr, _ = roc_curve(is_binder, -y)
    auc = roc_auc_score(is_binder, -y)
    b.plot(fpr, tpr, color=color, linewidth=2,
           label=f"binder < {threshold} nM (AUC {auc:.2f})")
b.plot([0, 1], [0, 1], "--", color="#777777")
b.set_xlabel("False positive rate")
b.set_ylabel("True positive rate")
b.set_title("Binder classification from predicted IC50", fontsize=12.5)
b.legend(loc="lower right", frameon=False)
b.grid(color=grid, linewidth=0.6)
fig.tight_layout()
fig.savefig(OUT_DIR / "plot2_roc.png", dpi=200)
plt.close(fig)

# Figure 4: Cumulative distribution by evidence group (no bar chart)
colors = {
    "T-cell+": "#d6336c",
    "binding+": "#1f77b4",
    "both": "#e09f00",
    "Both": "#e09f00",
}
display_names = {
    "T-cell+": "T-cell+ only",
    "binding+": "Binding+ only",
    "both": "Both",
    "Both": "Both",
}
fig, a = plt.subplots(figsize=(7.2, 5.6))
for name, group in df.groupby(COL_EVIDENCE):
    values = group["pred_log10_ic50"].dropna().to_numpy()
    if len(values) == 0:
        continue
    values = np.sort(values)
    cumulative = np.arange(1, len(values) + 1) / len(values)
    median_nm = np.median(group[COL_PRED_NM].dropna())
    color = colors.get(name, "#777777")
    label = display_names.get(name, str(name))
    a.plot(
        values, cumulative, color=color, linewidth=2.2,
        label=f"{label} (n={len(values):,}, median {median_nm:.0f} nM)"
    )
a.axvline(np.log10(500), color="#777777", linestyle="--", linewidth=1)
a.text(np.log10(500) + 0.05, 0.3, "500 nM",
       color="#555555", fontsize=10)
a.set_xlabel("Predicted IC50, log10 (nM)")
a.set_ylabel("Cumulative fraction of peptides")
a.set_title("Predicted binding by evidence group")
a.legend(frameon=False, loc="lower right", fontsize=10.5)
a.grid(color=grid, linewidth=0.6)
fig.tight_layout()
fig.savefig(OUT_DIR / "plot3_evidence_curves.png", dpi=200)
plt.close(fig)

print("\nSaved selected figures:")
































# fig, ax = plt.subplots(1, 3, figsize=(16, 4.8))
# ax[0].hexbin(x, y, gridsize=45, bins="log", mincnt=1)
# lim = [0, 4.8]
# ax[0].plot(lim, lim, "r--", label="y = x")
# xs = np.array([-1.8, 4.0]); ax[0].plot(xs, slope * xs + intercept, color="orange", label="regression fit")
# ax[0].set(xlabel="Experimental log10 IC50 (nM)", ylabel="Predicted log10 IC50 (nM)", title="Predicted vs experimental")
# ax[0].legend()

# is_binder = (x < np.log10(500)).astype(int)
# fpr, tpr, _ = roc_curve(is_binder, -y)
# ax[1].plot(fpr, tpr, label=f"AUC = {roc_auc_score(is_binder, -y):.2f}")
# ax[1].plot([0, 1], [0, 1], "k--")
# ax[1].set(xlabel="False positive rate", ylabel="True positive rate", title="ROC (binder < 500 nM)")
# ax[1].legend()

# for name, grp in df.groupby(COL_EVIDENCE):
#     s = np.sort(grp["pred_log10_ic50"])
#     ax[2].plot(s, np.arange(1, len(s) + 1) / len(s), label=name)
# ax[2].set(xlabel="Predicted log10 IC50 (nM)", ylabel="Cumulative fraction", title="By evidence group")
# ax[2].legend()

# plt.tight_layout()
# plt.savefig("../results/figures.png", dpi=200)
# print("\nSaved figures.png")