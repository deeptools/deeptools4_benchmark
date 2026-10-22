import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

diff = pd.read_table(snakemake.input.tsv)
tools = ["bamCoverage", "bamCompare", "multiBamSummary", "computeMatrix"]
prefix = {"bamCoverage": "bamcoverage", "bamCompare": "bamcompare",
          "multiBamSummary": "multibamsummary", "computeMatrix": "computematrix"}
hists = {}
for f in snakemake.input.hists:
    name = f.rsplit("/", 1)[1].removesuffix(".hist.tsv")
    hists[name] = pd.read_table(f)

colors = {"human": "#C76B95", "triticum": "#5B8FC7"}

fig, axes = plt.subplots(2, 2, figsize=(11, 8), sharey=True)
for ax, tool in zip(axes.ravel(), tools):
    sub = diff[diff["tool"] == tool]
    for _, r in sub.iterrows():
        h = hists[f"{prefix[tool]}_{r['sample']}"]
        total = h["count"].sum()
        nz = h.iloc[1:]  # first row holds the exact zeros
        mid = np.sqrt(np.where(nz["lower"] > 0, nz["lower"], nz["upper"]) * np.where(np.isinf(nz["upper"]), nz["lower"] * 2, nz["upper"]))
        org = "human" if "human" in r["sample"] else "triticum"
        ax.plot(mid, nz["count"] / total, color=colors[org], alpha=0.8, marker="o", ms=2.5, lw=1,
                label=r["sample"])
        ax.axvline(r["max_abs_diff"], color=colors[org], ls=":", lw=1)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_title(tool)
    ax.set_xlabel("|dt3 - dt4|")
    ax.legend(fontsize=6, frameon=False)
for ax in axes[:, 0]:
    ax.set_ylabel("fraction of values")
fig.suptitle("Element-wise differences between deepTools 3 and 4 (non-zero differences; dotted = max)", fontsize=10)
fig.tight_layout()
fig.savefig(snakemake.output.png, dpi=300)
fig.savefig(snakemake.output.pdf)
fig.savefig(snakemake.output.tiff, dpi=300)
