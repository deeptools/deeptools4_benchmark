import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

diff = pd.read_table(snakemake.input.tsv)
tools = ["bamCoverage", "bamCompare", "multiBamSummary", "computeMatrix"]
prefix = {"bamCoverage": "bamcoverage", "bamCompare": "bamcompare",
          "multiBamSummary": "multibamsummary", "computeMatrix": "computematrix"}


def read_hists(files, suffix):
    return {f.rsplit("/", 1)[1].removesuffix(suffix): pd.read_table(f) for f in files}


hists = read_hists(snakemake.input.hists, ".hist.tsv")
relhists = read_hists(snakemake.input.relhists, ".relhist.tsv")

colors = {"bamCoverage": "#C76B95", "bamCompare": "#5B8FC7", "multiBamSummary": "#E0A030", "computeMatrix": "#4DA37E"}
organisms = {"human": "Human", "triticum": "Triticum"}

STEP = 10 ** 0.25
ymin = 1e-11
rel_override = {"bamCompare": (hists, lambda x: np.expm1(np.minimum(x, 64) * np.log(2)))}
panels = [
    (hists, -8, 1, "abs. diff", {}),
    (relhists, -12, 0, "rel. diff", rel_override),
]


def run_fracs(tool, org, hs, override):
    sub = diff[(diff["tool"] == tool) & diff["sample"].str.contains(org)]
    fracs = []
    for _, r in sub.iterrows():
        key = f"{prefix[tool]}_{r['sample']}"
        grid = hs[key]
        if tool in override:
            src, fn = override[tool]
            h = src[key]
            hi = h["upper"].to_numpy()
            lo = np.concatenate([[0, hi[1] / STEP], hi[1:-1]])
            mid = fn(np.sqrt(lo * np.where(np.isinf(hi), lo * STEP, hi)))
            rows = np.concatenate([[0], np.searchsorted(grid["upper"].to_numpy()[1:], mid[1:]) + 1])
            counts = np.bincount(rows, weights=h["count"].to_numpy(), minlength=len(grid))
        else:
            counts = grid["count"].to_numpy()
        fracs.append(counts / counts.sum())
    return np.mean(fracs, axis=0), grid["upper"].to_numpy(), len(fracs)


def draw(ax, hs, org, first, last, override):
    zero_lo, zero_hi = 10.0 ** (first - 2), 10.0 ** (first - 1.5)
    for k, tool in enumerate(tools):
        frac, hi, n = run_fracs(tool, org, hs, override)
        lo = np.concatenate([[0, hi[1] / STEP], hi[1:-1]])
        mid = np.sqrt(lo[1:-1] * hi[1:-1])
        f = frac[1:-1]
        keep = f > 0
        x = np.concatenate([[mid[keep][0] / STEP], mid[keep], [mid[keep][-1] * STEP]])
        y = np.concatenate([[ymin], f[keep], [ymin]])
        ax.fill_between(x, ymin, y, color=colors[tool], alpha=0.2, lw=0)
        ax.plot(x, y, color=colors[tool], lw=1, label=f"{tool} (n={n})")
        slot = (zero_hi / zero_lo) ** (1 / len(tools))
        ax.bar(zero_lo * slot ** k, frac[0] - ymin, width=zero_lo * slot ** (k + 1) - zero_lo * slot ** k,
               bottom=ymin, align="edge", color=colors[tool], lw=0)
    ax.axvline(10.0 ** (first - 1.4), color="0.6", lw=0.8, ls=":")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylim(bottom=ymin, top=2)
    ax.set_xlim(zero_lo * 0.8, 10.0 ** (last + 0.1) if last >= 1 else 2)
    exps = list(range(first, last + 1, 2))
    ax.set_xticks([10 ** (first - 1.75)] + [10.0 ** e for e in exps])
    ax.set_xticklabels(["0"] + [f"$10^{{{e}}}$" for e in exps])
    ax.minorticks_off()


fig, axes = plt.subplots(2, 2, figsize=(14, 7), sharey=True)
for i, (org, title) in enumerate(organisms.items()):
    for j, (hs, first, last, xlabel, override) in enumerate(panels):
        ax = axes[i, j]
        draw(ax, hs, org, first, last, override)
        ax.set_title(title)
        if j == 0:
            ax.set_ylabel("fraction of values")
        if i == 1:
            ax.set_xlabel(xlabel)
        else:
            ax.tick_params(labelbottom=False)
handles, labels = axes[0, 0].get_legend_handles_labels()
fig.legend(handles, labels, loc="lower center", ncol=len(tools), fontsize=8, frameon=False)
fig.tight_layout(rect=(0, 0.04, 1, 1))
fig.savefig(snakemake.output.png, dpi=300)
fig.savefig(snakemake.output.pdf)
fig.savefig(snakemake.output.tiff, dpi=300)
