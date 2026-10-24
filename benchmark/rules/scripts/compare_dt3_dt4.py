import re

import numpy as np
import pandas as pd

mode = snakemake.params.mode
tool = snakemake.params.tool
sample = snakemake.params.sample

row = {"tool": tool, "sample": sample}

CHUNK_SIZE = 50_000_000
CANONICAL = re.compile(r"^(chr)?(\d{1,2}|X|Y|M|MT)$|^[1-7][ABD]$")


def is_canonical(name):
    return CANONICAL.match(str(name)) is not None


row["canonical_only"] = True
HIST_EDGES = np.logspace(-9, 4, 13 * 4 + 1)
REL_EDGES = np.logspace(-12, 0, 12 * 4 + 1)  # |dt3 - dt4| / max(|dt3|, |dt4|) is at most 1
THRESHOLDS = (1e-6, 1e-3, 1e-1)


class DiffStats:

    def __init__(self):
        self.n = 0
        self.n_missing_mismatch = 0
        self.n_nonzero = 0
        self.sx = self.sy = self.sxx = self.syy = self.sxy = 0.0
        self.max_abs = 0.0
        self.max_abs_loc = None
        self.max_rel = 0.0
        self.counts = np.zeros(len(HIST_EDGES) + 1, dtype=np.int64)  # under, bins..., over
        self.rel_counts = np.zeros(len(REL_EDGES) + 1, dtype=np.int64)
        self.n_zero = 0
        self.n_gt = {t: 0 for t in THRESHOLDS}

    def update(self, x, y, loc=lambda flat_idx: flat_idx):
        x = np.asarray(x, dtype=np.float64).ravel()
        y = np.asarray(y, dtype=np.float64).ravel()
        nx, ny = np.isnan(x), np.isnan(y)
        self.n_missing_mismatch += int((nx != ny).sum())
        keep = ~(nx | ny)
        if not keep.all():
            idx_map = np.flatnonzero(keep)
            x, y = x[keep], y[keep]
        else:
            idx_map = None
        if x.size == 0:
            return

        self.n += x.size
        self.sx += x.sum()
        self.sy += y.sum()
        self.sxx += (x * x).sum()
        self.syy += (y * y).sum()
        self.sxy += (x * y).sum()

        d = np.abs(x - y)
        k = int(d.argmax())
        if d[k] > self.max_abs:
            self.max_abs = float(d[k])
            flat = int(idx_map[k]) if idx_map is not None else k
            self.max_abs_loc = loc(flat)

        nz = d > 0
        n_nz = int(nz.sum())
        self.n_nonzero += n_nz
        self.n_zero += d.size - n_nz
        if n_nz:
            dn = d[nz]
            denom = np.maximum(np.abs(x[nz]), np.abs(y[nz]))
            rel = dn / denom
            self.max_rel = max(self.max_rel, float(rel.max()))
            self.counts += np.bincount(np.searchsorted(HIST_EDGES, dn), minlength=len(self.counts))
            self.rel_counts += np.bincount(np.searchsorted(REL_EDGES, rel), minlength=len(self.rel_counts))
            for t in THRESHOLDS:
                self.n_gt[t] += int((dn > t).sum())

    def finish(self):
        n = self.n
        if n == 0:
            pearson_r = float("nan")
        else:
            cov = self.sxy / n - (self.sx / n) * (self.sy / n)
            var_x = self.sxx / n - (self.sx / n) ** 2
            var_y = self.syy / n - (self.sy / n) ** 2
            denom = np.sqrt(var_x * var_y)
            pearson_r = float(cov / denom) if denom > 0 else float("nan")
        out = {
            "n_values": n,
            "pearson_r": pearson_r,
            "n_nan_mismatch": self.n_missing_mismatch,
            "n_nonidentical": self.n_nonzero,
            "frac_identical": self.n_zero / n if n else float("nan"),
            "max_abs_diff": self.max_abs,
            "max_abs_diff_loc": self.max_abs_loc,
            "max_rel_diff": self.max_rel,
        }
        for t in THRESHOLDS:
            out[f"n_diff_gt_{t:g}"] = self.n_gt[t]
        return out

    def histogram(self, relative=False):
        edges, counts = (REL_EDGES, self.rel_counts) if relative else (HIST_EDGES, self.counts)
        lo = np.concatenate([[0.0], [0.0], edges[:-1], [edges[-1]]])
        hi = np.concatenate([[0.0], [edges[0]], edges[1:], [np.inf]])
        counts = np.concatenate([[self.n_zero], counts])
        return pd.DataFrame({"lower": lo, "upper": hi, "count": counts})


def write_hist(stats):
    if "hist" in snakemake.output.keys():
        stats.histogram().to_csv(snakemake.output.hist, sep="\t", index=False)
    if "relhist" in snakemake.output.keys():
        stats.histogram(relative=True).to_csv(snakemake.output.relhist, sep="\t", index=False)


if mode == "bigwig":
    import pyBigWig

    def paired_chunks(bw3, bw4):
        for chrom, length in bw3.chroms().items():
            if chrom not in bw4.chroms() or not is_canonical(chrom):
                continue
            for start in range(0, length, CHUNK_SIZE):
                end = min(start + CHUNK_SIZE, length)
                v3 = np.array(bw3.values(chrom, start, end, numpy=True), dtype=np.float32)
                v4 = np.array(bw4.values(chrom, start, end, numpy=True), dtype=np.float32)
                yield chrom, start, v3, v4

    bw3 = pyBigWig.open(snakemake.input.dt3)
    bw4 = pyBigWig.open(snakemake.input.dt4)
    row["chroms_only_dt3"] = ",".join(c for c in bw3.chroms() if c not in bw4.chroms())
    row["chroms_only_dt4"] = ",".join(c for c in bw4.chroms() if c not in bw3.chroms())

    stats = DiffStats()
    for chrom, start, v3, v4 in paired_chunks(bw3, bw4):
        stats.update(v3, v4, loc=lambda i, c=chrom, s=start: f"{c}:{s + i}")
    bw3.close()
    bw4.close()

    row.update(stats.finish())
    write_hist(stats)


elif mode == "npz":
    d3 = pd.read_csv(snakemake.input.dt3, sep="\t", header=0, low_memory=False)
    d4 = pd.read_csv(snakemake.input.dt4, sep="\t", header=0, low_memory=False)

    row["shape_dt3"] = str(d3.shape)
    row["shape_dt4"] = str(d4.shape)
    row["shape_mismatch"] = d3.shape != d4.shape

    for d in (d3, d4):
        d.columns = ["chrom", "start", "end"] + [f"s{i}" for i in range(d.shape[1] - 3)]
        d["chrom"] = d["chrom"].astype(str)
    merged = d3.merge(d4, on=["chrom", "start", "end"], suffixes=("_3", "_4"))
    canon = merged["chrom"].map(is_canonical)
    row["n_noncanonical_bins_excluded"] = int((~canon).sum())
    merged = merged[canon]
    ncol = d3.shape[1] - 3
    v3 = merged[[f"s{i}_3" for i in range(ncol)]].to_numpy(dtype=np.float64)
    v4 = merged[[f"s{i}_4" for i in range(ncol)]].to_numpy(dtype=np.float64)

    row["n_bins_dt3"] = len(d3)
    row["n_bins_dt4"] = len(d4)
    row["n_common_bins"] = len(merged)

    coords = merged[["chrom", "start", "end"]].to_numpy()

    def loc(i):
        r, c = divmod(i, ncol)
        return f"{coords[r][0]}:{coords[r][1]}-{coords[r][2]} sample{c}"

    stats = DiffStats()
    stats.update(v3, v4, loc=loc)
    row.update(stats.finish())
    write_hist(stats)


elif mode == "matrix":
    from deeptools.heatmapper import heatmapper

    ROW_CHUNK = 20_000

    hm3 = heatmapper()
    hm3.read_matrix_file(snakemake.input.dt3)
    hm4 = heatmapper()
    hm4.read_matrix_file(snakemake.input.dt4)
    m3 = hm3.matrix.matrix
    m4 = hm4.matrix.matrix
    row["shape_dt3"] = str(m3.shape)
    row["shape_dt4"] = str(m4.shape)
    row["shape_mismatch"] = m3.shape != m4.shape

    n_rows = min(m3.shape[0], m4.shape[0])
    ncol = min(m3.shape[1], m4.shape[1])
    # rows are regions; keep canonical chromosomes only, and check both versions list the same regions
    r3, r4 = hm3.matrix.regions[:n_rows], hm4.matrix.regions[:n_rows]
    keep = np.array([is_canonical(r[0]) for r in r3], dtype=bool)
    row["n_noncanonical_rows_excluded"] = int((~keep).sum())
    row["n_region_mismatch"] = sum(
        (a[0], a[1], a[2]) != (b[0], b[1], b[2]) for a, b, k in zip(r3, r4, keep) if k
    )
    rows = np.flatnonzero(keep)
    stats = DiffStats()
    for a in range(0, len(rows), ROW_CHUNK):
        idx = rows[a:a + ROW_CHUNK]
        c3 = np.ma.filled(m3[idx, :ncol].astype(np.float64), np.nan)
        c4 = np.ma.filled(m4[idx, :ncol].astype(np.float64), np.nan)
        stats.update(c3, c4, loc=lambda i, idx=idx: f"row{idx[i // ncol]} col{i % ncol}")

    row.update(stats.finish())
    write_hist(stats)


elif mode == "bam":
    import os
    import pysam
    from collections import Counter
    from multiprocessing import Pool

    COUNTERS = ["n3", "n4", "n_identical", "n_only_dt3", "n_only_dt4",
                "n_flag_mismatch", "n_qlen_mismatch", "n_cigar_mismatch", "n_mapq_mismatch"]

    def groups(f, contig):
        """Yield [(key, info), ...] for each distinct start position on contig (coordinate order)."""
        cur, buf = None, []
        for read in f.fetch(contig):
            if read.is_unmapped:
                continue
            if read.reference_start != cur:
                if buf:
                    yield cur, buf
                cur, buf = read.reference_start, []
            buf.append((
                (read.query_name, bool(read.flag & 64)),
                (read.flag, read.query_length, read.cigarstring, read.mapping_quality),
            ))
        if buf:
            yield cur, buf

    def compare_group(b3, b4, c):
        """Match the reads at one start position, updating the counter dict c."""
        c3 = Counter(b3)
        c4 = Counter(b4)
        same = c3 & c4
        c["n_identical"] += sum(same.values())
        left3 = c3 - same
        left4 = c4 - same
        by_key4 = {}
        for (k, i), n in left4.items():
            by_key4.setdefault(k, []).extend([i] * n)
        for (k, i3), n in left3.items():
            for _ in range(n):
                cand = by_key4.get(k)
                if cand:
                    i4 = cand.pop()
                    c["n_flag_mismatch"] += i3[0] != i4[0]
                    c["n_qlen_mismatch"] += i3[1] != i4[1]
                    c["n_cigar_mismatch"] += i3[2] != i4[2]
                    c["n_mapq_mismatch"] += i3[3] != i4[3]
                else:
                    c["n_only_dt3"] += 1
        c["n_only_dt4"] += sum(len(v) for v in by_key4.values())

    def compare_contig(args):
        path3, path4, contig = args
        c = dict.fromkeys(COUNTERS, 0)
        with pysam.AlignmentFile(path3, "rb") as f3, pysam.AlignmentFile(path4, "rb") as f4:
            g3, g4 = groups(f3, contig), groups(f4, contig)
            a, b = next(g3, None), next(g4, None)
            while a is not None or b is not None:
                if b is None or (a is not None and a[0] < b[0]):
                    c["n3"] += len(a[1]); c["n_only_dt3"] += len(a[1]); a = next(g3, None)
                elif a is None or b[0] < a[0]:
                    c["n4"] += len(b[1]); c["n_only_dt4"] += len(b[1]); b = next(g4, None)
                else:
                    c["n3"] += len(a[1]); c["n4"] += len(b[1])
                    compare_group(a[1], b[1], c)
                    a, b = next(g3, None), next(g4, None)
        return contig, c
    sorted_paths = {}
    for key in ("dt3", "dt4"):
        src = snakemake.input[key]
        dst = f"{snakemake.output.tsv}.{key}.sorted.bam"
        pysam.sort("-@", str(snakemake.threads), "-m", "500M", "-o", dst, src)
        pysam.index("-c", dst)
        sorted_paths[key] = dst
    sorted3, sorted4 = sorted_paths["dt3"], sorted_paths["dt4"]

    with pysam.AlignmentFile(sorted3, "rb") as f3, pysam.AlignmentFile(sorted4, "rb") as f4:
        if f3.references != f4.references:
            raise RuntimeError("dt3 and dt4 BAMs have different reference sequences")
        order = sorted(((n, l) for n, l in zip(f3.references, f3.lengths) if is_canonical(n)),
                       key=lambda t: -t[1])
        jobs = [(sorted3, sorted4, name) for name, _ in order]

    total = dict.fromkeys(COUNTERS, 0)
    with Pool(snakemake.threads) as pool:
        for _, c in pool.imap_unordered(compare_contig, jobs):
            for k in COUNTERS:
                total[k] += c[k]

    for path in sorted_paths.values():
        os.remove(path)
        os.remove(path + ".csi")

    row["n_reads_dt3"] = total["n3"]
    row["n_reads_dt4"] = total["n4"]
    for k in COUNTERS[2:]:
        row[k] = total[k]


else:
    raise ValueError(f"unknown comparison mode: {mode}")

pd.DataFrame([row]).to_csv(snakemake.output.tsv, sep="\t", index=False)
