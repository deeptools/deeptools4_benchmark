import re
import numpy as np
import pandas as pd

_stem_re = re.compile(r'_(dt[34])_t(\d+)_rep(\d+)$')

def read_benchmark(file):
    stem = file.split('/')[-1].split('.')[0]
    m = _stem_re.search(stem)
    a = pd.read_table(file, sep='\t')
    a = a[['s', 'max_rss']]
    a['rep'] = int(m.group(3))
    a['version'] = m.group(1)
    a['threads'] = int(m.group(2))
    a['sample'] = stem[:m.start()]
    a['tool'] = file.split('/')[-2]
    return a

df = pd.concat([read_benchmark(f) for f in snakemake.input], ignore_index=True)

# memory usage change (dt4 vs dt3), averaged over all threads/reps per sample
mem = (
    df.groupby(['tool', 'sample', 'version'])['max_rss']
      .mean()
      .unstack('version')
      .dropna(subset=['dt3', 'dt4'])
)
mem['dt3_gb'] = mem['dt3'] / 1024
mem['dt4_gb'] = mem['dt4'] / 1024
mem['mem_diff_gb'] = (mem['dt4'] - mem['dt3']) / 1024
mem['mem_diff_pct'] = 100 * (mem['dt4'] - mem['dt3']) / mem['dt3']
mem_summary = mem.groupby('tool')[['dt3_gb', 'dt4_gb', 'mem_diff_pct', 'mem_diff_gb']].mean()
mem_summary.columns = ['mean_mem_dt3_gb', 'mean_mem_dt4_gb', 'mean_mem_diff_pct', 'mean_mem_diff_gb']

# runtime speedup (dt3 / dt4) over all thread counts:
time = (
    df.groupby(['tool', 'sample', 'threads', 'version'])['s']
      .mean()
      .unstack('version')
      .dropna(subset=['dt3', 'dt4'])
)
time['speedup'] = time['dt3'] / time['dt4']
time['log_speedup'] = np.log(time['speedup'])
grp = time.groupby('tool')

time = time.reset_index()
time['genome'] = time['sample'].str.extract(r'(human|triticum)', expand=False)
if time['genome'].isna().any():
    raise ValueError(f"Cannot infer genome for samples: {time.loc[time['genome'].isna(), 'sample'].unique()}")
tg = time.groupby(['tool', 'threads', 'genome'])['speedup'].agg(['min', 'max']).groupby('tool')

speed_summary = pd.DataFrame({
    'geomean_speedup': np.exp(grp['log_speedup'].mean()),
    'median_speedup': grp['speedup'].median(),
    'mean_min_speedup': tg['min'].mean(),
    'mean_max_speedup': tg['max'].mean(),
    'n_samples': time.groupby('tool')['sample'].nunique(),
    'n_threads': time.groupby('tool')['threads'].nunique(),
})

summary = mem_summary.join(speed_summary, how='outer').reset_index()
summary = summary[['tool', 'geomean_speedup', 'median_speedup', 'mean_min_speedup', 'mean_max_speedup',
                   'n_samples', 'n_threads',
                   'mean_mem_dt3_gb', 'mean_mem_dt4_gb', 'mean_mem_diff_pct', 'mean_mem_diff_gb']]
summary = summary.round(2)
summary.to_csv(snakemake.output.tsv, sep='\t', index=False)
