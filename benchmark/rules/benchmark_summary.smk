rule benchmark_summary:
  localrule: True
  input:
    alsieve = bench_inputs("alignmentsieve", alignmentSieve_samples, ALSIEVE_THREADS),
    bcov = bench_inputs("bamcoverage", bamCoverage_samples, BAMCOVERAGE_THREADS),
    bcom = bench_inputs("bamcompare", bamCompare_samples, BAMCOMPARE_THREADS),
    mbs = bench_inputs("multibamsummary", multibamSummary_samples, MULTIBAMSUMMARY_THREADS),
    cm = bench_inputs("computeMatrix", computeMatrix_samples, COMPUTEMATRIX_THREADS)
  output:
    tsv = 'results/benchmark_summary.tsv'
  conda: 'env/plotter.yaml'
  script:
    'scripts/benchmark_summary.py'

rule combine_numerical_diff:
  localrule: True
  input:
    expand("results/numerical_diff/alignmentsieve_{s}.tsv", s=alignmentSieve_samples.keys()),
    expand("results/numerical_diff/bamcoverage_{s}.tsv", s=bamCoverage_samples.keys()),
    expand("results/numerical_diff/bamcompare_{s}.tsv", s=VERIFY_bamCompare_samples.keys()),
    expand("results/numerical_diff/multibamsummary_{s}.tsv", s=VERIFY_multibamSummary_samples.keys()),
    expand("results/numerical_diff/computematrix_{s}.tsv", s=VERIFY_computeMatrix_samples.keys())
  output:
    tsv = 'results/numerical_diff.tsv'
  run:
    import pandas as pd
    df = pd.concat([pd.read_csv(f, sep="\t") for f in input], ignore_index=True)
    df.to_csv(output.tsv, sep="\t", index=False)

rule plot_numerical_diff:
  localrule: True
  input:
    tsv = 'results/numerical_diff.tsv',
    hists = lambda wc: (
      expand("results/numerical_diff/bamcoverage_{s}.hist.tsv", s=bamCoverage_samples.keys())
      + expand("results/numerical_diff/bamcompare_{s}.hist.tsv", s=VERIFY_bamCompare_samples.keys())
      + expand("results/numerical_diff/multibamsummary_{s}.hist.tsv", s=VERIFY_multibamSummary_samples.keys())
      + expand("results/numerical_diff/computematrix_{s}.hist.tsv", s=VERIFY_computeMatrix_samples.keys())
    ),
    relhists = lambda wc: (
      expand("results/numerical_diff/bamcoverage_{s}.relhist.tsv", s=bamCoverage_samples.keys())
      + expand("results/numerical_diff/bamcompare_{s}.relhist.tsv", s=VERIFY_bamCompare_samples.keys())
      + expand("results/numerical_diff/multibamsummary_{s}.relhist.tsv", s=VERIFY_multibamSummary_samples.keys())
      + expand("results/numerical_diff/computematrix_{s}.relhist.tsv", s=VERIFY_computeMatrix_samples.keys())
    ),
  output:
    png = 'results/numerical_diff_hist.png',
    pdf = 'results/numerical_diff_hist.pdf',
    tiff = 'results/numerical_diff_hist.tiff'
  conda: 'env/plotter.yaml'
  script:
    'scripts/plot_numerical_diff.py'
