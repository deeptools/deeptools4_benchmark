rule plot_benchmarks:
  localrule: True
  input:
    alsieve = bench_inputs("alignmentsieve", alignmentSieve_samples, ALSIEVE_THREADS),
    bcov = bench_inputs("bamcoverage", bamCoverage_samples, BAMCOVERAGE_THREADS),
    bcom = bench_inputs("bamcompare", bamCompare_samples, BAMCOMPARE_THREADS),
    mbs = bench_inputs("multibamsummary", multibamSummary_samples, MULTIBAMSUMMARY_THREADS),
    cm = bench_inputs("computeMatrix", computeMatrix_samples, COMPUTEMATRIX_THREADS)
  output:
    csv = 'results/performance.csv',
    png = 'results/performance.png',
    tiff = 'results/performance.tiff',
    pdf = 'results/performance.pdf'
  params:
    os = config['os']
  conda: 'env/plotter.yaml'
  threads: 2
  script:
    'scripts/plotter.py'

rule plot_benchmarks_suppfigs:
  localrule: True
  input:
    csv = 'results/performance.csv',
  output:
    runpng = 'results/performance_runtime.png',
    runpdf = 'results/performance_runtime.pdf',
    runtiff = 'results/performance_runtime.tiff',
    mempng = 'results/performance_memory.png',
    mempdf = 'results/performance_memory.pdf',
    memtiff = 'results/performance_memory.tiff'
  conda: 'env/plotter.yaml'
  threads: 2
  script:
    'scripts/benchmark_supfigs.py'
