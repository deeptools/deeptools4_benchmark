VERIFY_THREADS = 20

def verify_pair(tool, sample, ext):
    return {
        "dt3": f"output/verify/{tool}_{sample}.dt3.{ext}",
        "dt4": f"output/verify/{tool}_{sample}.dt4.{ext}",
    }

rule verify_alsieve_dt3:
  input: bam = lambda wc: f"bamfiles/{alignmentSieve_samples[wc.s]}.bam"
  output: bam = temp("output/verify/alignmentsieve_{s}.dt3.bam")
  threads: VERIFY_THREADS
  resources: mem_mb = 10000, runtime = 1440
  shell: "ulimit -n \"$(ulimit -Hn)\"; alignmentSieve_old -b {input.bam} -o {output.bam} -p {threads} --ATACshift"

rule verify_alsieve_dt4:
  input: bam = lambda wc: f"bamfiles/{alignmentSieve_samples[wc.s]}.bam"
  output: bam = temp("output/verify/alignmentsieve_{s}.dt4.bam")
  threads: VERIFY_THREADS
  resources: mem_mb = 10000, runtime = 1440
  shell: "alignmentSieve -b {input.bam} -o {output.bam} -p {threads} --ATACshift"

rule compare_alignmentsieve:
  input: unpack(lambda wc: verify_pair("alignmentsieve", wc.s, "bam"))
  output: tsv = "results/numerical_diff/alignmentsieve_{s}.tsv"
  params: mode = "bam", tool = "alignmentSieve", sample = lambda wc: wc.s
  threads: VERIFY_THREADS
  resources: mem_mb = 20000, runtime = 1440
  script: "scripts/compare_dt3_dt4.py"

rule verify_bamcoverage_dt3:
  input: bam = lambda wc: f"bamfiles/{bamCoverage_samples[wc.s]}.bam"
  output: bw = temp("output/verify/bamcoverage_{s}.dt3.bw")
  threads: VERIFY_THREADS
  resources: mem_mb = 10000, runtime = 1440
  shell: "bamCoverage_old -b {input.bam} -o {output.bw} --binSize 10 --normalizeUsing RPKM -p {threads}"

rule verify_bamcoverage_dt4:
  input: bam = lambda wc: f"bamfiles/{bamCoverage_samples[wc.s]}.bam"
  output: bw = temp("output/verify/bamcoverage_{s}.dt4.bw")
  threads: VERIFY_THREADS
  resources: mem_mb = 16000, runtime = 1440
  shell: "bamCoverage -b {input.bam} -o {output.bw} --binSize 10 --normalizeUsing RPKM -p {threads}"

rule compare_bamcoverage:
  input: unpack(lambda wc: verify_pair("bamcoverage", wc.s, "bw"))
  output:
    tsv = "results/numerical_diff/bamcoverage_{s}.tsv",
    hist = "results/numerical_diff/bamcoverage_{s}.hist.tsv"
  params: mode = "bigwig", tool = "bamCoverage", sample = lambda wc: wc.s
  resources: mem_mb = 20000, runtime = 1440
  script: "scripts/compare_dt3_dt4.py"

rule verify_bamcompare_dt3:
  input:
    bam = lambda wc: f"bamfiles/{bamCompare_samples[wc.s][0]}.bam",
    ctrlbam = lambda wc: f"bamfiles/{bamCompare_samples[wc.s][1]}.bam"
  output: bw = temp("output/verify/bamcompare_{s}.dt3.bw")
  threads: VERIFY_THREADS
  resources: mem_mb = 20000, runtime = 1440
  shell: "bamCompare_old -b1 {input.bam} -b2 {input.ctrlbam} -o {output.bw} --binSize 10 -p {threads}"

rule verify_bamcompare_dt4:
  input:
    bam = lambda wc: f"bamfiles/{bamCompare_samples[wc.s][0]}.bam",
    ctrlbam = lambda wc: f"bamfiles/{bamCompare_samples[wc.s][1]}.bam"
  output: bw = temp("output/verify/bamcompare_{s}.dt4.bw")
  threads: VERIFY_THREADS
  resources: mem_mb = 20000, runtime = 1440
  shell: "bamCompare -b1 {input.bam} -b2 {input.ctrlbam} -o {output.bw} --binSize 10 -p {threads}"

rule compare_bamcompare:
  input: unpack(lambda wc: verify_pair("bamcompare", wc.s, "bw"))
  output:
    tsv = "results/numerical_diff/bamcompare_{s}.tsv",
    hist = "results/numerical_diff/bamcompare_{s}.hist.tsv"
  params: mode = "bigwig", tool = "bamCompare", sample = lambda wc: wc.s
  resources: mem_mb = 20000, runtime = 1440
  script: "scripts/compare_dt3_dt4.py"

rule verify_multibamsummary_dt3:
  input: bam_files = lambda wc: [f"bamfiles/{i}.bam" for i in multibamSummary_samples[wc.s]]
  output:
    npz = temp("output/verify/multibamsummary_{s}.dt3.npz"),
    tsv = temp("output/verify/multibamsummary_{s}.dt3.tsv")
  threads: VERIFY_THREADS
  resources: mem_mb = 20000, runtime = 1440
  shell: "multiBamSummary_old bins -p {threads} -o {output.npz} --outRawCounts {output.tsv} -b {input.bam_files} --genomeChunkSize 10000000"

rule verify_multibamsummary_dt4:
  input: bam_files = lambda wc: [f"bamfiles/{i}.bam" for i in multibamSummary_samples[wc.s]]
  output:
    npz = temp("output/verify/multibamsummary_{s}.dt4.npz"),
    tsv = temp("output/verify/multibamsummary_{s}.dt4.tsv")
  threads: VERIFY_THREADS
  resources: mem_mb = 20000, runtime = 1440
  shell: "multiBamSummary bins -p {threads} -o {output.npz} --outRawCounts {output.tsv} -b {input.bam_files}"

rule compare_multibamsummary:
  input: unpack(lambda wc: verify_pair("multibamsummary", wc.s, "tsv"))
  output:
    tsv = "results/numerical_diff/multibamsummary_{s}.tsv",
    hist = "results/numerical_diff/multibamsummary_{s}.hist.tsv"
  params: mode = "npz", tool = "multiBamSummary", sample = lambda wc: wc.s
  resources: mem_mb = 40000, runtime = 1440
  script: "scripts/compare_dt3_dt4.py"

rule verify_computematrix_dt3:
  input:
    bw_files = lambda wc: [f"zenodo/{i}.bw" for i in computeMatrix_samples[wc.s]],
    gtf = lambda wc: f"zenodo/{wc.s.split('_')[1]}.gtf"
  output: npz = temp("output/verify/computematrix_{s}.dt3.npz")
  threads: VERIFY_THREADS
  resources: mem_mb = 20000, runtime = 1440
  shell: "computeMatrix_old reference-point -p {threads} -o {output.npz} -R {input.gtf} -S {input.bw_files}"

rule verify_computematrix_dt4:
  input:
    bw_files = lambda wc: [f"zenodo/{i}.bw" for i in computeMatrix_samples[wc.s]],
    gtf = lambda wc: f"zenodo/{wc.s.split('_')[1]}.gtf"
  output: npz = temp("output/verify/computematrix_{s}.dt4.npz")
  threads: VERIFY_THREADS
  resources: mem_mb = 40000, runtime = 1440
  shell: "computeMatrix reference-point -p {threads} -o {output.npz} -R {input.gtf} -S {input.bw_files}"

rule compare_computematrix:
  input: unpack(lambda wc: verify_pair("computematrix", wc.s, "npz"))
  output:
    tsv = "results/numerical_diff/computematrix_{s}.tsv",
    hist = "results/numerical_diff/computematrix_{s}.hist.tsv"
  params: mode = "matrix", tool = "computeMatrix", sample = lambda wc: wc.s
  resources: mem_mb = 60000, runtime = 1440
  script: "scripts/compare_dt3_dt4.py"
