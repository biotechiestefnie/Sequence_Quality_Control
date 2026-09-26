from snakemake.io import expand

samples = ["SRR25083113"]
reads = ["1", "2"]

rule all:
    input:
        # 1. fastqc assessment on raw baseline reads
        expand("results/fastqc/baseline/{sample}_{read}_baseline_fastqc.html", sample=samples, read=reads),
        # 2. fully processed and decontaminated alignment file
        expand("results/aligned/{sample}_processed.bam", sample=samples),
        # 3. index map for the processed alignment file
        expand("results/aligned/{sample}_processed.bam.bai", sample=samples),
        # 4. fastqc assessment on processed paired end reads
        expand("results/fastqc/processed/{sample}_processed_fastqc.html", sample=samples),
        # 5. parse fastqc files into single unified dashboard
        "results/multiqc_report.html"

# run quality control on raw fastq files and rename to baseline benchmarks
rule fastqc_raw:
    input:
        "data/{sample}_{read}.fastq"
    output:
        html = "results/fastqc/baseline/{sample}_{read}_baseline_fastqc.html",
        zip = "results/fastqc/baseline/{sample}_{read}_baseline_fastqc.zip"
    threads: 2
    shell:
        """
        fastqc -t {threads} {input} -o results/fastqc/baseline/ && \
        mv results/fastqc/baseline/{wildcards.sample}_{wildcards.read}_fastqc.html {output.html} && \
        mv results/fastqc/baseline/{wildcards.sample}_{wildcards.read}_fastqc.zip {output.zip}
        """

# step 1: trim adapters and low quality bases using trimmomatic
rule trim_paired:
    input:
        r1 = "data/{sample}_1.fastq",
        r2 = "data/{sample}_2.fastq"
    output:
        r1_paired = "results/trimmed/{sample}_1_paired.fastq.gz",
        r1_unpaired = "results/trimmed/{sample}_1_unpaired.fastq.gz",
        r2_paired = "results/trimmed/{sample}_2_paired.fastq.gz",
        r2_unpaired = "results/trimmed/{sample}_2_unpaired.fastq.gz"
    shell:
        """
        trimmomatic PE {input.r1} {input.r2} \
        {output.r1_paired} {output.r1_unpaired} \
        {output.r2_paired} {output.r2_unpaired} \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
        """

# step 2: remove prepackaged phix k-mer contaminants using bbduk.sh
rule bbduk_decon:
    input:
        r1 = "results/trimmed/{sample}_1_paired.fastq.gz",
        r2 = "results/trimmed/{sample}_2_paired.fastq.gz"
    output:
        r1_decon = "results/decon/{sample}_1_decon.fastq.gz",
        r2_decon = "results/decon/{sample}_2_decon.fastq.gz"
    shell:
        """
        bbduk.sh in1={input.r1} in2={input.r2} \
        out1={output.r1_decon} out2={output.r2_decon} \
        ref=phix k=31 hdist=1 rcomp=t tbo
        """

# step 3: align trimmed, decontaminated reads to indexed reference genome
rule align_paired:
    input:
        r1 = "results/trimmed/{sample}_1_decon.fastq.gz",
        r2 = "results/trimmed/{sample}_2_decon.fastq.gz",
        reference = "data/{sample}_reference.fa"
    output:
        "results/aligned/{sample}.bam"  # should this have a variable assignment like the others?
    threads: 4
    shell:
        """
        bwa mem -t {threads} {input.reference} {input.r1} {input.r2} | \
        samtools view -Sb - > {output}
        """

# remove pcr duplicates with samtools
rule dedup:
    input:
        "results/aligned/{sample}.bam"
    output:
        bam = "results/aligned/{sample}_processed.bam",
        metrics = "results/aligned/{sample}_processed_metrics.txt"
    threads: 4
    shell:
        """
        samtools sort -n -t {threads} {input} | \
        samtools fixmate -m - - | \
        samtools sort -t {threads} - | \
        samtools markdup -r -s -f {output.metrics} - {output.bam}
        """

# index processed bam file with samtools
rule index_bam:
    input:
        "results/aligned/{sample}_processed.bam"
    output:
        "results/aligned/{sample}_processed.bam.bai"
    shell:
        "samtools index {input}"

# run quality control on final processed alignment files
rule fastqc_processed:
    input:
        "results/aligned/{sample}_processed.bam"
    output:
        html = "results/fastqc/processed/{sample}_processed_fastqc.html",
        zip = "results/fastqc/processed/{sample}_processed_fastqc.zip"
    threads: 2
    shell:
        """
        fastqc -t {threads} {input} -o results/fastqc/processed/ && \
        mv results/fastqc/processed/{wildcards.sample}_processed_fastqc.html {output.html} && \
        mv results/fastqc/processed/{wildcards.sample}_processed_fastqc.zip {output.zip}
        """

# aggregate all statistics from raw baseline fastqc, final processed fastqc, and metrics
rule multiqc:
    input:
        raw = expand("results/fastqc/baseline/{sample}_{read}_baseline_fastqc.html", sample=samples, read=reads),
        processed_qc = expand("results/fastqc/processed/{sample}_processed_fastqc.html", sample=samples),
        metrics = expand("results/aligned/{sample}_processed_metrics.txt", sample=samples)
    output:
        "results/multiqc_report.html"
    shell:
        "multiqc results/ -o results/ -n multiqc_report.html --force"
