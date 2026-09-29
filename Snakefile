
# import required package
import os

configfile: "config/config.yaml"

sample = config["sample"]["accession"]
reference_name = config["sample"]["reference_name"]
num_threads = config["threads"]

rule all:
    input:
        f"results/multiqc/multiqc_report.html"

# rule 1: baseline fastqc
rule fastqc_raw:

    input:
        r1=f"data/{sample}_1.fastq",
        r2=f"data/{sample}_2.fastq"

    output:
        html1=f"results/qc/baseline/{sample}_1_fastqc.html",
        html2=f"results/qc/baseline/{sample}_2_fastqc.html"

    benchmark:
        f"benchmarks/fastqc_raw_{sample}.txt"

    log:
        f"logs/fastqc_raw_{sample}.log"

    threads: 2

    shell:
        """
        fastqc \
            -t {threads} \
            {input.r1} \
            {input.r2} \
            --outdir results/qc/baseline \
            > {log} 2>&1
        """

# rule 2: trim adapters & low quality bases
rule trim:

    input:
        r1=f"data/{sample}_1.fastq",
        r2=f"data/{sample}_2.fastq"

    output:
        r1=f"results/trimmed/{sample}_1_trimmed.fastq.gz",
        r2=f"results/trimmed/{sample}_2_trimmed.fastq.gz",
        u1=f"results/trimmed/{sample}_1_unpaired.fastq.gz",
        u2=f"results/trimmed/{sample}_2_unpaired.fastq.gz"

    params:
        adapters=os.path.join(
            os.environ["CONDA_PREFIX"],
            "share",
            "trimmomatic-0.40-0",
            "adapters",
            "TruSeq3-PE.fa"
        )

    benchmark:
        f"benchmarks/trimmomatic_{sample}.txt"

    log:
        f"logs/trimmomatic_{sample}.log"

    shell:
        """
        trimmomatic PE \
            {input.r1} \
            {input.r2} \
            {output.r1} \
            {output.u1} \
            {output.r2} \
            {output.u2} \
            ILLUMINACLIP:{params.adapters}:2:30:10 \
            LEADING:3 \
            TRAILING:3 \
            SLIDINGWINDOW:4:15 \
            MINLEN:36 \
            > {log} 2>&1
        """

# rule 3: remove phix contamination
rule bbduk:

    input:
        r1=f"results/trimmed/{sample}_1_trimmed.fastq.gz",
        r2=f"results/trimmed/{sample}_2_trimmed.fastq.gz"

    output:
        r1=f"results/trimmed/{sample}_1_decon.fastq.gz",
        r2=f"results/trimmed/{sample}_2_decon.fastq.gz"

    benchmark:
        f"benchmarks/bbduk_{sample}.txt"

    log:
        f"logs/bbduk_{sample}.log"

    shell:
        """
        bbduk.sh \
            in1={input.r1} \
            in2={input.r2} \
            out1={output.r1} \
            out2={output.r2} \
            ref=phix \
            k=31 \
            hdist=1 \
            > {log} 2>&1
        """

# rule 4: generate processed fastqc

rule fastqc_processed:

    input:
        r1=f"results/trimmed/{sample}_1_decon.fastq.gz",
        r2=f"results/trimmed/{sample}_2_decon.fastq.gz"

    output:
        html1=f"results/qc/processed/{sample}_1_fastqc.html",
        html2=f"results/qc/processed/{sample}_2_fastqc.html"

    benchmark:
        f"benchmarks/fastqc_processed_{sample}.txt"

    log:
        f"logs/fastqc_processed_{sample}.log"

    threads: 2

    shell:
        """
        fastqc \
            -t {threads} \
            {input.r1} \
            {input.r2} \
            --outdir results/qc/processed \
            > {log} 2>&1
        """

# rule 5: align reads to reference genome
rule align:

    input:
        r1=f"results/trimmed/{sample}_1_decon.fastq.gz",
        r2=f"results/trimmed/{sample}_2_decon.fastq.gz"

    output:
        bam=f"results/aligned/{sample}.bam"

    params:
        ref=f"references/{reference_name}.fa"

    benchmark:
        f"benchmarks/alignment_{sample}.txt"

    log:
        f"logs/alignment_{sample}.log"

    threads:
        num_threads

    shell:
        """
        bwa mem \
            -t {threads} \
            {params.ref} \
            {input.r1} \
            {input.r2} \
            2> {log} | \
        samtools view \
            -@ {threads} \
            -Sb - \
            > {output.bam}
        """

# rule 6: mark & remove pcr duplicates
rule dedup:

    input:
        bam=f"results/aligned/{sample}.bam"

    output:
        bam=f"results/aligned/{sample}_processed.bam",
        metrics=f"results/aligned/{sample}_markdup_metrics.txt"

    benchmark:
        f"benchmarks/dedup_{sample}.txt"

    log:
        f"logs/dedup_{sample}.log"

    threads:
        num_threads

    shell:
        """
        samtools sort \
            -n \
            -@ {threads} \
            {input.bam} | \
        samtools fixmate \
            -m \
            - \
            - | \
        samtools sort \
            -@ {threads} \
            - | \
        samtools markdup \
            -r \
            -s \
            -f {output.metrics} \
            - \
            {output.bam} \
            2> {log}
        """

# rule 7: build bam index
rule index_bam:

    input:
        bam=f"results/aligned/{sample}_processed.bam"

    output:
        bai=f"results/aligned/{sample}_processed.bam.bai"

    shell:
        """
        samtools index {input.bam}
        """

# rule 8: collect superficial mapping quality statistics
rule flagstat:

    input:
        bam=f"results/aligned/{sample}_processed.bam"

    output:
        f"results/aligned/{sample}_flagstat.txt"

    shell:
        """
        samtools flagstat {input.bam} > {output}
        """

# rule 9: collect comprehensive alignment stats
rule stats:

    input:
        bam=f"results/aligned/{sample}_processed.bam"

    output:
        f"results/aligned/{sample}_stats.txt"

    shell:
        """
        samtools stats {input.bam} > {output}
        """

# rule 10: combine all fastqc results into single dashboard
rule multiqc:

    input:
        f"results/aligned/{sample}_flagstat.txt",
        f"results/aligned/{sample}_stats.txt"

    output:
        f"results/multiqc/multiqc_report.html"

    log:
        "logs/multiqc.log"

    shell:
        """
        multiqc results \
            -o results/multiqc \
            --force \
            > {log} 2>&1
        """