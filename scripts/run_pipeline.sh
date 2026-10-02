#!/bin/bash

# stop execution on failure or unspecified variables
set -euo pipefail

# start runtime counter
start_time=$SECONDS

# create project subdirectories if not present
mkdir -p "raw"
mkdir -p "resources/ref_genome"
mkdir -p "resources/adapters"
mkdir -p "logs"
mkdir -p "benchmarks"

# configuration lookups for params
sample=$(python3 -c "import yaml; print(yaml.safe_load(open('config/config.yaml'))['sample']['accession'])")
genome_name=$(python3 -c "import yaml; print(yaml.safe_load(open('config/config.yaml'))['sample']['genome_name'])")
genome_url=$(python3 -c "import yaml; print(yaml.safe_load(open('config/config.yaml'))['sample']['genome_url'])")
adapter_type=$(python3 -c "import yaml; print(yaml.safe_load(open('config/config.yaml'))['trimming']['adapter_type'])")
adapter_url=$(python3 -c "import yaml; print(yaml.safe_load(open('config/config.yaml'))['trimming']['adapter_url'])")
threads=$(python3 -c "import yaml; print(yaml.safe_load(open('config/config.yaml'))['threads'])")

# specify run accession
echo "[info] Target Dataset: ${sample}"

# retrieve paired-end sequencing reads
if [ -f "raw/${sample}_1.fastq.gz" ] && [ -f "raw/${sample}_2.fastq.gz" ]; then
    echo "[info] Skipping: Compressed FASTQ Files Already Present"

elif [ -f "raw/${sample}_1.fastq" ] && [ -f "raw/${sample}_2.fastq" ]; then
    echo "[info] Compressing Existing FASTQ Files..."
    gzip "raw/${sample}_1.fastq"
    gzip "raw/${sample}_2.fastq"

else
    echo "[info] Downloading Forward Reads..."
    wget -qP "raw" \
        "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${sample:0:6}/013/${sample}/${sample}_1.fastq.gz"

    echo "[info] Downloading Reverse Reads..."
    wget -qP "raw" \
        "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${sample:0:6}/013/${sample}/${sample}_2.fastq.gz"
fi

# retrieve reference genome for alignment
if [ -f "resources/ref_genome/${genome_name}.fa" ]; then
    echo "[info] Skipping: Reference Genome Already Present"

elif [ -f "resources/ref_genome/${genome_name}.fa.gz" ]; then
    echo "[info] Decompressing Existing Reference Genome..."
    gunzip -k "resources/ref_genome/${genome_name}.fa.gz"

else
    echo "[info] Downloading Reference Genome..."
    wget -qO- "${genome_url}" | gunzip \
        > "resources/ref_genome/${genome_name}.fa"
fi

# download sequence adapters according to NGS run for trimming
if [ ! -f "resources/adapters/${adapter_type}.fa" ]; then
    echo "[info] Downloading ${adapter_type} Adapter Fasta..."
    wget -qO "resources/adapters/${adapter_type}.fa" "${adapter_url}"

else
    echo "[info] Adapter Fasta File Already Present"
fi

# build bwa index for ref genome if absent
if [ ! -f "resources/ref_genome/${genome_name}.fa.bwt" ]; then
    echo "[info] Building BWA Index..."
    bwa index "resources/ref_genome/${genome_name}.fa" 2>/dev/null

else
    echo "[info] Skipping: BWA Index Already Present"
fi

# dry run to confirm workflow integrity
echo "[info] Running Validation Dry Run..."
snakemake -n --rerun-incomplete

# execute qc & alignment workflow
echo "[info] Executing Live Multi-Core Workflow..."
snakemake \
    --cores "${threads}" \
    --printshellcmds \
    --rerun-incomplete

# calculate total runtime
duration=$((SECONDS-start_time))
minutes=$((duration/60))
seconds=$((duration%60))

# completion status, total runtime, dashboard location to terminal
echo "Pipeline Completed Successfully!"
echo "TOTAL WORKFLOW RUNTIME: ${minutes}m ${seconds}s"
echo "MultiQC Dashboard Located In 'deliverables/' Folder/"