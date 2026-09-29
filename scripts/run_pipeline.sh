#!/bin/bash

# stop execution on failure
set -euo pipefail

# start runtime counter
start_time=$SECONDS

# print status to terminal
echo "[info] Activating bioinfo Environment"
# activate environment with micromamba
eval "$(micromamba shell hook --shell bash)"
micromamba activate bioinfo

# create runtime directories
mkdir -p data
mkdir -p references
mkdir -p logs
mkdir -p benchmarks

# parse config values
sample=$(grep "accession:" config/config.yaml | awk '{print $2}')
reference_name=$(grep "reference_name:" config/config.yaml | awk '{print $2}')
reference_url=$(grep "reference_url:" config/config.yaml | awk '{print $2}')
threads=$(grep "threads:" config/config.yaml | awk '{print $2}')

# print sample name to terminal
echo "[info] Sample: ${sample}"

# download reads rule
# check if unzipped read files already present
if [ ! -f "data/${sample}_1.fastq" ]; then
    # print update to terminal
    echo "[info] Downloading Forward Reads"
    # extract gzipped fwd reads file
    wget -P data \
    ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${sample:0:6}/013/${sample}/${sample}_1.fastq.gz
    # print update to terminal
    echo "[info] Downloading Reverse Reads"
    # extract gzipped rvs reads file
    wget -P data \
    ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${sample:0:6}/013/${sample}/${sample}_2.fastq.gz
    # unzip fwd & rvs read files
    gunzip data/${sample}_1.fastq.gz
    gunzip data/${sample}_2.fastq.gz

else
    # print update to terminal
    echo "[info] Skipping: Reads Already Present"

fi

# download reference genome rule
# check if reference genome already present
if [ ! -f "references/${reference_name}.fa" ]; then
    # print update to terminal
    echo "[info] Downloading Reference Genome"
    # retrieve reference genome
    wget -qO- "${reference_url}" | gunzip > references/${reference_name}.fa

else
    # print reference genome already present to terminal
    echo "[info] Reference Already Present"

fi

# bwa index for reference genome rule
# check that index not present
if [ ! -f "references/${reference_name}.fa.bwt" ]; then
    # print index building status to terminal
    echo "[info] Building BWA Index"
    # build bwa index for reference
    bwa index references/${reference_name}.fa

else
    # print already present to terminal
    echo "[info] Skipping: BWA Index Already Present"

fi

# print dry run execution status to terminal
echo "[info] Running Dry Run"
# validate workflow
snakemake -n

# print execution status update to terminal
echo "[info] Executing Workflow"
# execute workflow
snakemake \
    --cores ${threads} \
    --printshellcmds

# calculate runtime
duration=$((SECONDS-start_time))
# convert to min/sec
minutes=$((duration/60))
seconds=$((duration%60))

# print completion status update to terminal
echo "Pipeline Completed Successfully"
echo "Runtime: ${minutes}m ${seconds}s"
echo "Multiqc Report: results/multiqc/multiqc_report.html"
