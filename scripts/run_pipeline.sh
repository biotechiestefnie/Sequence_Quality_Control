#!/bin/bash

# stop execution on failure
set -euo pipefail

# start runtime counter to documentation at terminal
start_time=$SECONDS

# step 1: activate vitual env with necessary tools
# print environment activation update
echo "[info] activating bioinfo environment"
#
eval "$(micromamba shell hook --shell bash)"
# activate virtual env with needed tools
micromamba activate bioinfo
# check if bbmap present in env

# install bbmap into bioinfo env for decon

# create runtime directories
mkdir -p data
mkdir -p references

# parse config values
sample=$(grep "accession:" config/config.yaml | awk '{print $2}')
reference_name=$(grep "reference_name:" config/config.yaml | awk '{print $2}')
reference_url=$(grep "reference_url:" config/config.yaml | awk '{print $2}')
threads=$(grep "threads:" config/config.yaml | awk '{print $2}')

# print run id to terminal
echo "[info] sample: ${sample}"

# step 2: retrieve & unpack sample files
# check that fwd fastq not already present
if [ ! -f "data/${sample}_1.fastq" ]; then
    # print update for downloading paired-end run files
    echo "[info] downloading paired-end reads"
    # retrieve fwd read file
    wget -P data \
    ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${sample:0:6}/013/${sample}/${sample}_1.fastq.gz
    # retrieve rvs read file
    wget -P data \
    ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${sample:0:6}/013/${sample}/${sample}_2.fastq.gz

    # unpack gzipped fastq fwd/rvs paired end run files
    gunzip data/${sample}_1.fastq.gz
    gunzip data/${sample}_2.fastq.gz

fi

# step 3: retrieve corresponding reference genome
# ensure reference not already present
if [ ! -f "references/${reference_name}.fa" ]; then
    # print reference extraction update to terminal
    echo "[info] downloading reference genome"
    # pull reference genome
    wget -qO- "${reference_url}" | gunzip > references/${reference_name}.fa

fi

# step 4: generate index for reference genome
# check that index doesnt already exist
if [ ! -f "references/${reference_name}.fa.bwt" ]; then
    # print bwa index update to terminal
    echo "[info] building bwa index"
    # build bwa index
    bwa index references/${reference_name}.fa

fi

# step 5: perform dry run to check file system before execution
# print dry run update to terminal
echo "[info] Running Snakemake Dry Run"
# validate snakemake workflow with dry run
snakemake -n

# step 6: run project workflow
# print execution initiation update to threshold
echo "[info] Executing Workflow"
# execute quality control/alignment workflow
snakemake --cores ${threads}

# calculate workflow runtime for output to terminal
duration=$((SECONDS-start_time))

minutes=$((duration/60))
seconds=$((duration%60))

# print workflow completion status update to terminal
echo "Pipeline Completed Successfully"
echo "Runtime: ${minutes}m ${seconds}s"  # print workflow runtime
echo "Multiqc Report: results/multiqc/multiqc_report.html"  # print dashboard location
