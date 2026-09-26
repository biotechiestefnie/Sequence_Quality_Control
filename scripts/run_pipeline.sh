#!/bin/bash

# exit immediately if any command fails or variables unset
set -euo pipefail

# global portability configuration variables
SAMPLE="SRR25083113"
REF_URL="https://nih.gov"
# i think i should also have a variable called READS=[1, 2] or something like that here!

# start run timer
start_time=$SECONDS

# print update to terminal indicating initiation of execution
echo "Starting Bioinformatics Quality Control & Alignment"


# step 1: initialize and activate bioinfo environment with micromamba
echo "[INFO] Activating bioinfo environment..."
eval "$(micromamba shell hook --shell bash)"
micromamba activate bioinfo

# step 2: ensure core data directory exists within root
mkdir -p data

# step 3a: download sequencing data if both gzipped and fastq files missing
if [ ! -f "data/${SAMPLE}_1.fastq" ] && [ ! -f "data/${SAMPLE}_1.fastq.gz" ]; then
    echo "[INFO] Downloading forward sequencing reads..."
    wget -P data/ ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${SAMPLE:0:6}/013/${SAMPLE}/${SAMPLE}_1.fastq.gz

    echo "[INFO] Downloading reverse sequencing reads..."
    wget -P data/ ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${SAMPLE:0:6}/013/${SAMPLE}/${SAMPLE}_2.fastq.gz
else
    echo "[INFO] Sequencing reads or archive files already exist. Skipping download."
fi

# step 3b: extract sequencing data only if gzipped files exist but no fastq files
if [ -f "data/${SAMPLE}_1.fastq.gz" ] && [ ! -f "data/${SAMPLE}_1.fastq" ]; then
    echo "[INFO] Found compressed archives without text files. Unzipping sequencing files..."
    gunzip data/${SAMPLE}_1.fastq.gz
    gunzip data/${SAMPLE}_2.fastq.gz
else
    echo "[INFO] Uncompressed sequencing reads are already present. Skipping extraction."
fi

# step 4: download reference genome assembly matching variable template
if [ ! -f "data/${SAMPLE}_reference.fa" ]; then
    echo "[INFO] Downloading reference genome assembly for alignment template..."
    wget -qO- "${REF_URL}" | gunzip > data/${SAMPLE}_reference.fa
    echo "[INFO] Reference genome downloaded successfully."
else
    echo "[INFO] Matching reference genome template already exists. Skipping download."
fi

# step 5: index reference genome for bwa alignment if indices missing
if [ ! -f "data/${SAMPLE}_reference.fa.bwt" ]; then
    echo "[INFO] Building BWA reference genome index map layers..."
    bwa index data/${SAMPLE}_reference.fa
else
    echo "[INFO] BWA index map layers already exist. Skipping indexing."
fi

# step 6: run snakemake dry-run validation map
echo -e "\n[INFO] Running Snakemake validation dry-run..."
snakemake -n

# step 7: execute pipeline with 4 CPU cores
echo -e "\n[INFO] Executing live Snakemake workflow..."
snakemake --cores 4

# calculate total run time
end_time=$SECONDS
total_duration=$((end_time - start_time))

# convert seconds to minutes and seconds
minutes=$((total_duration / 60))
seconds=$((total_duration % 60))


echo "Pipeline Finished Successfully!"
echo "Master QC dashboard: results/multiqc_report.html"
echo "Mapped Genome BAM: results/aligned/${SAMPLE}_processed.bam"
echo "TOTAL RUNTIME: ${minutes}m ${seconds}s"

