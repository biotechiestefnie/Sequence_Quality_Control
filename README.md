#  Quality Control Assessments of Raw Versus Processed Paired End WGS Run SRR25083113 for *E. coli* 

This project contains an automated and portable bioinformatics pipeline for WGS quality
control utilizing Snakemake for workflow management. It processes raw paired-end 
next-generation sequencing (NGS) genomic DNA-seq data from *Escherichia coli* (SRR25083113) 
using the following de facto standard tools: FastQC, Trimmomatic, Samtools, BBmap, BWA, 
and MultiQC.  

## Pipeline Steps:

**Shell Script:**

1. Activate Virtual Environment: Command Micromamba to activate the bioinfo environment, 
   which contains the bioinformatics tools needed for this project.
2. Extract Sequence Files and Reference Genome: Download the *E. coli* O157:H7 bacterial 
   genome (SRR25083113) stored as two FastQ files (forward and reverse paired end reads) 
   from ENA and the K-12 MG1655 reference genome from the NCBI RefSeq FTP server to the local 
   system.
3. Index Reference Genome: Construct positional coordinate memory index maps of the reference
   genome sequence utilizing BWA to provide a rapid lookup template that allows the short-read 
   aligner to efficiently map and resolve millions of sequencing fragments against the 
   chromosome loci.

**Snakemake Script:**

1. Assess Raw Reads: Perform a quality control initial assessment with FastQC to identify base-
   calling drops, adapter contamination, or tile anomalies for baseline comparisons.
2. Trim & Filter: Purge poor-quality sequences, cut adapter read-throughs, and discard reads 
   that drop below the specified 36 bp length threshold with Trimmomatic.
3. Map Reads: Coordinate clean fragments to the designated chromosome loci with BWA.
4. Remove PCR Duplicates: Find and scrub co-localized read amplification duplicates to normalize
   dynamic alignment metrics using Samtools.
5. Index Aligned BAM File: Generate a coordinate-sorted tracking index map of the final 
   deduplicated sequence alignment file utilizing Samtools, creating a downstream .bai 
   companion file that allows variant-calling software and genome visualization browsers like 
   IGV to perform rapid coordinate jumps and instantly retrieve specific read regions without 
   searching the entire binary file.
6. Assess Final Alignment: Evaluate post-deduplication data health within the binary BAM layout
   using FastQC again, producing plots to compare against baseline run on raw reads.
7. Aggregate Dashboards: Gather the cross-workflow QC logs into a centralized file layout for 
   seamless reporting using MultiQC.

---

## Directory Structure

The workspace must be organized according to the structure below for the pipeline to resolve 
the file patterns correctly:

```text
Sequence_Quality_Control/
│
├── config/
│   └── config.yaml
│
├── data/
│   ├──
│   ├──
│   ├──
│   ├──
│   ├──
│   ├──
│   └──
│
├── references/
│   ├── ecoli_k12.fa
│   ├── ecoli_k12.fa.amb
│   ├── ecoli_k12.fa.ann
│   ├── ecoli_k12.fa.bwt
│   ├── ecoli_k12.fa.pac
│   └── ecoli_k12.fa.sa
│
├── envs/
│   └── bioinfo.yaml
│
├── scripts/
│   └── run_pipeline.sh
│
├── results/
│   ├── aligned/
│   ├── benchmarks/
│   ├── multiqc/
│   ├── qc/
│   │   ├── baseline/
│   │   └── processed/
│   └── trimmed/
│
├── benchmarks/
│
├── Snakefile
│
├── README.md
│
├── workflow.png
│
├── final_report.html
│
└── .gitignore

```

---

## Methodology & Bioinformatics Tools:

The Snakemake pipeline orchestrates five industry-standard bioinformatics engines:

### 1. FastQC (FastQ and BAM sequence quality assessments)

* Purpose: Baseline analysis of raw sequence quality and processed, aligned quality assessment
  for comparison.

* Application: Run at two distinct points within the Snakemake workflow. The first pass scans 
  the raw FastQ data files to establish baselines of raw Phred scores, GC content skew, and 
  adapter contamination status. The second pass assesses the quality of the final sequences
  located within the BAM file that have undergone trimming, filtering, alignment to the 
  reference genome, deduplication, and indexing to extract sequence lines and confirm data 
  quality improvements after processing.

* Background: FastQC is an industry-standard, open-source bioinformatics software designed 
  to perform quality control (QC) checks on high-throughput sequencing data. Developed by 
  Simon Andrews at Babraham Bioinformatics in Java, it provides a quick, modular overview of 
  sequencing runs, allowing researchers to catch technical issues, biases, or contamination 
  before spending hours or days on downstream genomic analysis. The program accepts the 
  following file formats: FastQ, SAM, and BAM files (any variant). Detailed documentation 
  for FastQC can be found at: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/

### Trimmomatic (PE - Quality Trimming)

* Purpose: Initial read cleaning, adapter stripping, and low-quality base clipping for artifact
  removal.

* Application: Run as the primary processing filter immediately after raw quality control 
  checks. The tool takes the raw forward and reverse paired-end FastQ files from the data 
  folder, scans them using a sliding-window approach, and evaluates them based on the 
  Snakemake script's established thresholds. It completely drops low-quality read endings, 
  crops remaining adapter read-through segments, and filters out short, uninformative sequences
  to generate clean, paired survival fragments optimized for high-performance genome mapping.

* Background: Trimmomatic is a fast, multithreaded command-line software application developed
  by Anthony Bolger and Björn Usadel at the Usadel Lab (Aachen University). Written natively 
  in Java, it can process both single-end and paired-end data, but it was specifically designed
  to handle paired-end sequencing datasets produced by Illumina platforms. When processing 
  paired-end reads, it evaluates sequences in a synchronized fashion to maintain the exact 
  sequence correspondence of read pairs and avoid tracking errors. Prior to Trimmomatic, when
  one read in a pair passed the quality thresholds but its partner failed, trimmers could no 
  longer track them together. To mitigate this issue, Trimmomatic splits the data into four 
  distinct output files: two synchronized files containing the surviving forward and reverse 
  read pairs, and two separate files containing the discarded, unpaired forward and reverse 
  reads whose corresponding partners failed. Detailed documentation for Trimmomatic can be 
  found at: http://www.usadellab.org/cms/?page=trimmomatic.

### BWA (Burrows-Wheeler Aligner)

* Purpose: High-performance genomic mapping and short-read sequence alignment tracking against
  a reference coordinate template.

* Application: Run immediately after Trimmomatic trimming and filtering stage. The tool 
  accesses the surviving, fully synchronized forward and reverse read pairs 
  (r1_paired.fastq.gz and r2_paired.fastq.gz), mapping them with the BWA MEM alignment 
  algorithm against the indexed *Escherichia coli* K-12 reference genome (reference.fa). 
  This algorithm is the newest and most recommended algorithm for high-quality queries; 
  it performs faster and more accurate seed-and-extend alignments. BWA MEM was designed 
  to handle longer Illumina reads from 70 bp up to 1 Mbp and supports split alignments
  for structural variations. It finds the exact position where each sequence fragment 
  belongs on the bacterial chromosome and streams those structural coordinates directly 
  into Samtools to build the initial binary alignment map (.bam) file.

* Background: The Burrows-Wheeler Aligner (BWA) is an industry-standard, open-source software 
  suite designed for aligning low-divergence DNA and RNA sequencing reads against a large 
  reference genome. Developed by Heng Li and Richard Durbin at the Wellcome Trust Sanger 
  Institute, it is written in C and applies the Burrows-Wheeler Transform to optimize 
  processing speeds and minimize system memory requirements while efficiently aligning 
  reads and handling mismatches or gaps. Although in theory BWA works with arbitrarily 
  long reads, its performance is degraded on long reads, especially when the sequencing 
  error rate is high. Furthermore, BWA always requires the full read to be aligned, but 
  longer reads are more likely to be interrupted by structural variations or mis-assemblies 
  in the reference genome, which will fail BWA. Therefore, if running the program with long 
  reads, it is important to ensure the quality of the reference genome is optimized, or 
  consider dividing long reads into multiple shorter fragments and aligning them separately 
  with BWA, then joining the partial alignments to achieve the complete alignment of long 
  reads to a suboptimal quality reference genome. Detailed documentation for BWA can be 
  found at: https://bio-bwa.sourceforge.net/bwa.shtml.

### Samtools (Sorting, Coordinate Processing, and Deduplication)

* Purpose: High-throughput stream processing, file formatting, coordinate sorting, and 
  structural artifact removal for sequence alignment files.

* Application: Applied as a multi-stage execution pipeline immediately following the bwa mem 
  mapping step. The program takes the initial uncompressed alignment stream as an input and 
  utilizes Samtools view with the -Sb parameters to compress the human-readable data into a 
  binary alignment map (BAM) file. It then processes this file sequentially through 
  name-sorting, mate-pair fixing, coordinate-sorting, and deduplication. The parameters 
  used in this project filter out PCR duplicate reads, reduce coverage depth biases caused 
  by amplification artifacts, and index the final clean file. The tool outputs the final, 
  deduplicated alignment map, a companion index map, and a text summary of duplication 
  statistics into the results folder within the project directory.

* Background: Samtools is a foundational, open-source software suite designed specifically for
  interacting with and manipulating high-throughput sequencing alignments stored in SAM, BAM, 
  or CRAM file formats. Originally developed by Heng Li at the Wellcome Trust Sanger Institute 
  and maintained by the global HTSlib community, it is engineered in C to process massive 
  genomics datasets as a high-speed data stream using standard Unix piping structures. In 
  paired-end workflows, a name-sorting step is required so that Samtools fixmate can read 
  matching read names consecutively, verifying mate-pair orientations, template lengths, and 
  positional score tags. The data must then be re-sorted by genomic coordinates to structure 
  the file for positional matching. The Samtools markdup algorithm utilizes specific mate-
  score (ms) and mate-CIGAR (MC) tags appended by fixmate to compare read locations and 
  orientation; if multiple fragments match identical coordinates, Samtools preserves the 
  read with the highest physical base quality scores and completely removes or flags the 
  remaining copies as analytical duplicates to protect downstream variant-calling steps 
  from frequency inflation. Detailed documentation for Samtools can be found at: 
  https://htslib.org.

### MultiQC (Workflow Reporting Aggregation)

* Purpose: High-throughput reporting compilation, metadata parsing, and cross-workflow 
  quality control aggregation.

* Application: Executed as the ultimate master diagnostic stage of the pipeline, immediately
  following post-processing FastQC sequence quality assessment. The program targets the root 
  results/ folder as its primary input, automatically scanning all deep subdirectories for 
  raw console outputs, tool telemetry data, and tracking text documents. It parses the 
  separate text arrays generated by the project's raw FastQC rules inside the Snakemake file, 
  the Trimmomatic execution logs, the Samtools markdup metric tracking files, and the final 
  post-processed BAM FastQC assets. The tool utilizes the --force flag to clear out outdated 
  cache dependencies and compresses this cross-workflow metrics footprint into a single, 
  interactive dashboard webpage file located inside the results folder.

* Background: MultiQC is an open-source, Python-based reporting application designed to 
  resolve the diagnostic bottleneck associated with analyzing multi-stage high-throughput 
  sequencing datasets. Developed by Philip Ewels at SciLifeLab (Stockholm) and maintained 
  natively alongside the Seqera Platform group, the software does not perform standalone 
  quality control assessments. Instead, it searches specified target directories recursively 
  to parse distinct output streams, standard errors, and log syntax matrices belonging to 
  over 150 supported bioinformatics tools. By normalizing and aggregating disparate run 
  properties such as read truncation values from Trimmomatic, sequence duplication curves 
  from FastQC, and percentage mapping depths from Samtools into a single web browser layout, 
  it allows researchers to instantly isolate project trends, catch batch anomalies, and 
  evaluate processing improvements at a single glance. Detailed documentation for MultiQC 
  can be found at the Official Seqera MultiQC Portal: https://seqera.io/multiqc/.
---

## Workflow

![Workflow](workflow.png)

---

## Instructions for Running this Project

### Prerequisites:

While the project was designed to be executed with Rosetta on a 32Gb Apple Silicon M2 Max 
running macOS Golden Gate Version 27.0, arm64 configurations of all the bioinformatics tools 
in this project can be downloaded specifically configured to be executed without requiring 
Rosetta, and this will be necessary for all future macOS version releases. However, this 
approach requires each tool to be downloaded separately, as the `bioinfo` tool bundle was only 
designed for Intel-based configurations (osx-64).

For those intending to follow my precise workflow, Micromamba is implemented to activate the 
`bioinfo` virtual environment and tool bundle, sourced from the Biostar Handbook.This 
environment was chosen for its centralized bioinformatics tool accessibility, as it natively 
contains all the programs utilized in this project.

For those forgoing `bioinfo` and planning on downloading each tool separately, any conda
application can be used to activate a virtual environment for running the project. It will
then be necessary to replace the environment activation command on line 12 of the Bash script
with the following, which creates and activates a new custom workspace, then downloads all
required tools into this workspace and activates it so the project will run seamlessly:

```bash
# replace line 12 of run_pipeline.sh with this command block to create your own virtual
# environment, download the tools required for this project, and activate the new environment
conda create -n custom_bioinfo \
  -c conda-forge \
  -c bioconda \
  python=3.11 snakemake fastqc trimmomatic bwa samtools multiqc -y && \
conda activate custom_bioinfo

```

For Intel-based Mac systems, Linux systems, or Windows Subsystem for Linux (WSL), these 
configuration translation layers are not applicable. Only Micromamba and `bioinfo` will be 
required to run the project on your system without making any script alterations. 
Alternatively, if you choose not to run Micromamba and `bioinfo`, the standalone conda 
setup instructions above can be applied without modification, utilizing a conda application 
of your preference to create a custom environment with all the tools required for this project.

### Step 1: Run the Master Automation Script

This script activates the `bioinfo` virtual environment, downloads and unzips the FastQ.gz forward 
and reverse read files, retrieves the reference genome and builds the bwa index, and launches 
the Snakemake file. The command to execute the master automation script is called from the
root directory in the local terminal as follows:

```bash
./run_pipeline.sh
```

### Step 2: Evaluate Results

Click on the generated output documents located inside the results folder, which will launch
inside your browser. Inspect the sequencing metrics and coverage parameters:

* Master Interactive Dashboard: Open `results/multiqc_report.html`

* Final Processed Genomic Mappings: File located at `results/aligned/srr25083113_dedup.bam`


## Portability: Application for Different SRR Datasets

If you wish to run this workflow with a different paired-end run and corresponding reference 
genome, simply follow the modifications below:

### Step #1: Swap

Drop your new paired-end raw FastQ.gz files into the `data/` directory following the `_1.fastq.gz` and `_2.fastq.gz` naming convention, then update the samples variable at the top of the `Snakefile` script to include your new sample ID string within the array brackets:

```python
samples = ["YOUR_SAMPLE_SRR"] # swap your unique dataset run id here
reads = ["1", "2"]
```

**Step #2:** If your new dataset maps to a different reference organism, look inside the `run_pipeline.sh` master shell script at **Line 47** (under Step 4). Replace the NCBI FTP URL inside the `curl -L` command with the download link corresponding to your new reference genome template:

```bash
# Line 47: Swap out the URL inside the quotes for your new assembly reference
curl -L "https://NEW_NCBI_FTP_REFERENCE_URL.fna.gz" | gunzip > data/reference.fa
```

**IMPORTANT!!!:** This specific pipeline configuration is built strictly for paired-end sequencing datasets and cannot execute a single-read dataset without manual alterations to the core Trimmomatic and BWA rule frameworks inside the `Snakefile`.


---

## Troubleshooting

### Stuck Folder Locks (`LockException`)
If the execution is manually canceled using `Ctrl + C` or interrupted by a system crash, 
Snakemake will leave a hidden safety lock on your working folder to prevent data corruption. 
Subsequent attempts to run the pipeline will result in a `LockException` error message.

To resolve this and release the folder lock, execute the following command in your terminal 
before launching the master shell script again:

```bash
snakemake --unlock
```

to change species later, place : sample: SRRXXXXXXX

reference_name: hg38

reference_url: https://...

into config.yaml, without touching the script

In the deduplication step, we introduced the tool called dedupe.sh from BBMapLinks to an external site.. Again, there are a number of tools out there for decontamination, but we will stick with BBMap right now and use a tool called bbduk.sh.

Similar to dedupe.sh, the documentation is all over the place in bbduk.sh. Sometimes there isn't enough information. Othertimes, it feels like there is too much information. Feel free to use bbduk.sh --help if you want, but I suggest using the US DoE Joint Genome Institute guideLinks to an external site..

BBDuk, part of bbmap, is a versatile, high-performance tool for cleaning up sequencing data—chiefly by removing contaminants using kmer-based strategies. The heart of BBDuk's contaminant removal method is the use of kmers, the same fixed-length substrings that are central to genome assembly techniques. This parallel is crucial: just as kmers connect sequencing reads to reconstruct entire genomes, they also allow BBDuk to compare and match reads against known contaminant sequences with precision and speed.

Kmers: Fundamental to Assembly and Filtering
In genome assembly, short reads are split into kmers (e.g., 31 bases long), enabling the assembler to build graphs that depict how reads overlap by their shared subsequences. By linking reads through shared kmers, assemblers reconstruct long contiguous genomic regions from short fragments. BBDuk harnesses this same technology—but instead of connecting overlaps, it checks whether any read shares a kmer with entries in a separate reference of known contaminants. This simple, fast matching process effectively identifies and removes unwanted sequences.

How BBDuk Removes Contaminants
A common step in some library prep is the use of "spike-ins." A spike-in is a technical control library added into sequencing libraries to compensate for low base diversity, support the clustering algorithm, and control for sequencing accuracy. One of the most common Illumina spike-ins is PhiX Control Library v3 (known as "PhiX").

These controls are helpful for a number of reasons and are part of normal sequencing prep. However, the bioinformaticist likely has no control over their use. One of the negatives of using spike-ins is that the artificial sequence that is added is not the sequence of the sample. If left unchecked, these sequences can contaminate sequencing data.
Reference Preparation: BBDuk requires a FASTA file containing the sequences to remove, such as host DNA, environmental contaminants, adapters, or spike-ins.
Kmer-Based Filtering: Both the sequencing reads and reference contaminants are broken into kmers of chosen length 
 (often 23–31). BBDuk slides a window of size 
 over each sequence, building a kmer library for every contaminant sequence. When scanning the reads, it compares the read’s kmers to the contaminant kmers, allowing users to specify the number of allowed mismatches with the hdist parameter. This is conceptually identical to how assembly algorithms match kmer overlaps between reads, but in BBDuk the objective is to detect sequence identity rather than continuity.
Filtering Logic and Modes:
If a read contains a match to a reference kmer, BBDuk can:
Remove the entire read from the output (default behavior)
Output matching reads to a separate file for inspection
Mask the matching segment within the read, typically replacing it with an 'N'
Customization and Sensitivity Tuning: Users control parameters such as kmer size (k), mismatch allowance (hdist for Hamming distance or edist for edit distance), and the minimal number of matching kmers required for a read to be considered a contaminant (mkh).
Paired-End Handling: BBDuk keeps or discards paired reads together, ensuring consistency for downstream analyses.
Additional Features
Multiple Trimming/Filtering Options: BBDuk can also perform adapter trimming, quality filtering, length and entropy filtering, sequence masking, barcode filtering, and more along with contaminant removal.
Memory and Performance: Efficiently uses available memory and threads, scaling to large datasets and complex references.
Reporting: Generates statistics outlining how many reads matched each contaminant sequence, supporting thorough quality assessments.
Tuning and Special Considerations
Specificity vs. Sensitivity: Adjusting kmer length and mismatch tolerance tunes the trade-off between missing mutated contaminants (with long kmers or few mismatches) and risking false positives (with short kmers or more mismatches).
Complex Reference Management: The number of unique kmers in large references (like whole genomes) affects memory usage linearly and processing time exponentially, especially when mismatches are allowed.
Multiple Filtering Strategies: While BBDuk supports a single kmer-based operation per run, BBDuk2 allows multiple such operations concurrently.
Parameter Fine-Tuning: Options such as mink (minimal kmer size for flexible detection), maskmiddle (masking the middle of kmers to enhance sensitivity), and rcomp (searching both forward and reverse complement sequences) provide further flexibility for diverse biological scenarios.
BBDuk’s contaminant filtering exemplifies the power of kmers, binding sequencing quality control with core assembly principles. By leveraging the rapid, graph-based search methods established in assembly, BBDuk delivers accurate, scalable, and customizable contaminant removal—paving the way for cleaner and more reliable sequencing data.
To remove PhiX (or any other) contaminants, you need a reference to read from. Thankfully, since this is such a common task, the author of BBMap has prepackaged bbduk.sh it with multiple well-known contaminant reference files. Additionally, you don't even need to use a random path to access them, like we did for trimmomatic. We can just use the phix keyword.

Another quality of life (QoL) feature of bbduk.sh is that if it is given an interleaved fastq file, it will automatically detect it and handle it appropriately. Kind of nice because it saves some typing.


DEDUPLICATING
A duplicate read is usually an artifact of the PCR process. During PCR, the same DNA fragment can be copied (or duplicated) more than once. The premise of duplicate removal is predicated on the belief that some of these reads may contain mutations caused by the PCR process or inordinately affect the occurrence of duplicated alleles (Ebbert et al., 2016
 It is actually a regular practice when working with RNA to never remove duplicates because it is ultimately viewed as a loss of "good" data. Recent research suggests that, given the way we sequence things now, duplicate removal has very little impact on the results of the research.

That said, "deduping" is a common practice in many workflows and pipelines. Because of that, we will introduce the process.


TRIMMING ADAPTERS AND LOW QUALITY BASE CALLS

It cannot be said enough how integral FastQC is to the genomics, transcriptomics, and epigenomics communities. While there are a handful of other tools out there for evaluating datasets for QC purposes
FastQC is somewhat the de facto standard for most workflows
For the most part, one of the first preprocessing steps in sequencing is to "trim" your reads. That is, remove the special adapters used by the sequencer for its isolation, amplification, priming, etc. In many short-read sequencing technologies, reads that are shorter than expected are sequenced along with the "correct length" reads. This can often cause the sequencer to "read through" the read and add the attached adapter in with the sequencing (source: https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-1069-7
Another use case of trimming is if a researcher wants to enrich their data only with high-quality data. They can do this by trimming away "bad" data from their reads. For example, a FastQC result with a "Quality scores across all bases" plot that looks like the one below suggests that there is a lot of "bad" data in the tail ends of the reads.
This image came from a tutorial on Cleaning and Filtering Reads by Data Carpentry at https://datacarpentry.github.io/wrangling-genomics/02-quality-control.htmlLinks to an external site..

This can tend to happen during sequencing for a lot of reasons, but for the most part, quality scores tend to degrade as the sequencing run progresses. Nevertheless, researchers can trim the "bad" data off the reads.
trimmomatic is super flexible in its ability to trim reads. It has settings for just about everything and comes prepackaged with known adapters for multiple sequencing technologies.

The biggest issue with trimmomatic is that it is a Java program. This means it doesn't act like a "normal" CLI program that you are probably used to using. Additionally, the authors used non-standard argument parsing conventions, so pay close attention to the command templates that you will be provided with.
The third one is "tricky" in that it isn't super intuitive. trimmomatic needs a list of known adapters, and since trimmomatic comes prepackaged with them, we have to look in the trimmomatic directory for them. Furthermore, since trimmomatic doesn't work like a "normal" program, we need to know where it is as well. Since we don't "know" where that directory is, we will leverage the module command for show us:
TRIM_ADAPT=$(dirname "$(which trimmomatic-0.39.jar)")


FASTQC

For each FASTQ file that is input to FastQC, two output files are generated.

The first is an HTML file, a self-contained document with various graphs embedded within it. Each of the graphs evaluates different aspects of our data quality, which we will discuss in more detail in this lesson.
Alongside the HTML file is a zip file (with the same name as the HTML file, but with .zip added to the end). This file contains the different plots from the report as separate image files, but also contains data files that are designed to be easily parsed to allow for a more detailed and automated evaluation of the raw data on which the QC report is built.
FastQC has a really well-documented manual pageLinks to an external site. with detailed explanationsLinks to an external site. about every plot in the report.

Within our report, a summary of all of the modules is given on the left-hand side of the report. Don’t take the yellow “WARNING”s and red “FAIL”s too seriously; they should be interpreted as flags for modules to check out.
The first module gives the basic statistics for the sample. Generally, it is a good idea to keep track of the total number of reads sequenced for each sample and to make sure the read length and %GC content are as expected.
Per base sequence quality
One of the most important analysis modules is the “Per base sequence quality” plot. This plot provides the distribution of quality scores at each position in the read across all reads. The y-axis gives the quality scores, while the x-axis represents the position in the read. The color coding of the plot denotes what are considered high, medium, and low quality scores.

This plot can alert us to whether there were any problems occurring during sequencing and whether we might need to contact the sequencing facility.
For example, the box plot at nucleotide 1 shows the distribution of quality scores for the first nucleotide of all reads in the Mov10_oe_1 sample. The yellow box represents the 25th and 75th percentiles, with the red line as the median. The whiskers are the 10th and 90th percentiles. The blue line represents the average quality score for the nucleotide. Based on these metrics, the quality scores for the first nucleotide are quite high, with nearly all reads having scores above 28.

The quality scores appear to drop going from the beginning toward the end of the reads. For reads generated by Illumina sequencing, this is not alarming, and there are known causes for this drop in quality.

For Illumina sequencing, the quality of the nucleotide base calls is related to the signal intensity and purity of the fluorescent signal. Low intensity fluorescence or the presence of multiple different fluorescent signals can lead to a drop in the quality score assigned to the nucleotide. Due to the nature of sequencing-by-synthesis, there are some drops in quality that can be expected, but other quality issues can be indicative of a problem at the sequencing facility.

We will now explore different quality issues arising from the sequencing-by-synthesis used by Illumina, both expected and unexpected.

Expected Error Profiles
As sequencing progresses from the first cycle to the last cycle, we often anticipate a drop in the quality of the base calls. This is often due to signal decay and phasing as the sequencing run progresses.

Signal decay: As sequencing proceeds, the fluorescent signal intensity decays with each cycle, yielding decreasing quality scores at the 3’ end of the read. This is due to:

Degrading fluorophores
A proportion of the strands in the cluster are not being elongated
Therefore, the proportion of signal being emitted continues to decrease with each cycle.
Phasing: As the number of cycles increases, the signal starts to blur as the cluster loses synchronicity, also yielding a decrease in quality scores at the 3’ end of the read. As the cycles progress, some strands experience random failure of nucleotides to incorporate due to:

Incomplete removal of the 3’ terminators and fluorophores
Incorporation of nucleotides without effective 3’ terminators
Worrisome Error Profiles
Overclustering: Sequencing facilities can overcluster the flow cells, which results in small distances between clusters and an overlap in the signals. The two clusters can be interpreted as a single cluster with mixed fluorescent signals being detected, decreasing signal purity, and generating lower quality scores across the entire read.
Instrumentation breakdown: Sequencing facilities can occasionally have issues with the sequencing instruments during a run. Any sudden drop in quality or a large percentage of low-quality reads across the read could indicate a problem at the facility. Examples of such issues are shown below, including a manifold burst, cycles lost, and a read 2 failure. For such data, the sequencing facility should be contacted for resolution, if possible.
The “Per sequence quality scores” plot gives you the average quality score on the x-axis and the number of sequences with that average on the y-axis. We hope the majority of our reads have a high average quality score with no large bumps at the lower quality values.
The next plot gives the “Per base sequence content”, which always gives a FAIL for RNA-seq data. This is because the first 10-12 bases result from the ‘random’ hexamer priming that occurs during RNA-seq library preparation. This priming is not as random as we might hope, giving an enrichment in particular bases for these initial nucleotides.
The “Per sequence GC content” plot gives the GC distribution over all sequences. Generally is a good idea to note whether the GC content of the central peak corresponds to the expected % GC for the organismLinks to an external site.. Also, the distribution should be normal unless over-represented sequences (sharp peaks on a normal distribution) or contamination with another organism (broad peak).


The mammalian genome is characterized by its high spatial heterogeneity in base composition. The average GC content of a 100-kb fragment of the human genome can be as low as 35% or as high as 60%, a range that is twice as wide as that typically observed in teleostean fishes, for instance 
The next module explores the number of duplicated sequences in the library. This plot can help identify a low complexity library, which could result from too many cycles of PCR amplification or too little starting material. For RNA-seq, we don’t normally do anything to address this in the analysis, but if this were a pilot experiment, we might adjust the number of PCR cycles, amount of input, or amount of sequencing for future libraries.
The “Overrepresented sequences” table is another important module, as it displays the sequences (at least 20 bp) that occur in more than 0.1% of the total number of sequences. This table aids in identifying contamination, such as vector or adapter sequences. If the %GC content was off in the above module, this table can help identify the source. If not listed as a known adapter or vector, it can help to BLAST the sequence to determine the identity.
I encourage you to look at the report for the full set of reads and note how the QC results differ when using the entire dataset.