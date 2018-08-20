# Alignment workflow on our NGS data

## ChIP-SEQ workflow

### Files

* `Snakefile` : main script for snakemake, call all RULES and config file.
* `config.json` : config file which contain the list and location of all ChIP-SEQ and references files.
* `RULES/fastqc.rules` : All cmd from FastQC program (quality control)
* `RULES/bwa.rules` : All cmd from bwa (alignment to reference genome)
* `RULES/samtools.rules` : All cmd from samtools (sam/bam manipulation)
* `RULES/bigwig.rules` : R script call by snakemake to compute coverage in BigWig format.

### Minimal configuration

* Open config file (`config.json`) and type a directory where your data will be written (default `PROCESSED`).
* Indicate the directory were the raw files (`fastq` format) are stored (default `RAW`).
* Indicate the location, the name and the extension of the reference genome (default `""`,`female.hg19`,`.fa`).

### Run the workflow

Cmd : `snakemake -s Snakefile`.

### Rulegraph

![rulegraph](https://github.com/LegubeDNAREPAIR/HistoneMapping/blob/master/script/alignment/rulegraph_chipseq.png)

## BLESS workflow

### Files

* `Snakefile_BLESS` : main script for snakemake, call all RULES and config file.
* `config_BLESS.json` : config file which contain the list and location of all ChIP-SEQ and references files.
* `RULES/bwa_paired.rules` : All cmd from bwa for paired-end data (alignment to reference genome)
* `RULES/samtools.rules` : All cmd from samtools (sam/bam manipulation)
* `RULES/bless.rules` : R script call by snakemake to merge paired-end reads into fragments, remove small and non-unique fragments. Then compute coverage in BigWig format.

### Minimal configuration

* Open config file (`config.json`) and type a directory where your data will be written (default `BLESS_JUNE`).
* Indicate the directory were the raw files (`fastq` format) are stored (default `RAW`).
* Indicate the location, the name and the extension of the reference genome (default `""`,`female.hg19`,`.fa`).

### Run the workflow

Cmd : `snakemake -s Snakefile_BLESS`.

### Rulegraph

![rulegraph](https://github.com/LegubeDNAREPAIR/HistoneMapping/blob/master/script/alignment/rulegraph_BLESS.png)
