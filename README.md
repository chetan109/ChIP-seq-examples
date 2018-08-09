# Comprehensive mapping of histone modifications at DNA Double Strand Breaks deciphers repair pathway chromatin signatures.

## Overview

Scripts used to generate all papers figures and some raw data, like AsiSI search motif and DSB classes (HR/NHEJ).

## System requirements

### Alignment

The workflow was written with snakemake (python), and need some genomic tools : 

* `bwa 0.7.12-r1039`.
* `samtools 1.4`.
* `FastQC v0.11.5`.
* `R version 3.4.0` :
	- `library(Rsamtools)`.
        - `library(GenomicAlignments)`.
        - `library(rtracklayer)`.
        - `library(BSgenome.Hsapiens.UCSC.hg19)`.



