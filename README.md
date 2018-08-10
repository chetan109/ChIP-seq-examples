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
  * `library(Rsamtools)`.
  * `library(GenomicAlignments)`.
  * `library(rtracklayer)`.
  * `library(BSgenome.Hsapiens.UCSC.hg19)`.

### Figure generation

the scripts were written with `R`, and need some packages : 

* `library(rtracklayer)`.
* `library(BSgenome.Hsapiens.UCSC.hg19)`.
* `library(ggplot2)`.
* `library(reshape2)`.
* `library(dplyr)`.
* `library(BSgenome.Hsapiens.UCSC.hg19)`.
* `library(Homos.sapiens)`.

| Script                    | Description                                                                          | Figures                           |
|---------------------------|--------------------------------------------------------------------------------------|-----------------------------------|
| get_asi_location.R        | Search AsiSI motif across genome.                                                    | Table S1                          |
| get_repair_classes.R      | Compute HR and NHEJ classes based on RAD51/XRCC4 and RAD51/Lig4 ratios.              |                                   |
| index_BLESS.R             | Plot sorted BLESS value for each AsiSI location.                                     | Figure S1B                        |
| make_bamcompare_profile.R | Plot LogRatio coverage for given bigwig.                                             | Figures 2B,3B,6B,S2D              |
| make_boxplot.R            | Plot the distribution for one/two bigwig with one/two bed.                           | Figures 1C,4B,6C,7B,7E,S4E        |
| make_circle_plot.R        | Compute and plot enrichment for HR sites for several windows in multiple bigwigs.    | Figures 4F,S5E                    |
| make_heatmap_sorted.R     | Plot sorted heatmap for given bigwig.                                                | Figures 1B,4C,7D,S4D,S6C          |
| make_jitter_plot.R        | Plot a jittered distribution for a given bigwig with one bed file.                   | Figure S1D                        |
| make_log_boxplot.R        | Plot the log ratio distribution of multiple bigwigs for one bed file.                | Figures 2A, 3A, 5A, 5B            |
| make_profile.R            | Plot average profile for a given bigwig with one bed file.                           | Figures 1D,2B,4A,4D,6B,7C,S6A,S6F |
| make_profil_gene_cat.R    | Plot average profile for multiple bigwig over human genes sorted by expression level | Figure S2B                        |
| make_vs_boxplot.R         | Plot the distribution of multiple bigwigs for two bed files                          | Figures 4E,S5D                    |
