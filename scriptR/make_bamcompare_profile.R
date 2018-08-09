#Use BamCompare before to perform a logRatio between 2 bigwig (in or case pOHT/mOHT)
# c="FK2_GATC"
# dir=""
# pref="GATC_Merged_0316_"
# mOHT="FK2_minus"
# pOHT="FK2_OHT"
# suff=".rmdups.bam"
# 
# b1=$dir$pref$pOHT$suff
# b2=$dir$pref$mOHT$suff
# echo $c
# bamCompare -b1 $b1 -b2 $b2 -o $c"_bin50bp_log2ratio.bw" -of=bigwig --scaleFactorsMethod=readCount --ratio=log2 -bs=50 -p=32


#functions
source("src/Functions.R")
#packages
require("ggplot2")
require("reshape2")
require("rtracklayer")
library("BSgenome.Hsapiens.UCSC.hg19")
seqlens = seqlengths( Hsapiens );

SITES <- import.bed("BLESS_80best_JunFragPE_Rmdups_pm500bp.bed")
#Get value for heatmap
bw.file <- "FK2_GATC_bin50bp_log2ratio.bw"
window <- 5000000
span <- 50000
wig <- import.bw(bw.file,as = "RleList")
signal <- computeProfileLog( bed = SITES, wig = wig, w = window, span = span, seqlens = seqlens );
dat.plot <- data.frame(Windows=seq( -window, window - span + 1, span ),
                       value=signal)

whichcolor = c("#B8312F","#2969B0")

dat.plot$pos <- 0
dat.plot[dat.plot$value > 0 ,]$pos <- dat.plot[dat.plot$value > 0 ,]$value
dat.plot$neg <- 0
dat.plot[dat.plot$value < 0 ,]$neg <- dat.plot[dat.plot$value < 0 ,]$value
gp <- ggplot( na.omit( dat.plot ), aes( Windows, pos ) ) +
    labs( list( title = "", x = "", y = "" ) ) +
    geom_line( colour = whichcolor[1] ) +
    geom_ribbon( aes( ymin = 0, ymax = pos ),fill = whichcolor[1], alpha = 0.5 ) +
    geom_line( aes( Windows, neg ),colour = whichcolor[2] ) +
    geom_ribbon( aes( ymin = 0, ymax = neg ), fill = whichcolor[2], alpha = 0.5 ) +
    geom_hline( aes( yintercept = 0 ) ) +
    theme_classic(base_size = 18)
print( gp )
