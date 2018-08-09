#functions
source("src/Functions.R")
#packages
require("ggplot2")
require("reshape2")
require(rtracklayer)
library("BSgenome.Hsapiens.UCSC.hg19")
seqlens = seqlengths( Hsapiens );
SITES <- import.bed("ASIsites_hg19.bed")
#Get value for jitter_plot
bw.file <- "HWNTLBGX3_H3K9me3_DIVA_17s005723-1-1_Clouaire_lane117s005723_sequence_normalized.bw"
names(bw.file) <- "H3K9me3"
wig <- import.bw(bw.file,as = "RleList")
window <- 5000

signal <- compute1ValPerSite(bed=SITES,wig=wig,w=window,seqlens = seqlens,fun = "mean")

dat.plot <- data.frame(SITES = SITES$name,Value=signal,File=names(bw.file))

p <- ggplot(dat.plot,aes(x=File,y=Value)) + geom_jitter() + theme_classic()
print(p)
#Add some point from another bed
SITES <- import.bed("BLESS_80best_JunFragPE_Rmdups_pm500bp.bed")

signal <- compute1ValPerSite(bed=SITES,wig=wig,w=window,seqlens = seqlens,fun = "mean")
new.dat.plot <- data.frame(SITES = SITES$name,Value=signal,File=names(bw.file))


p <- p + geom_jitter(data=new.dat.plot,aes(x=File,y=Value),col="red")
print(p)
