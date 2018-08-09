#functions
source("src/Functions.R")
#packages
require("dplyr")
require("ggplot2")
require("reshape2")
require("rtracklayer")
library("BSgenome.Hsapiens.UCSC.hg19")
seqlens = seqlengths( Hsapiens );
SITES <- import.bed("ASIsites_hg19.bed")
#Get value for jitter_plot
bw.file <- "/mnt/volumes/NAS1/DATA/ChIP-Seq/BLESS/CHIP-SEQ/Jun2015/PROCESSED/ALIGNED_PAIRED/WIGGLE/trimmed_BLESS_U2OSpTam.sorted_fragments_800bp_rmdups_normalized.bw"
names(bw.file) <- "BLESS"
wig <- import.bw(bw.file,as = "RleList")
window <- 500

signal <- compute1ValPerSite(bed=SITES,wig=wig,w=window,seqlens = seqlens,fun = "mean")

dat.plot <- data.frame(SITES = SITES$name,Value=signal) %>% arrange(Value) %>% mutate(Index=1:length(SITES))

p <- ggplot(dat.plot,aes(x=Index,y=Value)) + geom_point() + theme_classic()
print(p)

plot(asiBLESS$BLESS,col=ifelse(asiBLESS$best_80=="YES","red","black"),pch=20,ylab="Average BLESS Seq",main="Bless Distribution for Asi sites")
abline(v=(length(asiBLESS$BLESS)-80),col="red")
abline(v=(length(asiBLESS$BLESS)-174),col="blue",lty="dashed")
legend("topleft",legend=c("outliers","best 80"),col=c("blue","red"),lty=c(2,1))


#Add some point from another bed
SITES2 <- import.bed("BLESS_80best_JunFragPE_Rmdups_pm500bp.bed")
new.dat.plot <- dat.plot %>% filter(SITES %in% SITES2$name)


p <- p + geom_point(data=new.dat.plot,aes(x=Index,y=Value),col="red") + 
    geom_vline(xintercept = max(dat.plot$Index)-80,col="red") +
    geom_vline(xintercept = max(dat.plot$Index)-174,col="blue",linetype="dashed")
print(p)
