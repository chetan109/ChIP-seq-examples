#functions
source("src/Functions.R")
#packages
require("ggplot2")
require("reshape2")
require(rtracklayer)
library("BSgenome.Hsapiens.UCSC.hg19")
seqlens = seqlengths( Hsapiens );

SITES <- import.bed("ASIsites_hg19.bed")
##Sort by specific chip-seq value
bw <- "trimmed_BLESS_U2OSpTam.sorted_fragments_800bp_rmdups_normalized.bw"
wig <- import.bw(bw,as = "RleList")
my.w <- 500
count <- compute1ValPerSite(bed=SITES,wig=wig,w=my.w,seqlens = seqlens,fun="mean")
asi_order = as.data.frame(SITES) %>% mutate(value = count) %>% arrange(count) %>% pull(name)


#Get value for heatmap
bw.file <- "trimmed_BLESS_U2OSpTam.sorted_fragments_800bp_rmdups_normalized.bw"
window <- 5000
span <- 1
wig <- import.bw(bw.file,as = "RleList")
#HM and profile
signal <- computeProfile(bed=SITES,wig=wig,w=window,span=span,seqlens = seqlens,method = "mean")
colnames(signal) <- seq( -window, window - span + 1, span )
rownames(signal) <- SITES$name

signal <- melt(signal)
colnames(signal) <- c("Sites","Window","Value")
signal$Sites <- factor(signal$Sites,levels = asi_order) #Order by specific value


background = "white"
highvalues = "#551A8B"


p = ggplot(signal,aes(Window,Sites)) + 
    geom_tile(aes(fill = Value)) + 
    scale_fill_gradient(low = background,high = highvalues,na.value=background) + 
    labs(list(title = "", x = "", y = "")) +
    theme_classic()+
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())

print(p)