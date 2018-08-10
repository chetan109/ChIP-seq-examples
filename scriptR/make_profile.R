#functions
source("src/Functions.R")
#packages
require("dplyr")
require("ggplot2")
require("reshape2")
require("rtracklayer")
library("BSgenome.Hsapiens.UCSC.hg19")
seqlens = seqlengths( Hsapiens );


##Code for one bed and two bigwigs
SITES <- import.bed("BLESS_80best_JunFragPE_Rmdups_pm500bp.bed")
#Get value for profile
bw.files <- c("trimmed_BLESS_U2OSpTam.sorted_fragments_800bp_rmdups_normalized.bw","trimmed_BLESS_U2OS-Tam.sorted_fragments_800bp_rmdups_normalized.bw")
names(bw.files) <- c("+4OHT","-4OHT")
window <- 5000
span <- 1

dat.plot <- NULL
for(n in names(bw.files)){
    
    wig <- import.bw(bw.files[[n]],as = "RleList")
    signal <- computeProfile(bed=SITES,wig=wig,w=window,span=span,seqlens = seqlens,method = "mean") %>% colMeans()
    sub.dat.plot <- data.frame(Window=seq( -window, window - span + 1, span ),
                               Value=signal,
                               Type=n)
    if(is.null(dat.plot)){
        dat.plot <- sub.dat.plot
    }else{
        dat.plot <- rbind(dat.plot,sub.dat.plot)
    }
}

##Code for two bed and one bigwig
SITES <- list("HR"=import.bed("BLESS_HR_JunFragPE_Rmdups_pm500bp.bed"),
              "NHEJ"=import.bed("BLESS_NHEJ_JunFragPE_Rmdups_pm500bp.bed")
)
#Get value for profile
bw.file <- "trimmed_BLESS_U2OSpTam.sorted_fragments_800bp_rmdups_normalized.bw"
names(bw.file) <- "+4OHT"


wig <- import.bw(bw.file,as = "RleList")

window <- 5000
span <- 1

dat.plot <- NULL
for(n in names(SITES)){
    
    
    signal <- computeProfile(bed=SITES[[n]],wig=wig,w=window,span=span,seqlens = seqlens,method = "mean") %>% colMeans()
    sub.dat.plot <- data.frame(Window=seq( -window, window - span + 1, span ),
                               Value=signal,
                               Type=n)
    if(is.null(dat.plot)){
        dat.plot <- sub.dat.plot
    }else{
        dat.plot <- rbind(dat.plot,sub.dat.plot)
    }
}

##Plot for each case

p <- ggplot(dat.plot,aes(x=Window,y=Value,color=Type)) +
    geom_line() +
    #geom_ribbon( aes( ymin = 0, ymax = Value ,fill=Type), alpha = 0.5 ) + #Add if you want shape under the plot line
    theme_classic(base_size = 18)
print(p)
