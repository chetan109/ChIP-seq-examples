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
#Get value for boxplot
bw.files <- c("trimmed_BLESS_U2OSpTam.sorted_fragments_800bp_rmdups_normalized.bw","trimmed_BLESS_U2OS-Tam.sorted_fragments_800bp_rmdups_normalized.bw")
names(bw.files) <- c("+4OHT","-4OHT")
window <- 500

dat.plot <- NULL
for(n in names(bw.files)){
    
    wig <- import.bw(bw.files[[n]],as = "RleList")
    signal <- compute1ValPerSite(bed=SITES,wig=wig,w=window,seqlens = seqlens,fun = "sum")
    sub.dat.plot <- data.frame(Value=signal,
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
#Get value for boxplot
bw.file <- "trimmed_BLESS_U2OSpTam.sorted_fragments_800bp_rmdups_normalized.bw"


wig <- import.bw(bw.file,as = "RleList")

window <- 500

dat.plot <- NULL
for(n in names(SITES)){
    
    
    signal <- signal <- compute1ValPerSite(bed=SITES[[n]],wig=wig,w=window,seqlens = seqlens,fun = "sum")
    sub.dat.plot <- data.frame(Value=signal,
                               Type=n)
    if(is.null(dat.plot)){
        dat.plot <- sub.dat.plot
    }else{
        dat.plot <- rbind(dat.plot,sub.dat.plot)
    }
}

##Plot for each case
B = boxplot(data=dat.plot,Value~Type,outline=FALSE)
points(B$group, B$out, type = "p", pch=20,cex = 1)

##Specific boxplot for 2 bed and 2 dataset (Fig1C and Fig7E)
SITES <- list("uncut"=import.bed("ASIsites_hg19.bed"),
              "cut"=import.bed("BLESS_80best_JunFragPE_Rmdups_pm500bp.bed")
)
#Uncut sites are AsiSI not in best 80 bless list
SITES[["uncut"]] <- SITES[["uncut"]][!SITES[["uncut"]]$name %in% SITES[["cut"]]$name]

bw.files <- c("trimmed_BLESS_U2OSpTam.sorted_fragments_800bp_rmdups_normalized.bw","trimmed_BLESS_U2OS-Tam.sorted_fragments_800bp_rmdups_normalized.bw")
names(bw.files) <- c("+4OHT","-4OHT")
window <- 500

dat.plot <- NULL
for(t in names(SITES)){
    for(n in names(bw.files)){
        
        wig <- import.bw(bw.files[[n]],as = "RleList")
        signal <- compute1ValPerSite(bed=SITES[[t]],wig=wig,w=window,seqlens = seqlens,fun = "sum")
        sub.dat.plot <- data.frame(Value=signal,
                                   Type=t,
                                   File=n)
        if(is.null(dat.plot)){
            dat.plot <- sub.dat.plot
        }else{
            dat.plot <- rbind(dat.plot,sub.dat.plot)
        }
    }
}

dat.plot <- dat.plot %>% mutate(log10val = log10(Value))
dat.plot$Type <- factor(dat.plot$Type,levels = c("uncut","cut"))
dat.plot$File <- factor(dat.plot$File,levels = c("-4OHT","+4OHT"))

mycols <- c("#BEBEBE","#BEBEBE","#FDBECD","#FDBECD")
boxplot(data=dat.plot,log10val~File+Type,outline=FALSE,names=rep(c("-4OHT","+4OHT"),2),col = mycols)
