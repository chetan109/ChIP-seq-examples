### packages
library("BSgenome.Hsapiens.UCSC.hg19")
library("rtracklayer")
library("dplyr")
source("src/Functions.R")
#Compute classes HR/NHEJ based on RAD51 and XRCC4 dataset (Aymard et al, 2014. E-MTAB-1241)
bless80 <- import.bed("BLESS_80best_JunFragPE_Rmdups_pm500bp.bed")
bw.files <- c("RAD51_hg19_normalized.bw",
              "XRCC4_hg19_normalized.bw",
              "HVL7VBGX5_Lig4_2_OHT_18s000869-1-1_Clouaire_lane118s000869_sequence_normalized.bw")
names(bw.files) <- c("RAD51","XRCC4","Lig4")
windows <- c(2000,250,500)
names(windows) <- names(bw.files)

cat.dat <- NULL
for(n in names(bw.files)){
    wig <- import.bw(bw.files[n],as = "RleList")
    my.w <- windows[n]
    count <- compute1ValPerSiteForBLESS(bed=bless80,wig=wig,w=my.w,seqlens = seqlens,fun="mean")
    if(is.null(cat.dat)){
        cat.dat <- data.frame("Value"=count)
        colnames(cat.dat) <- n
    }else{
        oldcol <- colnames(cat.dat)
        cat.dat <- cbind(cat.dat,count)
        colnames(cat.dat) <- c(oldcol,n)
    }
}

#Compute ratio for each
cat.dat <- cat.dat %>%
    mutate(Ratio_RAD51.XRCC4 = (RAD51+0.01)/(XRCC4+0.01)) %>%
    mutate(Ratio_RAD51.Lig4 = (RAD51+0.01)/(Lig4+0.01))

#Extract classes for RAD51/XRCC4 ratio
ratio.pos <- cat.dat %>% mutate(pos = row_number())%>% arrange(desc(Ratio_RAD51.XRCC4))
#Extract HR sites
HR.sites <- ratio.pos  %>% dplyr::slice(1:30) %>% pull(pos)
HR.sites <- bless80[HR.sites] %>% sortSeqlevels() %>% sort()

export.bed(HR.sites,"BLESS_HR_JunFragPE_Rmdups_pm500bp.bed")
#Extract NHEJ sites
NHEJ.sites <- ratio.pos  %>% dplyr::slice(51:80) %>% pull(pos)
NHEJ.sites <- bless80[NHEJ.sites] %>% sortSeqlevels() %>% sort()

export.bed(NHEJ.sites,"BLESS_NHEJ_JunFragPE_Rmdups_pm500bp.bed")

#Extract classes for RAD51/Lig4 ratio
ratio.pos <- cat.dat %>% mutate(pos = row_number())%>% arrange(desc(Ratio_RAD51.Lig4))
#Extract HR sites
HR.sites <- ratio.pos  %>% dplyr::slice(1:30) %>% pull(pos)
HR.sites <- bless80[HR.sites] %>% sortSeqlevels() %>% sort()

export.bed(HR.sites,"BLESS_HR_JunFragPE_Rmdups_pm500bp_Lig4.bed")
#Extract NHEJ sites
NHEJ.sites <- ratio.pos  %>% dplyr::slice(51:80) %>% pull(pos)
NHEJ.sites <- bless80[NHEJ.sites] %>% sortSeqlevels() %>% sort()

export.bed(NHEJ.sites,"BLESS_NHEJ_JunFragPE_Rmdups_pm500bp_Lig4.bed")
