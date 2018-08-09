#functions
source("src/Functions.R")
#packages
require("ggplot2")
require("reshape2")
require("rtracklayer")
library("Homo.sapiens")

SITES <- import.bed("BLESS_80best_JunFragPE_Rmdups_pm500bp.bed")


#Specific order

FILES <- c(
    "C8LVDACXX_H2AZac_DIVA_16s000910-1-1_Clouaire_lane116s000910_sequence_normalized.bw",
    "H4K16ac.clean_normalized.bw",
    "H3K56ac.clean_normalized.bw",
    "H3K36me3.clean_normalized.bw",
    "H3K4me2_-OHT.clean_normalized.bw",
    "H3K9me2_normalized.bw",
    "HWNTLBGX3_H3K9me3_DIVA_17s005723-1-1_Clouaire_lane117s005723_sequence_normalized.bw",
    "HN25TBGXX_H1_minus_DIVA_16s000414-1-1_Clouaire_lane116s000414_sequence_normalized.bw",
    "HFJNHBGX5_H4K20me2_DIVA_17s005866-1-1_Clouaire_lane117s005866_sequence_normalized.bw",
    "H3K79me2.clean_normalized.bw",
    "C8LVDACXX_H3_DIVA_16s000906-1-1_Clouaire_lane116s000906_sequence_normalized.bw",
    "H3K36me2.clean_normalized.bw",
    "C8LVDACXX_H4K12ac_DIVA_16s000908-1-1_Clouaire_lane116s000908_sequence_normalized.bw",
    "HMNY7BGXX_H2AZ_DIVA_16s000418-1-1_Clouaire_lane116s000418_sequence_normalized.bw",
    "H2Bub.clean_normalized.bw",
    "H4K20me1Mono_normalized.bw",
    "H4S1P_normalized.bw",
    "HMNY7BGXX_MacroH2A_DIVA_16s000416-1-1_Clouaire_lane116s000416_sequence_normalized.bw",
    "HMNY7BGXX_H2BK120ac_DIVA_16s000420-1-1_Clouaire_lane116s000420_sequence_normalized.bw"
    )

prots <- c( "H2AZac",
            "H4K16ac",
            "h3k56ac",
            "h3k36me3",
            "H3K4me2",
            "H3K9me2",
            "H3K9me3",
            "H1",
            "H4K20me2",
            "h3k79me2",
            "H3",
            "h3k36me2",
            "H4K12ac",
            "H2AZ",
            "h2bub",
            "H4K20me1Mono",
            "H4S1P",
            "MacroH2A",
            "H2BK120ac"
)
names(FILES) <- prots





hg19.genes <- genes(Homo.sapiens)
seqlens <- seqlengths(Homo.sapiens)


profs <- list()
for(n in prots){
    #Read the bigwig
    myWig = import.bw(FILES[[n]],as="RleList")
    #Compute profil
    profs[[n]] <- computeProfileGenes(bed = hg19.genes, wig = myWig, seqlens = seqlens , fun="mean",get="all")
}


polII.dir <- "polII_normalized_hg19.bw"

myWig = import.bw(polII.dir,as="RleList")
profs[["POLII"]] <- computeProfileGenes(bed = hg19.genes, wig = myWig, seqlens = seqlens , fun="mean",get="all")


POLII.cat.name <- c("Low","Mid","High")
POLII.cat <- bin.var(x = rowMeans(profs[["POLII"]]),bins = 3,method = "proportions",labels = POLII.cat.name)
table(POLII.cat)

subs.profs <- list()

for(cat in POLII.cat.name){
    subs.profs[[cat]] <- list()
    for(n in names(profs)){
        if(n=="POLII"){
            next
        }
        my.val <- profs[[n]][which(POLII.cat == cat),]
        subs.profs[[cat]][[n]] <- my.val
    }
}


dat.plot <- data.frame("value"=NULL,"Windows"=NULL,"Histone"=NULL,"Condition"=NULL,"Cat"=NULL)

for(n in names(profs)){
    if(n=="POLII"){
        next
    }
    for(cat in POLII.cat.name){
        dat.plot <- rbind(dat.plot,
                          data.frame("value"=colMeans(subs.profs[[cat]][[n]]),"Windows"=c(1:130),"Histone"=Name[[n]],"Condition"=Cond[[n]],"Cat"=cat))
    }
}

mycolor <- c("#27ae60","#2980b9","#c0392b")
p<- ggplot(dat.plot,aes(Windows,value,colour=Cat)) +
    geom_line() +
    geom_vline(xintercept = 15,linetype="longdash") +
    geom_vline(xintercept = 115,linetype="longdash") +
    scale_colour_manual(values=mycolor) +
    scale_fill_manual(values=mycolor) +
    scale_x_continuous(name = 'Position',
                       breaks = c(1,15,115,130),
                       labels = c("TSS-3000bp", 'TSS', 'TES', 'TES+3000bp')
    ) +
    facet_wrap(~Histone+Condition,ncol=3,scales="free_y") +
    theme_classic()

print(p)
