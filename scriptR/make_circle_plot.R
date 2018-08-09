#functions
source("src/Functions.R")
#packages
require("ggplot2")
require("reshape2")
require("rtracklayer")
library("BSgenome.Hsapiens.UCSC.hg19")
seqlens = seqlengths( Hsapiens );

SITES <- list("HR"=import.bed("BLESS_HR_JunFragPE_Rmdups_pm500bp.bed"),
              "NHEJ"=import.bed("BLESS_NHEJ_JunFragPE_Rmdups_pm500bp.bed")
)
#Specific order

FILES <- list(
    c("C8LVDACXX_H2AZac_DIVA_16s000910-1-1_Clouaire_lane116s000910_sequence_normalized.bw","C8LVDACXX_H2AZac_OHT_DIVA_16s000911-1-1_Clouaire_lane116s000911_sequence_normalized.bw"),
    c("H4K16ac.clean_normalized.bw","H4K16ac_OHT.clean_normalized.bw"),
    c("H3K56ac.clean_normalized.bw","H3K56ac_OHT.clean_normalized.bw"),
    c("H3K36me3.clean_normalized.bw","H3K36me3_OHT.clean_normalized.bw"),
    c("H3K4me2_-OHT.clean_normalized.bw","H3K4me2_OHT.clean_normalized.bw"),
    c("H3K9me2_normalized.bw","H3K9me2_OHT_normalized.bw"),
    c("HWNTLBGX3_H3K9me3_DIVA_17s005723-1-1_Clouaire_lane117s005723_sequence_normalized.bw","HWNTLBGX3_H3K9me3_OHT_DIVA_17s005724-1-1_Clouaire_lane117s005724_sequence_normalized.bw"),
    c("HN25TBGXX_H1_minus_DIVA_16s000414-1-1_Clouaire_lane116s000414_sequence_normalized.bw","HN25TBGXX_H1_OHT_DIVA_16s000415-1-1_Clouaire_lane116s000415_sequence_normalized.bw"),
    c("HFJNHBGX5_H4K20me2_DIVA_17s005866-1-1_Clouaire_lane117s005866_sequence_normalized.bw","HFJNHBGX5_H4K20me2_OHT_17s005867-1-1_Clouaire_lane117s005867_sequence_normalized.bw"),
    c("H3K79me2.clean_normalized.bw","H3K79me2_OHT.clean_normalized.bw"),
    c("C8LVDACXX_H3_DIVA_16s000906-1-1_Clouaire_lane116s000906_sequence_normalized.bw","C8LVDACXX_H3_OHT_DIVA_16s000907-1-1_Clouaire_lane116s000907_sequence_normalized.bw"),
    c("H3K36me2.clean_normalized.bw","H3K36me2_OHT.clean_normalized.bw"),
    c("C8LVDACXX_H4K12ac_DIVA_16s000908-1-1_Clouaire_lane116s000908_sequence_normalized.bw","C8LVDACXX_H4K12ac_OHT_DIVA_16s000909-1-1_Clouaire_lane116s000909_sequence_normalized.bw"),
    c("HMNY7BGXX_H2AZ_DIVA_16s000418-1-1_Clouaire_lane116s000418_sequence_normalized.bw","HMNY7BGXX_H2AZ_OHT_DIVA_16s000419-1-1_Clouaire_lane116s000419_sequence_normalized.bw"),
    c("H2Bub.clean_normalized.bw","H2Bub_OHT.clean_normalized.bw"),
    c("H4K20me1Mono_normalized.bw","H4K20me1Mono_OHT_normalized.bw"),
    c("H4S1P_normalized.bw","H4S1P_OHT_normalized.bw"),
    c("HMNY7BGXX_MacroH2A_DIVA_16s000416-1-1_Clouaire_lane116s000416_sequence_normalized.bw","HMNY7BGXX_MacroH2A_OHT_DIV_16s000417-1-1_Clouaire_lane116s000417_sequence_normalized.bw"),
    c("HMNY7BGXX_H2BK120ac_DIVA_16s000420-1-1_Clouaire_lane116s000420_sequence_normalized.bw","HMNY7BGXX_H2BK120ac_OHT_DI_16s000421-1-1_Clouaire_lane116s000421_sequence_normalized.bw")
)


prots <- c("h3k36me3",
           "h2bub",
           "h3k79me2",
           "H4K20me1Mono",
           "H3K4me2",
           "H3",
           "H4K12ac",
           "H4K20me2",
           "H2BK120ac",
           "H2AZac",
           "h3k56ac",
           "MacroH2A",
           "H4K16ac",
           "h3k36me2",
           "H2AZ",
           "H4S1P",
           "H1",
           "H3K9me2",
           "H3K9me3"
)


names(FILES) <- prots


windows <- list("5000bp"=5000,"2000bp"=2000,"1000bp"=1000,"500bp"=500)
inter=0.05
maxi=0.01
Sizes <- c(1,5,10)
#List creation : LogFC, pOHT and mOHT profile signal for each Histone with specific windows
rlpOHT = rlmOHT = rllogFC = list();
for( n in names(FILES) ){
    rlpOHT[[n]] = rlmOHT[[n]] = rllogFC[[n]] = list();
}




#For each Histone
for(n in prots){
    #Read the bigwig
    myWig = readData( dir = "", fnpOHT = FILES[[n]][2], fnmOHT = FILES[[n]][1] );
    for(bedname in names(SITES)){
        for(w in names(windows)){
            window = windows[[w]]
            #compute and stock logFC
            tmpPOHT = compute1ValPerSite( bed = SITES[[bedname]], wig = myWig[["pOHT"]], w = window, seqlens );
            tmpMOHT = compute1ValPerSite( bed = SITES[[bedname]], wig = myWig[["mOHT"]], w = window, seqlens );
            tmpLog = log2( tmpPOHT / tmpMOHT );
            tmpLog[tmpPOHT == 0 & tmpMOHT == 0] = 0;
            tmpLog[tmpPOHT != 0 & tmpMOHT == 0] = max( tmpLog[is.finite( tmpLog )] ) + 1;
            tmpLog[tmpPOHT == 0 & tmpMOHT != 0] = min( tmpLog[is.finite( tmpLog )] ) - 1;
            rlpOHT[[n]][[bedname]][[w]] = tmpPOHT;
            rlmOHT[[n]][[bedname]][[w]] = tmpMOHT;
            rllogFC[[n]][[bedname]][[w]] = tmpLog;
        }
    }
}


#Compute data for HR/NHEJ comparison
dat.plot <- data.frame("Histones"=NULL,"Windows"=NULL,"Type"=NULL,"P.val"=NULL)
for(n in prots){
    for(w in names(windows)){
        dat.HR <- rlmOHT[[n]][["HR"]][[w]]
        dat.NHEJ <- rlmOHT[[n]][["NHEJ"]][[w]]
        #Point radius computation
        pval.1 <- wilcox.test(dat.HR,dat.NHEJ)$p.value
        Size <- NULL
        if(pval.1>inter){
            Size <- Sizes[1]
            Type <- "None"
            dat.plot <- rbind(dat.plot,data.frame("Histones"=n,"Windows"=w,"Type"=Type,"P.val"=Size))
            next;
        }else {
            if(pval.1 > maxi){
                Size <- Sizes[2]
            }else {
                Size <- Sizes[3]
            }
        }
        #Point color computation
        pval.2 <- wilcox.test(dat.HR,dat.NHEJ,alternative = "greater")$p.value
        if(pval.2 > inter){
            pval.3 <- wilcox.test(dat.HR,dat.NHEJ,alternative = "less")$p.value
            if(pval.3 > inter){
                Type <- "None"
            }else {
                Type <- "NHEJ"
            }
            
        }else{
            Type <- "HR"
        }
        dat.plot <- rbind(dat.plot,data.frame("Histones"=n,"Windows"=w,"Type"=Type,"P.val"=Size))
    }
}

dat.plot$Windows <- factor(dat.plot$Windows, levels = c("500bp","1000bp","2000bp","5000bp"))
dat.plot$Histones <- factor(dat.plot$Histones, levels = rev(levels(dat.plot$Histones)))
dat.plot$Type <- factor(dat.plot$Type, levels = c("HR","NHEJ","None"))


whichcolor = c("#F5AB35","#049372")

p <- ggplot(dat.plot,aes(Windows,Histones,colour=Type,size = P.val)) +
    geom_point() +
    scale_colour_manual(values=c(whichcolor,"black"),drop=FALSE,guide = guide_legend(title = "Enriched in")) +
    scale_size(range = c(1,10),breaks = Sizes,labels = c("NS","0.01<p.val<0.05","<=0.01"),guide = guide_legend(title = "P value (HR vs NHEJ)")) +
    theme_bw() +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          #legend.position="none",
          text = element_text(size=16),
          strip.background = element_rect(colour = "white", fill = "white"),
          panel.spacing = unit(0, "lines")
    ) 

print(p)

#Compute data for pOHT/mOHT comparison with specific bed
type <- "HR"
dat.plot <- data.frame("Histones"=NULL,"Windows"=NULL,"Type"=NULL,"P.val"=NULL)
for(n in prots){
    for(w in names(windows)){
        dat.pOHT <- rlpOHT[[n]][[type]][[w]]
        dat.mOHT <- rlmOHT[[n]][[type]][[w]]
        #Point radius computation
        pval.1 <- wilcox.test(dat.pOHT,dat.mOHT)$p.value
        Size <- NULL
        if(pval.1>inter){
            Size <- Sizes[1]
            Type <- "None"
            dat.plot <- rbind(dat.plot,data.frame("Histones"=n,"Windows"=w,"Type"=Type,"P.val"=Size))
            next;
        }else {
            if(pval.1 > maxi){
                Size <- Sizes[2]
            }else {
                Size <- Sizes[3]
            }
        }
        #Point color computation
        pval.2 <- wilcox.test(dat.pOHT,dat.mOHT,alternative = "greater")$p.value
        if(pval.2 > inter){
            pval.3 <- wilcox.test(dat.pOHT,dat.mOHT,alternative = "less")$p.value
            if(pval.3 > inter){
                Type <- "None"
            }else {
                Type <- "mOHT"
            }
            
        }else{
            Type <- "pOHT"
        }
        dat.plot <- rbind(dat.plot,data.frame("Histones"=n,"Windows"=w,"Type"=Type,"P.val"=Size))
    }
}
dat.plot$Windows <- factor(dat.plot$Windows, levels = c("500bp","1000bp","2000bp","5000bp"))
dat.plot$Histones <- factor(dat.plot$Histones, levels = rev(levels(dat.plot$Histones)))
dat.plot$Type <- factor(dat.plot$Type, levels = c("pOHT","mOHT","None"))

whichcolor = c("#B8312F","#2969B0")
p <- ggplot(dat.plot,aes(Windows,Histones,colour=Type,size = P.val)) +
    geom_point() +
    scale_colour_manual(values=c(whichcolor,"black"),drop=FALSE,guide = guide_legend(title = "Enriched in")) +
    scale_size(range = c(1,10),breaks = Sizes,labels = c("NS","0.01<p.val<0.05","<=0.01"),guide = guide_legend(title = "P value (pOHT vs mOHT)")) +
    theme_bw() +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          #legend.position="none",
          text = element_text(size=16),
          strip.background = element_rect(colour = "white", fill = "white"),
          panel.spacing = unit(0, "lines")
    ) 
print(p)
