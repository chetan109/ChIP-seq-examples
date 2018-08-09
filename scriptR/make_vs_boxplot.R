#functions
source("src/Functions.R")
#packages
require("ggplot2")
require("reshape2")
require(rtracklayer)
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
    c("HMNY7BGXX_H2BK120ac_DIVA_16s000420-1-1_Clouaire_lane116s000420_sequence_normalized.bw","HMNY7BGXX_H2BK120ac_OHT_DI_16s000421-1-1_Clouaire_lane116s000421_sequence_normalized.bw"),
    c("NG-8008_FK2_minus_lib76330_3801_7_1_normalized.bw","NG-8008_FK2_OHT_lib76331_3801_6_1_normalized.bw")
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

windows <- list("2kb"=2000)

#List creation : LogFC, pOHT and mOHT profile signal for each Histone with specific windows
rlmOHT = list();
for( n in names(FILES) ){
    rlmOHT[[n]] = list();
}


#For each Histone
for(n in prots){
    #Read the bigwig
    myWig = readData( dir = "", fnpOHT = FILES[[n]][2], fnmOHT = FILES[[n]][1] );
    for(bedname in names(SITES)){
        tmpMOHT = compute1ValPerSite( bed = SITES[[bedname]], wig = myWig[["mOHT"]], w = window, seqlens );
        rlmOHT[[n]][[bedname]] = tmpMOHT;
    }
}
##PLOT
dat.plot <- data.frame("Histones"=NULL,"HR"=NULL,"NHEJ"=NULL)
for(n in prots){
    dat.HR <- rlmOHT[[n]][["HR"]]
    dat.NHEJ <- rlmOHT[[n]][["NHEJ"]]
    dat.plot <- rbind(dat.plot,data.frame("Histones"=n,"HR"=dat.HR,"NHEJ"=dat.NHEJ))
    
}

sub.dat.plot <- melt(dat.plot)
sub.dat.plot$Log <- log2(sub.dat.plot$value)

whichcolor = c("#F5AB35","#049372")
B = boxplot(Log~variable+Histones,data=sub.dat.plot,xaxt = "n",outline = FALSE ,col = rep(whichcolor,19))
points(B$group, B$out, type = "p", pch=20,cex = 0.55)
axis(1,las = 3,at = seq(1.5,(19*2),by=2),labels = unique(sub.dat.plot$Histones))
