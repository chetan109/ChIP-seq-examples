### packages
library("BSgenome.Hsapiens.UCSC.hg19")
library("rtracklayer")
library("dplyr")
#### On assigne le genome à une variable
hg19 = BSgenome.Hsapiens.UCSC.hg19

#### Le motif recherché pour ASI
ASI = DNAString("GCGATCGC")

#### On enregistre le nom des differentes séquences (chromosomes) contenu dans le génome
chromosomes = seqnames(hg19)



file = "ASIsites_hg19.bed"

asi = data.frame()

for( chromosome in paste0("chr",c(1:22,"X")) ) {
    chr_hg19 = hg19[[chromosome]]
    cc = matchPattern(ASI, chr_hg19, max.mismatch=0)
    taille = length(cc)
    chr_names = rep(chromosome,taille)
    data = data.frame(seqnames=chr_names,start=start(cc),end=end(cc),check.names=FALSE)
    asi = rbind(asi,data)
}


asi <- asi%>% GRanges()
export.bed(asi,file)
