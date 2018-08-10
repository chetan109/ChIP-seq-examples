### On charge la library pour h19
library("BSgenome.Hsapiens.UCSC.hg19")

#### On assigne le genome à une variable
hg19 = BSgenome.Hsapiens.UCSC.hg19

#### Le motif recherché pour ASI
ASI = DNAString("GCGATCGC")

#### On enregistre le nom des differentes séquences (chromosomes) contenu dans le génome
chromosomes = seqnames(hg19)



file = "ASIsites_hg19.bed"

full_data = data.frame()

for( chromosome in paste0("chr",c(1:22,"X")) ) {
    chr_hg19 = hg19[[chromosome]]
    cc = matchPattern(ASI, chr_hg19, max.mismatch=0)
    taille = length(cc)
    chr_names = rep(chromosome,taille)
    data = data.frame(chromosome_name=chr_names,start=start(cc),end=end(cc),check.names=FALSE)
    full_data = rbind(full_data,data)
}


#Dans l'ancien jeu de données, le site 33 est une duplication du site 35, pour garder les mêmes noms, on passe de SITE32 à SITE34 directement
names = c(paste("SITE",1:32,sep=""),paste("SITE",34:1212,sep=""))

full_data = cbind(full_data,names)

write.table(full_data, file=file, quote=FALSE, sep="\t",row.names=FALSE, col.names=FALSE)
