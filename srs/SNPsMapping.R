setwd("Desktop/Entropy/")

################################# EXTRACT PATIENTS ############################################
status   = read.csv("Final_Cohort.csv",sep=";")
patients = read.table("ADNI3_PLINK_Final.fam", sep=" ")

adni3 = merge(patients,status,by.x="V2",by.y="Subject.ID")

write.table(adni3[adni3$Research.Group=="AD",c(2,1)],"AD_patients.txt",quote=F,col.names =F,row.names = F,sep="\t")
write.table(adni3[adni3$Research.Group=="CN",c(2,1)],"CN_patients.txt",quote=F,col.names =F,row.names = F,sep="\t")
write.table(adni3[adni3$Research.Group=="MCI",c(2,1)],"MCI_patients.txt",quote=F,col.names =F,row.names = F,sep="\t")

path = "/Users/guglielmo/Desktop/Università/Magistrale/Extra/WinterSchool/Tools/plink_mac_20191028/plink"

name = c("AD","CN","MCI")
for (i in name){
cmd <- sprintf(" --bfile ADNI3_PLINK_Final --keep %s_patients.txt --make-bed --out %s" 
               , i, i)
  system(paste0(path,cmd))
}

ADNI3 = read.table("ADNI3_PLINK_Final.bim",sep="\t")

# Sample snps
snps = paste0(ADNI3$V1,":",ADNI3$V4)
snps = snps[substr(snps,1,1)!="0"]
snps = unique(snps) # Facciamo unique perchè circa 10000 snps hanno posizioni duplicate (da verificare i nomi)

# Sample genes from GTF
gtfFile <- 'Gencode/gencode.v36lift37.annotation.gtf'
gtfFile = '../Gencode/exon_protein_coding.txt'
gtf <- read.table(gtfFile, header=F, stringsAsFactors=F, sep='\t') 
gtf.gene <- gtf[, c(1,4,5)]


gene.names <- unlist(lapply(gtf[, 9], function(x) {
  y <- strsplit(x, ';')[[1]][1] # prendo il 1 valore che sarebbe ENSEMBLE del gene, con 3 selezioni i SYMBOL
  gsub('gene_id ','', y)
}))

# prendiamo gli id degli esoni invece dei geni per non avere duplicati
#gene.exon <- unlist(lapply(gtf[gtf[,3]=="exon", 9], function(x) {
#  y = strsplit(x, "exon_id")
 # z = strsplit(y[[1]][2],";")[[1]][1]
#  z = gsub(" ","",z)
#}))

gtf.gene$Gene_id = gene.names
gtf.gene.unique = unique(gtf.gene)
gene_names_unique = unique(gene.names) 

gene.reduce = data.frame(V1=as.character("NA"),
                         V4=as.numeric("NA"),
                         V5=as.numeric("NA"),
                         Gene_id=as.character("NA"))


for (i in gene_names_unique){
  df = gtf.gene.unique[gtf.gene.unique$Gene_id==i,]
  In.ir <- IRanges(df$V4, df$V5)
  out.ir <- reduce(In.ir)
  df2 = cbind(V1=df$V1[1],as.data.frame(out.ir),Gene_id=i)
  colnames(df2) = c("V1","V4","V5","Width","Gene_id")
  gene.reduce = rbind(gene.reduce,df2[,-4])
}

gene.reduce = gene.reduce[-1,]
write.csv(gene.reduce,"gene.reduce.csv",quote=F,row.names = F)

gene.reduce = read.csv("../gene.reduce.csv", header = T)
#rownames(gtf.gene) <- gene.exon

# Define a few helper functions

string2range <- function(pos, delim=' ', region=TRUE) {
  posp <- as.data.frame(do.call(rbind, strsplit(pos, delim)))
  posp[,1] <- posp[,1]
  posp[,2] <- as.numeric(as.character(posp[,2]))
  if(region) {
    posp[,3] <- as.numeric(as.character(posp[,3]))
  } else {
    posp[,3] <- posp[,2]
  }
  return(posp)
}

range2GRanges <- function(df) {
  require(GenomicRanges)
  require(IRanges)
  gr <- GenomicRanges::GRanges(
    seqnames = df[,1],
    ranges=IRanges(start = df[,2], end = df[,3])
  )
  return(gr)
}

# convert SNPs to GRanges
snps.ranges <- string2range(snps, delim=":", region=FALSE)
head(snps.ranges)

snps.granges <- range2GRanges(snps.ranges)
seqlevelsStyle(snps.granges) <- "NCBI" 
names(snps.granges) <- snps
snps.granges

# convert genes to GRanges
row.names(gene.reduce) = paste0("exon", 1:nrow(gene.reduce))
gtf.granges <- range2GRanges(gene.reduce)
names(gtf.granges) <- row.names(gene.reduce) # usiamo i nome degli esoni indicizzati
seqlevelsStyle(gtf.granges) <- "NCBI" # leva il chr nei cromosomi lasciando solo il nome
gtf.granges = renameSeqlevels(gtf.granges,seqlevels(snps.granges)[-25]) 
gtf.granges

# Now that we have our two GRanges objects, we can easily overlap them using GenomicRanges::findOverlaps.
r1 <- snps.granges
r2 <- gtf.granges
overlap <- GenomicRanges::findOverlaps(r1, r2)
# make vector of SNPs to genes
hits <- names(r2)[subjectHits(overlap)]
names(hits) <- names(r1)[queryHits(overlap)]
hits
Rle(subjectHits(overlap))

length(hits)

library(tibble)
library(dplyr)
results = enframe(hits)
gene.reduce$Exon_id = rownames(gene.reduce)
merging = merge(results,gene.reduce,by.x="value",by.y="Exon_id")
merging2 = merging[,c(2,6)]
#ADNI4 = ADNI3[ADNI3 %in% merging2$name,]

#merging3 = merge(merging2, ADNI3[,c(2,7)],by.x="name",by.y="V7")
results2 = merging2 %>% group_by(Gene_id) %>% summarise(n = n(), SNPsNames = toString(unique(name)))

colnames(results2) = c("Genes","SNPsCount","SNPsNames")
results2$Genes <- gsub("\\..*","", results2$Genes)
write.table(results2,"../SnpToGene.txt",quote=F,row.names = F,sep="\t")

# ******* COUNTS SNPS PER CATEGORY PATIENTS # **************

# Reduce the recoded plink file to have just a set snps equal to the hits set 
ADNI3_reduced = ADNI3[ADNI3$V7 %in% names(hits),]
write.table(ADNI3_reduced$V2,"../SNPS_to_keep.txt",row.names = F,quote = F,col.names = F)

for (i in name){
  cmd <- sprintf(" --bfile %s --extract ../SNPS_to_keep.txt --recodeA --out Plink_reduced/%s_recoded" , i, i)
  system(paste0(path,cmd))
}

AD_recoded = read.table("Plink_reduced/AD_recoded.raw",sep=" ",header = T)
AD_recoded = AD_recoded[,-c(1,3:6)]
colnames(AD_recoded)[2:ncol(AD_recoded)] = unlist(lapply(colnames(AD_recoded)[2:ncol(AD_recoded)],
                                                  function(x) substr(x,0, nchar(x)-2)))  

results = read.table("SnpToGene.txt",sep="\t",header=T)
results$SNPsNames_Converter = "NA"

ADNI3 = read.table("Subjects/ADNI3_PLINK_Final.bim",sep="\t")
ADNI3$V7 = paste0(ADNI3$V1,":",ADNI3$V4)

for (i in 1:nrow(results)){
  snps = strsplit(results$SNPsNames[i],",")
  snps = lapply(snps,function(x) gsub(" ","",x))[[1]]
  ids = ADNI3[ADNI3$V7 %in% snps, ]
  #ids = ids[!duplicated(ids$V7),]
  results$SNPsNames_Converter[i] = paste(ids$V2,collapse=", ") 
}

write.table(results,"SnpToGene_RSID_2.0.txt",quote=F,row.names = F,sep="\t")

