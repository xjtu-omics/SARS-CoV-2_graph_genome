library(ggplot2)
library(data.table)
rm(list = ls())
allSNP <- fread("./SNP_anno_type", header = F, stringsAsFactors = F)
colnames(allSNP) <- c("ref", "position", "type", "numSeq", "anno_type", "acid")
snpLT100 <- allSNP[which(allSNP$numSeq > 100), ]
# snpLT100 <- rbind(snpLT100, data.frame(position = 241, type = "C->T", numSeq = 10, seqID = ""))
snpLT100 <- snpLT100[order(snpLT100$position), ]

sp <- ggplot(snpLT100, aes(x=position, y=numSeq, shape = anno_type, color=type)) + 
      geom_linerange(aes(x = position, ymin = 0, 
                         ymax = numSeq) , 
                     size = 0.5, colour = "gray", linetype = "dashed") +
      geom_point(size = 3.5, position="jitter") + 
      scale_color_manual(values=c("chocolate1", "chartreuse3", "brown3", "deepskyblue3",
                                  "darkorchid2", "gray33")) +
      xlab("position") + 
      ylab("number of sequence") + 
      scale_x_continuous(limits = c(1, 29903), breaks = c(1, 29903)) + 
      scale_y_continuous(limits = c(0, 1700)) +
      theme_classic() + 
      theme(axis.title = element_text(size = 12),
            axis.text = element_text(size = 12), 
            legend.text = element_text(size = 10),
            legend.position = "right")
print(sp)


seqInfo <- fread("./sequence_info", header = F, stringsAsFactors = F)
theSeqNumOfLocation <- as.data.frame(table(seqInfo$V2))
colnames(theSeqNumOfLocation) <- c("location", "number")
for(i in 1 : nrow(snpLT100)){
  oneSeqIDs <- snpLT100$seqID[i]
  theSplitSeqIDs <- strsplit(oneSeqIDs, "/")
  theEles <- theSplitSeqIDs[[1]]
  theLocation <- theEles[seq(from = 2, to = length(theEles), by = 3)]
  theLocationTable <- as.data.frame(table(theLocation))
  colnames(theLocationTable) <- c("location", paste(snpLT100$position[i], "_", snpLT100$type[i], "_", snpLT100$numSeq[i], sep = ''))
  theSeqNumOfLocation <- merge.data.frame(theSeqNumOfLocation, theLocationTable, by = "location", all = T)
}
theSeqNumOfLocation[which(is.na(theSeqNumOfLocation), 2)] <- 0
write.table(theSeqNumOfLocation, col.names = T, row.names = F, sep = "\t", quote = F, file = "./SNPLT100_diff_country.txt")


annoType <- fread("./SNP_anno_type", header = F, stringsAsFactors = F)
colnames(annoType) <- c("ref", "position", "type", "numSeq", "anno_type")
annoTypeFreq <- as.data.frame(table(annoType$anno_type))


hotSNPLocation <- read.table("./SNP_LT_100_diff_country.txt", header = F, stringsAsFactors = F, row.names = 1, sep='\t')
hotSNPLocation <- as.matrix(hotSNPLocation)
theSNPInfo <- hotSNPLocation[1, ]
mainlandInfo <- rownames(hotSNPLocation)[-1]
hotSNPLocation <- hotSNPLocation[-1, ]
hotSNPLocation <- matrix(as.numeric(hotSNPLocation), nrow = dim(hotSNPLocation)[1], byrow = F)
totalSeqNum <- colSums(hotSNPLocation)
theFC <- hotSNPLocation
theP <- hotSNPLocation
thePropation <- hotSNPLocation
theFCDF <- data.frame(FC = 0, pvalue = 0, mainland = "", sites = "")
theFCDF <- theFCDF[-1, ]
for(i in 1 : nrow(hotSNPLocation)){
  totalNum <- totalSeqNum[1]
  oneRegionNum <- hotSNPLocation[i, 1]
  for(j in 2 : ncol(hotSNPLocation)){
    oneSNPTotal <- totalSeqNum[j]
    oneSNPRegionNum <- hotSNPLocation[i, j]
    oneFC <- (oneSNPRegionNum / oneSNPTotal) / (oneRegionNum / totalNum)
    theFC[i, j] <- oneFC
    thePropation[i, j] <- oneSNPRegionNum / oneRegionNum
    x <- rbind(c(oneSNPRegionNum, oneRegionNum - oneSNPRegionNum), 
               c(oneSNPTotal - oneSNPRegionNum, (totalNum - oneRegionNum) - (oneSNPTotal - oneSNPRegionNum)))
    fishRes <- fisher.test(x, alternative = "greater")
    theP[i, j] <- fishRes$p.value
    theFCDF <- rbind(theFCDF, data.frame(FC = oneFC, pvalue = fishRes$p.value,
                                         mainland = mainlandInfo[i],
                                         sites = theSNPInfo[j])
                     )
  }
}

theLabel <- rep("", nrow(theFCDF))
theLabel[theFCDF$pvalue < 0.01] <- "*"
theLabel[theFCDF$pvalue < 0.0000001] <- "**"
theFCDF <- cbind(theFCDF, label = theLabel)

ggplot(data=theFCDF, aes(x=sites, y=FC, fill=mainland)) +
  geom_bar(stat="identity", position=position_dodge()) + 
  geom_hline(yintercept = 1, linetype = "dashed") + 
  scale_fill_manual(values=c("black", "chocolate1", 
                              "brown3", "deepskyblue3",
                              "darkorchid2", "chartreuse3")) +
  geom_text(aes(label=label), vjust=-0.1, color="black",
            position = position_dodge(0.9), size=3.5)+
  theme_classic() + 
  xlab("") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 13), 
        axis.text.x = element_text(angle = 60, hjust = 1),
        legend.text = element_text(size = 12),
        legend.position = "right")



theFC <- theFC[, -1]
rownames(theFC) <- mainlandInfo
colnames(theFC) <- theSNPInfo[-1]
theP  <- theP[, -1]
rownames(theP) <- mainlandInfo
colnames(theP) <- theSNPInfo[-1]
thePropation <- thePropation[, -1]
# thePBinary <- -log10(theP)
thePBinary <- theP
thePBinary[which(theP < 0.01, 2)] <- 1
thePBinary[which(theP >= 0.01, 2)] <- 0
library(pheatmap)

pheatmap(theFC, cluster_cols = F)

# thePBinary[which(thePBinary > 10, 2)] <- 10
rownames(thePBinary) <- mainlandInfo
colnames(thePBinary) <- theSNPInfo[-1]
pheatmap(thePBinary, cluster_cols = F, cluster_rows = F, 
         angle_col = 45, legend_breaks = c(0.25,0.75), 
         legend_labels = c(0, 1), 
         color = c("white", "red"), fontsize = 12, legend = F)

rownames(thePropation) <- mainlandInfo
colnames(thePropation) <- theSNPInfo[-1]
pheatmap(thePropation, cluster_cols = F, 
         cluster_rows = F, 
         display_numbers = TRUE,
         fontsize = 12, 
         fontsize_number = 12, 
         number_color = "black", 
         angle_col = 90, 
         color = colorRampPalette(c("white", "red"))(100)
)

#############################################
## all SNPs, date associations
#############################################
library(doBy)
library(ggplot2)
rm(list = ls())
theSNPWithDate <- fread("allSNP_non_summary", header = F, stringsAsFactors = F)
colnames(theSNPWithDate) <- c("position", "REF", "ALT", "ACC", "date")
SNPNumPerSeq <- as.data.frame(table(theSNPWithDate$ACC))
colnames(SNPNumPerSeq) <- c("ACC", "Number")
ACCWithDate <- data.frame(ACC = theSNPWithDate$ACC, date = theSNPWithDate$date)
ACCWithDate <- unique.data.frame(ACCWithDate)
SNPWithDateNumber <- merge.data.frame(x = SNPNumPerSeq, y = ACCWithDate, by = "ACC")
meanSNPPerDay <- summaryBy(Number ~ date, data = SNPWithDateNumber, FUN = "mean")
meanSNPPerDay <- meanSNPPerDay[c(-4,-89,-90,-91), ]
meanSNPPerDay <- meanSNPPerDay[order(meanSNPPerDay$date), ]
dates <- as.character(meanSNPPerDay$date)
meanSNPNumber <- data.frame(date = as.Date(dates), number = meanSNPPerDay$`Number."mean"`)
meanSNPNumber <- meanSNPNumber[order(meanSNPNumber$date), ]

sp <- ggplot(meanSNPNumber, aes(x=date, y=number)) + 
  geom_point(size = 3.5, position="jitter") + 
  stat_smooth(method=lm, level=0.95)+
  xlab("Collection Date") + 
  ylab("Mean number of SNPs") + 
  scale_x_date(labels = date_format("%m-%d-%Y"), breaks='5 days') +
  theme_classic() + 
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 13), 
        axis.text.x = element_text(angle = 60, hjust = 1))
print(sp)
theCor <- lm(number ~ date, data = meanSNPNumber)
summary(theCor)

###################### snp increase with date in different continent
library(doBy)
library(ggplot2)
rm(list = ls())
theSNPWithDate <- fread("allSNP_non_summary", header = F, stringsAsFactors = F)
colnames(theSNPWithDate) <- c("position", "REF", "ALT", "ACC", "date")
theSeqInfo <- fread("./sequence_info", header = F, stringsAsFactors = F)
locationMapping <- fread("./location_mapping", header = F, stringsAsFactors = F)
continents <- c()
for(i in 1 : nrow(theSNPWithDate)){
  oneacc <- theSNPWithDate$ACC[i]
  location <- theSeqInfo$V2[which(theSeqInfo$V5 == oneacc)]
  if(length(location) >= 2){
    message(paste(i, "some errors"))
  }
  continent <- locationMapping$V2[which(locationMapping$V1 == location)]
  continents <- c(continents, continent)
}
theSNPWithDate <- cbind(theSNPWithDate, continent = continents)
SNPNumPerSeq <- as.data.frame(table(theSNPWithDate$ACC))
colnames(SNPNumPerSeq) <- c("ACC", "Number")
ACCWithDate <- data.frame(ACC = theSNPWithDate$ACC, date = theSNPWithDate$date, continent = theSNPWithDate$continent)
ACCWithDate <- unique.data.frame(ACCWithDate)
SNPWithDateNumber <- merge.data.frame(x = SNPNumPerSeq, y = ACCWithDate, by = "ACC")

uniqueContinetns <- unique(continents)
for(oneContinent in uniqueContinetns){
  oneContinentData <- SNPWithDateNumber[which(SNPWithDateNumber$continent == oneContinent), ]
  meanSNPPerDay <- summaryBy(Number ~ date, data = oneContinentData, FUN = "mean")
  # meanSNPPerDay <- meanSNPPerDay[c(-4,-89,-90,-91), ]
  meanSNPPerDay <- meanSNPPerDay[order(meanSNPPerDay$date), ]
  dates <- as.character(meanSNPPerDay$date)
  meanSNPNumber <- data.frame(date = as.Date(dates), number = meanSNPPerDay$`Number."mean"`)
  meanSNPNumber <- meanSNPNumber[order(meanSNPNumber$date), ]
  
  sp <- ggplot(meanSNPNumber, aes(x=date, y=number)) + 
    ggtitle(oneContinent) +
    geom_point(size = 3.5, position="jitter") + 
    stat_smooth(method=lm, level=0.95)+
    xlab("Collection Date") + 
    ylab("Mean number of SNPs") + 
    scale_x_date(labels = date_format("%m-%d-%Y"), breaks='5 days') +
    theme_classic() + 
    theme(axis.title = element_text(size = 14),
          axis.text = element_text(size = 13), 
          axis.text.x = element_text(angle = 60, hjust = 1))
  print(sp)
  theCor <- lm(number ~ date, data = meanSNPNumber)
  summary(theCor)
}




library(data.table)
allMutation <- fread("./all_mutation_unique_LT2Seq.vcf", header = T, stringsAsFactors = F)
# allMutation <- allMutation[-1, ]
allMutation <- allMutation[order(allMutation$POS), ]
fwrite(allMutation, file = "./all_mutation_unique_LT2Seq_sort.vcf", row.names = F, col.names = T, sep = "\t", quote = F)


library(ggplot2)
theCatSnpSummary <- fread("../theCat_SNP_Summary.tab", header = F, stringsAsFactors = F)
SnpInfo <- theCatSnpSummary[, c(1, 2, 3)]
sp <- ggplot(SnpInfo, aes(x=V3)) + 
  geom_histogram(bins = 20, 
                 colour="black", 
                 fill="gray") +
  xlab("Frequency") + 
  ylab("Number of SNVs") + 
  scale_x_log10() + 
  theme_classic() + 
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 13))
print(sp)


##########################################
## the Anno snp summary
##########################################
library(data.table)
theSummaryData <- fread("./theCat_SNP_unique.anno.vcf.summary.table", header = F, stringsAsFactors = F)
colnames(theSummaryData) <- c("position", "type", "annoType", "location", "date")
theSummaryData <- theSummaryData[order(theSummaryData$position), ]
theSummaryData$date <- as.Date(theSummaryData$date)
locationMapping <- fread("./location_mapping", header = F, stringsAsFactors = F, sep = "\t")
colnames(locationMapping) <- c("location", "mainland")
theSummaryData <- cbind(theSummaryData, mainland = rep("", nrow(theSummaryData)))
for(i in 1 : nrow(locationMapping)){
  onelocation <- locationMapping$location[i]
  mainland    <- locationMapping$mainland[i]
  
  theSummaryData$mainland[which(theSummaryData$location == onelocation)] <- mainland
}

# positions <- theSummaryData$position

sp <- ggplot(theSummaryData, aes(x=date, y=position, color = annoType)) + 
  geom_point(size = 2, position="jitter", shape = 1) + 
  # stat_smooth(method=lm, level=0.95)+
  xlab("Collection Date") + 
  ylab("SNP position") + 
  # scale_x_continuous(breaks = c(1,2)) +
  theme_classic() + 
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 12), 
        axis.text.x = element_text(angle = 60),
        legend.text = element_text(size = 10),
        legend.position = "right")
print(sp)


sp <- ggplot(theSummaryData, aes(x=mainland, y=position, color = annoType)) + 
  geom_point(size = 2, position="jitter", shape = 1) + 
  # stat_smooth(method=lm, level=0.95)+
  xlab("") + 
  ylab("SNP position") + 
  # scale_x_continuous(breaks = c(1,2)) +
  theme_classic() + 
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 12), 
        axis.text.x = element_text(angle = 60, hjust = 1),
        legend.text = element_text(size = 10),
        legend.position = "right")
print(sp)



rm(list = ls())
s23404SeqInfo <- fread("./S_occured_seq", header = F, stringsAsFactors = F)
seqInDays <- as.data.frame(table(s23404SeqInfo$V6))
allSeqInfo <- fread("./sequence_info", header = F, stringsAsFactors = F)
allSeqInDays <- as.data.frame(table(allSeqInfo$V6))
colnames(seqInDays) <- c("date", "number_with_S_mut")
colnames(allSeqInDays) <- c("date", "all_number")
theSeqInDays <- merge.data.frame(x = allSeqInDays, y = seqInDays, by = "date", all = T)
theSeqInDays$number_with_S_mut[which(is.na(theSeqInDays$number_with_S_mut))] <- 0
theRatio <- theSeqInDays$number_with_S_mut / theSeqInDays$all_number
theSeqInDays <- cbind(theSeqInDays, ratio = theRatio)
theSeqInDays <- theSeqInDays[-c(4, 5, 31, 61), ]

theSeqInDays$date <- as.Date(theSeqInDays$date)
theSeqInDays <- theSeqInDays[order(theSeqInDays$date), ]

sp <- ggplot(theSeqInDays, aes(x=date, y=ratio)) + 
  geom_point(size = 3.5, position="jitter") + 
  # stat_smooth(method=lm, level=0.95)+
  xlab("Collection Date") + 
  ylab("Ratio of genome sequences with SNP \n at locus 23,403") + 
  scale_x_date(labels = date_format("%m-%d-%Y"), breaks='10 days') +
  theme_classic() + 
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 11), 
        axis.text.x = element_text(angle = 60, hjust = 1))
print(sp)



rm(list = ls())
s241SeqInfo <- fread("./241_SNP", header = F, stringsAsFactors = F)
seqInDays <- as.data.frame(table(s241SeqInfo$V6))
allSeqInfo <- fread("./sequence_info", header = F, stringsAsFactors = F)
allSeqInDays <- as.data.frame(table(allSeqInfo$V6))
colnames(seqInDays) <- c("date", "number_with_S_mut")
colnames(allSeqInDays) <- c("date", "all_number")
theSeqInDays <- merge.data.frame(x = allSeqInDays, y = seqInDays, by = "date", all = T)
theSeqInDays$number_with_S_mut[which(is.na(theSeqInDays$number_with_S_mut))] <- 0
theRatio <- theSeqInDays$number_with_S_mut / theSeqInDays$all_number
theSeqInDays <- cbind(theSeqInDays, ratio = theRatio)
theSeqInDays <- theSeqInDays[-c(4, 5, 31, 61), ]

theSeqInDays$date <- as.Date(theSeqInDays$date)
theSeqInDays <- theSeqInDays[order(theSeqInDays$date), ]

sp <- ggplot(theSeqInDays, aes(x=date, y=ratio)) + 
  geom_point(size = 3.5, position="jitter") + 
  # stat_smooth(method=lm, level=0.95)+
  xlab("Collection Date") + 
  ylab("Ratio of sequences with SNP at locus 241") + 
  scale_x_date(labels = date_format("%Y-%m-%d"), breaks='10 days') +
  theme_classic() + 
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 11), 
        axis.text.x = element_text(angle = 60, hjust = 1))
print(sp)


rm(list = ls())
library(scales)
s241SeqInfo <- fread("./SNP_14408_seqinfo", header = F, stringsAsFactors = F)
seqInDays <- as.data.frame(table(s241SeqInfo$V6))
allSeqInfo <- fread("./sequence_info", header = F, stringsAsFactors = F)
allSeqInDays <- as.data.frame(table(allSeqInfo$V6))
colnames(seqInDays) <- c("date", "number_with_S_mut")
colnames(allSeqInDays) <- c("date", "all_number")
theSeqInDays <- merge.data.frame(x = allSeqInDays, y = seqInDays, by = "date", all = T)
theSeqInDays$number_with_S_mut[which(is.na(theSeqInDays$number_with_S_mut))] <- 0
theRatio <- theSeqInDays$number_with_S_mut / theSeqInDays$all_number
theSeqInDays <- cbind(theSeqInDays, ratio = theRatio)
theSeqInDays <- theSeqInDays[-c(4, 5, 31, 61), ]

theSeqInDays$date <- as.Date(theSeqInDays$date)
theSeqInDays <- theSeqInDays[order(theSeqInDays$date), ]

theSeqInDays <- theSeqInDays[which(theSeqInDays$all_number >= 10), ]

sp <- ggplot(theSeqInDays, aes(x=date, y=ratio)) + 
  geom_point(size = 3.5, position="jitter") + 
  # stat_smooth(method=lm, level=0.95)+
  xlab("Collection Date") + 
  ylab("Ratio of sequences with SNP at locus 241") + 
  scale_x_date(labels = date_format("%Y-%m-%d"), breaks='10 days') +
  theme_classic() + 
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 11), 
        axis.text.x = element_text(angle = 60, hjust = 1))
print(sp)


#########################################################
#### random background of different SNV type
rm(list = ls())
setwd("~/2019-CnoV/SARS-COV-2_sequence/variance/random_annotation/")
library(data.table)
refPos2Type <- fread("../ref_genome.fasta.pos2type", header = F, stringsAsFactors = F, sep = "\t")
colnames(refPos2Type) <- c("id", "pos", "nuc")
snvTypeNum <- fread("./snv_type.txt", header = F, stringsAsFactors = F, sep = "\t")
colnames(snvTypeNum) <- c("ref", "alt", "num")

vcfHeader <-"##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description=\"All filters passed\">
##samtoolsVersion=1.9+htslib-1.9
##samtoolsCommand=samtools mpileup -uf /home/xiaofei/data/2019_ncov/variation/ref_genome.fasta sarsCov2.bam
##reference=file:///home/xiaofei/data/2019_ncov/variation/ref_genome.fasta
##contig=<ID=Wuhan-Hu-1,length=29903>
##ALT=<ID=*,Description=\"Represents allele(s) other than observed.\">
##INFO=<ID=INDEL,Number=0,Type=Flag,Description=\"Indicates that the variant is an INDEL.\">
##INFO=<ID=IDV,Number=1,Type=Integer,Description=\"Maximum number of reads supporting an indel\">
##INFO=<ID=IMF,Number=1,Type=Float,Description=\"Maximum fraction of reads supporting an indel\">
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Raw read depth\">
##INFO=<ID=VDB,Number=1,Type=Float,Description=\"Variant Distance Bias for filtering splice-site artefacts in RNA-seq data (bigger is better)\",Version=\"3\">
##INFO=<ID=RPB,Number=1,Type=Float,Description=\"Mann-Whitney U test of Read Position Bias (bigger is better)\">
##INFO=<ID=MQB,Number=1,Type=Float,Description=\"Mann-Whitney U test of Mapping Quality Bias (bigger is better)\">
##INFO=<ID=BQB,Number=1,Type=Float,Description=\"Mann-Whitney U test of Base Quality Bias (bigger is better)\">
##INFO=<ID=MQSB,Number=1,Type=Float,Description=\"Mann-Whitney U test of Mapping Quality vs Strand Bias (bigger is better)\">
##INFO=<ID=SGB,Number=1,Type=Float,Description=\"Segregation based metric.\">
##INFO=<ID=MQ0F,Number=1,Type=Float,Description=\"Fraction of MQ0 reads (smaller is better)\">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"List of Phred-scaled genotype likelihoods\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
##INFO=<ID=ICB,Number=1,Type=Float,Description=\"Inbreeding Coefficient Binomial test (bigger is better)\">
##INFO=<ID=HOB,Number=1,Type=Float,Description=\"Bias in the number of HOMs number (smaller is better)\">
##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count in genotypes for each ALT allele, in the same order as listed\">
##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">
##INFO=<ID=DP4,Number=4,Type=Integer,Description=\"Number of high-quality ref-forward , ref-reverse, alt-forward and alt-reverse bases\">
##INFO=<ID=MQ,Number=1,Type=Integer,Description=\"Average mapping quality\">
##INFO=<ID=SEQNUM,Number=1,Type=Integer,Description=\"SNP frequency\">
##INFO=<ID=sEQID,Number=1,Type=String,Description=\"Occured Seq IDs\">
##bcftools_callVersion=1.9+htslib-1.9
##bcftools_callCommand=call -mv; Date=Tue Apr  7 14:30:23 2020
##FILTER=<ID=LowQual,Description=\"Set if true: %QUAL<20 || DP>100\">
##bcftools_filterVersion=1.9+htslib-1.9
##bcftools_filterCommand=filter -s LowQual -e '%QUAL<20 || DP>100' sarsCov2.var.raw.vcf; Date=Tue Apr  7 14:30:23 2020
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sarsCov2.bam"
vcfOtherString <- "30.4183	PASS	DP=1;SGB=-0.379885;MQ0F=0;AC=2;AN=2;DP4=0,0,1,0;MQ=60;	GT:PL	1/1:60,3,0"

totalRandomNum <- 100
nuc <- c("A", "T", "C", "G")
for(i in 1 : totalRandomNum){
  message(i)
  write.table(vcfHeader, file = paste("theRandomSNVs/random_", i, ".vcf", sep = ''), append = F, quote = F, row.names = F, col.names = F)
  oneVCF <- data.frame(chr = "", pos = 0, id = ".", ref = "", alt = "", other = "")
  oneVCF <- oneVCF[-1, ]
  for(oneNuc in nuc){
    theSNVNum <- snvTypeNum[which(snvTypeNum$ref == oneNuc), ]
    thePositons <- refPos2Type[which(refPos2Type$nuc == oneNuc), ]
    theTotNum <- sum(theSNVNum$num)
    theRandomPos <- sample(thePositons$pos, size = theTotNum, replace = F)

    offset <- 0
    for(snvi in 1 : nrow(theSNVNum)){
      onePos <- theRandomPos[(1 + offset) : (theSNVNum$num[snvi] + offset)]
      offset <- theSNVNum$num[snvi] + offset
      chrs <- thePositons$id[which(thePositons$pos %in% onePos)]
      poss <- thePositons$pos[which(thePositons$pos %in% onePos)]
      id <- rep(".", length(onePos))
      ref <- rep(theSNVNum$ref[snvi], length(onePos))
      alt <- rep(theSNVNum$alt[snvi], length(onePos))
      
      oneVCF <- rbind(oneVCF, data.frame(chr = chrs, pos = poss, id = id, ref = ref, alt = alt, 
                                         other = rep(vcfOtherString, length(onePos))))
    }
  }
  oneVCF <- oneVCF[order(oneVCF$pos), ]
  write.table(oneVCF, file = paste("theRandomSNVs/random_", i, ".vcf", sep = ''), append = T, quote = F, row.names = F, col.names = F, sep = "\t")
}


theRandomTypeRes <- data.frame()
for(i in 1 : 100){
  message(i)
  oneRandomSNVType <- fread(paste("./theRandomSNVs/random_", i, ".anno.vcf.summary.table", sep = ''), header = F, stringsAsFactors = F)
  theRes <- as.data.frame(table(oneRandomSNVType$V4))
  if(nrow(theRandomTypeRes) == 0){
    theRandomTypeRes <- theRes
  }else{
    theRandomTypeRes <- merge.data.frame(x = theRandomTypeRes, y = theRes, by = "Var1", all = T)
  }
}
theRandomTypeRes[which(is.na(theRandomTypeRes), 2)] <- 0
theRandomNumbers <- theRandomTypeRes[, -1]
means <- rowMeans(as.matrix(theRandomNumbers))
sds   <- rowSds(as.matrix(theRandomNumbers))
theRandomRes <- data.frame(type = theRandomTypeRes$Var1, mean = means, sd = sds)



## random frequency background
rm(list = ls())
setwd("~/2019-CnoV/SARS-COV-2_sequence/variance/random_annotation/")
library(data.table)
totalSNVCount <- 18797
totalSeqNum <- 3160
refSeqLen <- 29303
allSNPNoSummary <- fread("../allSNP_non_summary", header = F, stringsAsFactors = F)
seqSNPNum <- as.data.frame(table(allSNPNoSummary$V4))

theRandomFreq <- data.frame(pos = 1 : refSeqLen)
for(i in 1 : 500){
  message(i)
  theRandomCount <- rep(0, refSeqLen)
  for(j in 1 : nrow(seqSNPNum)){
    oneSeqSNPNum <- seqSNPNum$Freq[j]
    oneRandomIndex <- sort(sample(1:refSeqLen, size = oneSeqSNPNum, replace = F))
    oneRandomCount <- rep(0, refSeqLen)
    oneRandomCount[oneRandomIndex] <- 1
    theRandomCount = theRandomCount + oneRandomCount
  }
  theRandomFreq <- cbind(theRandomFreq, theRandomCount)
}
theRandomFreq[which(theRandomFreq == 0, 2)] <- NA
theRandomFreqNum <- as.matrix(theRandomFreq[, -1])
means <- rowMeans(theRandomFreqNum, na.rm = T)
sds   <- rowSds(theRandomFreqNum, na.rm = T)
theRandPosFreq <- data.frame(pos = 1:refSeqLen, mean = means, sd = sds)
theRandPosFreq <- theRandomFreqNum[order(theRandomFreqNum$pos), ]
maxMean <- theRandPosFreq$mean[which.max(theRandPosFreq$mean)]
maxSD   <- theRandPosFreq$sd[which.max(theRandPosFreq$sd)]

library(ggplot2)
refSeqLen <- 29303
theCatSnpSummary <- fread("../theCat_SNP_Summary.tab", header = F, stringsAsFactors = F)
SnpInfo <- theCatSnpSummary[, c(1, 2, 3)]
SnpInfo <- SnpInfo[order(SnpInfo$V1), ]
pvalues <- 1 - pnorm(SnpInfo$V3, maxMean, maxSD)
SnpInfo <- cbind(SnpInfo, pvalues)

significantHighSNV <- SnpInfo[which(SnpInfo$pvalues < 0.0001), ]
colnames(significantHighSNV) <- c("position", "type", "num", "pvalues")
regionNum <- c()
for(i in 1 : 30){
  num <- length(significantHighSNV$position[which(significantHighSNV$position >= (i - 1) * 1000 & 
                                                    significantHighSNV$position <= i * 1000)])
  regionNum <- c(regionNum, num)
}

randomRegionNum <- matrix(0, nrow = 500, ncol = length(regionNum))
for(ri in 1 : 500){
  rpos <- sample(1:refSeqLen, size = nrow(significantHighSNV), replace = F)
  oneRRegionNum <- c()
  for(i in 1 : 30){
    num <- length(rpos[which(rpos >= (i - 1) * 1000 & 
                             rpos <= i * 1000)])
    oneRRegionNum <- c(oneRRegionNum, num)
  }
  randomRegionNum[ri, ] <- oneRRegionNum
}
randomRegionNumMean <- colMeans(randomRegionNum)
randomRegionNumSd   <- colSds(randomRegionNum)
pvalues <- 1 - pnorm(regionNum, randomRegionNumMean, randomRegionNumSd)

theResults <- data.frame(SNVRealNum = regionNum, background_mean = randomRegionNumMean, background_sd = randomRegionNumSd, pvalue = pvalues)
write.table(theResults, file = "bin_enriche_high_freq_SNV.txt", col.names = T, row.names = F, sep = "\t")



