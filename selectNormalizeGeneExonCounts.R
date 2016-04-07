      
library(edgeR)
library(WGCNA)
library(biomaRt)
library(plotrix)

library(foreach)
library(doMC)
registerDoMC()
library(proxy)

getDoParWorkers()
options(cores=14)
getDoParWorkers()

setwd("/home/dan/workDir/alexShell")

source("/home/dan/workDir/functionDefinitions.R")
try(dir.create("resultsCoexpr"), silent = T)
try( dir.create("figuresCoexpr"), silent = F)
try(dir.create("resultsCoSplicEx"), silent = T)
try( dir.create("figuresCoSplicEx"), silent = F)

geneReadsRaw=as.matrix(read.csv("data/RNASeq019_SH_mm10_gene_reads_not_normalized.csv", header=F, row.names=1))
geneColnames=as.matrix(read.csv("data/geneColnames.csv", header=F))
colnames(geneReadsRaw)=geneColnames

sample_info=read.csv("data/sampleInfo.csv", header=T, row.names=1)
#colnames(geneReadsRaw)= colnames(sample_info)



# divide the data in different groups
HSCC_H=geneReadsRaw[,sample_info[,"Line"]=="H" ]
HSCC_L=geneReadsRaw[,sample_info[,"Line"]=="L" ]

groupSelection=c(rep("HSCC_H",dim(HSCC_H)[2]),rep("HSCC_L",dim(HSCC_L)[2]))
groupSelection =factor(groupSelection)

d=DGEList(counts= cbind(HSCC_H, HSCC_L), group= groupSelection)
d <- estimateTagwiseDisp(d)
de.tgw <- exactTest(d, dispersion="tagwise") 
de.calls <- decideTestsDGE(de.tgw, p=0.05, adjust="BH")
sum(de.calls)
resultsDEtotal=cbind(de.tgw$table, de.calls)

# select genes with logCPM > 0 for further inclusion in network construction
resultsDEtotal=resultsDEtotal[resultsDEtotal[,"logCPM"]>0,]
write.csv(resultsDEtotal, file="resultsCoexpr/resultsDEtotal.csv")
resultsDEtotal=read.csv("resultsCoexpr/resultsDEtotal.csv")
rownames(resultsDEtotal)=resultsDEtotal[,1]

###################################################################################

# calculate edgeR normalization factors and normalize the data - use all data not just selected
UQnormFactors=calcNormFactors(geneReadsRaw, method=c("upperquartile"))

effectiveLibrarySizes= UQnormFactors*colSums(geneReadsRaw)
meanEffLibSize=mean(effectiveLibrarySizes)
countNormFactor= meanEffLibSize/effectiveLibrarySizes

normalizedGeneCountsUQ=0* geneReadsRaw

for (sample in 1:dim(normalizedGeneCountsUQ)[2]){	
	normalizedGeneCountsUQ[,sample]= geneReadsRaw[, sample]* countNormFactor[sample]	
}

geneReads=normalizedGeneCountsUQ[rownames(resultsDEtotal), ]
##################################################
# # find genes with high number of reads

connCounts=softConnectivity(t(geneReads), power=6)

quantileConn=quantile(connCounts, seq(0, 1, 0.1))  
geneReadsHighConn=geneReads[connCounts>quantileConn[6],]
geneNamesHighConn=rownames(geneReadsHighConn)
###################################################################################
exonCounts=read.csv("data/RNASeq019_SH_mm10_exon_reads_not_normalized.csv", header=F, row.names=1)
exonColnames=as.matrix(read.csv("data/exonColnames.csv", header=F))
colnames(exonCounts)=exonColnames

splitIDs=mapply(strsplit, rownames(exonCounts), MoreArgs=list(split="_", fixed = FALSE, perl = FALSE, useBytes = FALSE))
exonGeneName=unlist(lapply(splitIDs, "[[", 4))
exon_start=unlist(lapply(splitIDs, "[[", 2))

sum(exonCounts[exonGeneName=="Drd2",])
sum(geneReadsHighConn[geneNamesHighConn=="Drd2",])

exon_unique_id=mapply(paste, exonGeneName, exon_start, MoreArgs=list(sep="_"))
names(exonGeneName)=exon_unique_id

rownames(exonCounts)= exon_unique_id
sampleNames=as.vector(colnames(exonCounts))

normExonCounts=0* exonCounts
for (sample in 1:dim(exonCounts)[2]){
  normExonCounts[,sample]= exonCounts[, sample]* countNormFactor[sample]	
}
normExonCounts =round(normExonCounts)
rownames(normExonCounts)=rownames(exonCounts)


# select exons from genes with at least 1 CPM
exonCountsHighCounts=normExonCounts[which(exonGeneName %in% rownames(resultsDEtotal)),]
exonGeneNamesHighCounts=exonGeneName[exonGeneName %in% rownames(resultsDEtotal)]
  
canberraListExons=foreach (geneName = rownames(resultsDEtotal), .inorder=T, .verbose = T) %dopar% {
  #geneName=geneNames[i]
  currExonCounts= exonCountsHighCounts[which(exonGeneNamesHighCounts==geneName),]	
  if (is.null(dim(currExonCounts))){
    exonDistMatrix=as.matrix(dist(as.matrix(currExonCounts), method="canberra"))
  } else {
    exonDistMatrix=as.matrix(dist(t(as.matrix(currExonCounts)), method="canberra"))
  }
  colnames(exonDistMatrix)=exonColnames
  rownames(exonDistMatrix)=exonColnames
  exonDistMatrix
  
}
names(canberraListExons)=rownames(resultsDEtotal)

########################################################################################################


########################################################################################################

save(canberraListExons, file="data/canberraListExons.RData")
load("data/canberraListExons.RData")

nGenes=length(canberraListExons)
gene_indexes=1:nGenes

# reformat the data so one can use WGCNA adjacency function to construct CoSplicEx adjacency matrix
lengthVector=length(as.vector(as.dist(canberraListExons[[1]])))
distData=matrix(data=0, nrow=lengthVector, ncol=length(canberraListExons))
colnames(distData)=names(canberraListExons)

for(gene in names(canberraListExons)) {
  distData[,gene]=as.vector(as.dist(canberraListExons[[gene]]))
} 

adjCoSplicEx_large=adjacency(distData,power=6)


save(adjCoSplicEx_large, file="data/adjCoSplicEx_large.RData")

load("data/adjCoSplicEx_large.RData")
#just in case ...
adjCoSplicEx_large[is.na(adjCoSplicEx_large)]=0
diag(adjCoSplicEx_large)=1
colnames(adjCoSplicEx_large)=rownames(adjCoSplicEx_large)


connCoSplicEx=rowSums(adjCoSplicEx_large)
quantileConnExons=quantile(connCoSplicEx, probs = seq(0, 1, 0.1))  

geneNamesHighCoSplicExConn=names(canberraListExons)[connCoSplicEx>quantileConnExons[6]]
################################################################################################3
selectedGeneCounts=geneReadsHighConn

canberraListSelected=canberraListExons[geneNamesHighCoSplicExConn]
adjCoSplicEx=adjCoSplicEx_large[geneNamesHighCoSplicExConn,geneNamesHighCoSplicExConn]

exonGeneNameSelected=geneNamesHighCoSplicExConn
selectedExonCounts=normExonCounts[which(exonGeneName %in% geneNamesHighCoSplicExConn),]

########################################################################################################
canberraListExonsHigh=foreach (geneName = rownames(adjCoSplicEx), .inorder=T, .verbose = T) %dopar% {
  #geneName=geneNames[i]
  
  currExonCounts= exonCountsHighCounts[which(exonGeneNamesHighCounts==geneName),sample_info[,"Line"]=="H"]	
  if (is.null(dim(currExonCounts))){
    exonDistMatrix=as.matrix(dist(as.matrix(currExonCounts), method="canberra"))
  } else {
    exonDistMatrix=as.matrix(dist(t(as.matrix(currExonCounts)), method="canberra"))
  }
#   colnames(exonDistMatrix)=exonColnames
#   rownames(exonDistMatrix)=exonColnames
  exonDistMatrix
  
}
names(canberraListExonsHigh)=rownames(adjCoSplicEx)
nGenes=length(canberraListExonsHigh)
gene_indexes=1:nGenes

# reformat the data so one can use WGCNA adjacency function to construct CoSplicEx adjacency matrix
lengthVector=length(as.vector(as.dist(canberraListExonsHigh[[1]])))
distData=matrix(data=0, nrow=lengthVector, ncol=length(canberraListExonsHigh))
colnames(distData)=names(canberraListExonsHigh)

for(gene in names(canberraListExonsHigh)) {
  distData[,gene]=as.vector(as.dist(canberraListExonsHigh[[gene]]))
} 

adjCoSplicEx_High=adjacency(distData,power=6)
adjCoSplicEx_High[is.na(adjCoSplicEx_High)]=0

canberraListExonsLow=foreach (geneName = rownames(adjCoSplicEx), .inorder=T, .verbose = T) %dopar% {
  #geneName=geneNames[i]
  currExonCounts= exonCountsHighCounts[which(exonGeneNamesHighCounts==geneName),sample_info[,"Line"]=="L"]	
  if (is.null(dim(currExonCounts))){
    exonDistMatrix=as.matrix(dist(as.matrix(currExonCounts), method="canberra"))
  } else {
    exonDistMatrix=as.matrix(dist(t(as.matrix(currExonCounts)), method="canberra"))
  }
#   colnames(exonDistMatrix)=exonColnames
#   rownames(exonDistMatrix)=exonColnames
  exonDistMatrix
  
}
names(canberraListExonsLow)=rownames(adjCoSplicEx)
nGenes=length(canberraListExonsLow)
gene_indexes=1:nGenes

# reformat the data so one can use WGCNA adjacency function to construct CoSplicEx adjacency matrix
lengthVector=length(as.vector(as.dist(canberraListExonsLow[[1]])))
distData=matrix(data=0, nrow=lengthVector, ncol=length(canberraListExonsLow))
colnames(distData)=names(canberraListExonsLow)

for(gene in names(canberraListExonsLow)) {
  distData[,gene]=as.vector(as.dist(canberraListExonsLow[[gene]]))
} 

adjCoSplicEx_Low=adjacency(distData,power=6)
adjCoSplicEx_Low[is.na(adjCoSplicEx_Low)]=0


########################################################################################################
########################################################################################################




save(selectedGeneCounts, canberraListSelected,adjCoSplicEx,adjCoSplicEx_Low,adjCoSplicEx_High,selectedExonCounts, exonGeneNameSelected,file="data/selectedData.RData")
load("data/selectedData.RData")



############################################################################################
