library(foreach)
library(doMC)
registerDoMC()
library(multtest)
library(WGCNA)
library("org.Mm.eg.db")
library(biomaRt)
library(GOstats)
library("org.Mm.eg.db")
library("edgeR")
library(vegan)
library(ncf)

library(lawstat)

#source("http://bioconductor.org/biocLite.R")
# biocLite("org.Mm.eg.db")
 #biocLite("GOstats")

getDoParWorkers()
options(cores=27)
getDoParWorkers()

#enableWGCNAThreads(nThreads = 27)

setwd("/home/dan/workDir/alexShell")

source("/home/dan/workDir/functionDefinitions.R")
try(dir.create("resultsCoexpr"), silent = T)
try( dir.create("figuresCoexpr"), silent = F)


load("data/selectedData.RData")
sample_info=read.csv("data/sampleInfo.csv", header=T, row.names=1)

geneNames=rownames(selectedGeneCounts)

# divide the data in different groups
HSCC_H=selectedGeneCounts[,sample_info[,"Line"]=="H" ]
HSCC_L=selectedGeneCounts[,sample_info[,"Line"]=="L" ]

########################################################################################################################################
adjConsensus=adjacency(t(selectedGeneCounts), power=1)

powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold.fromSimilarity(adjConsensus, powerVector = powers, verbose = 5, moreNetworkConcepts=T)

plotNetConstruction(sft)
quartz.save("figuresCoexpr/netConstructionCoexpr.tif", type="tif", bg="white", dpi=300)
quartz.save("figuresCoexpr/netConstructionCoexpr.jpg", type="jpg", bg="white")

softPowerCoexpr=6
adjCoexpr=adjConsensus^softPowerCoexpr
adjCoexpr[is.na(adjCoexpr)]=0

connCoexpr=rowSums(adjCoexpr)
hist(connCoexpr, 100)


hierADJCoexpr = hclust(as.dist(1-adjCoexpr),method="average");

# code below might need modifications to make the number of modules between ~ 15 and 50, and to make number of grey genes less than ~ 2-3000
hybridCoexpr=cutreeHybrid(dendro = hierADJCoexpr, distM=1-adjCoexpr, cutHeight = 0.9995, minClusterSize = 100, deepSplit = 4, maxCoreScatter = NULL, minGap = NULL, maxAbsCoreScatter = NULL, minAbsGap = NULL, pamStage = TRUE, pamRespectsDendro = F, useMedoids = FALSE,  respectSmallClusters = TRUE, verbose = 2, indent = 0)

colorsCoexpr = labels2colors(hybridCoexpr$labels)
names(colorsCoexpr)=geneNames
table(colorsCoexpr)
length(table(colorsCoexpr))
modulesCoexpr=names(table(colorsCoexpr))
sum(colorsCoexpr=="grey")


fileConnSummary<-file("resultsCoexpr/SummaryResultsCoexpr.txt",  open="at")

writeLines(paste("Number modules ", length(table(colorsCoexpr)), sep=','), fileConnSummary)
writeLines(paste("Number grey genes  ",  sum(colorsCoexpr=="grey"), sep=','), fileConnSummary)

close(fileConnSummary)

adj_HSCC_H=adjacency(t(HSCC_H), power=softPowerCoexpr)
adj_HSCC_L=adjacency(t(HSCC_L), power=softPowerCoexpr)

names(colorsCoexpr)=geneNames

save(adjCoexpr, colorsCoexpr,modulesCoexpr, adj_HSCC_H,adj_HSCC_L, file="data/adjCoexprModules.RData")
load("data/adjCoexprModules.RData")
########################################################################################################


neuronsList=read.csv("data/CahoyNeurons.csv", header=TRUE)
neuronsSymbols= neuronsList[,"Gene.Name"]

astrosList=read.csv("data/CahoyAstros.csv", header=TRUE)
astrosSymbols= astrosList[,"Gene.Name"]

oligosList=read.csv("data/CahoyOligos.csv", header=TRUE)
oligosSymbols= oligosList[,"Gene.Name"]

moduleEnrichmentNeurons = moduleEnrichment (colorsCoexpr, neuronsSymbols)
moduleEnrichmentAstros = moduleEnrichment (colorsCoexpr, astrosSymbols)
moduleEnrichmentOligos = moduleEnrichment (colorsCoexpr, oligosSymbols)

cellTypeEnrichment=round(cbind(moduleEnrichmentNeurons,moduleEnrichmentAstros, moduleEnrichmentOligos),4)
colnames(cellTypeEnrichment)=c("Neurons", "Astros", "Oligos")
rownames(cellTypeEnrichment)=modulesCoexpr


# fileConnSummary<-file("resultsCoexpr/SummaryResultsCoexpr.txt",  open="at")
# writeLines(paste("Modules enriched in neuronal cell types\n "), fileConnSummary)
# close(fileConnSummary)

write.csv(cellTypeEnrichment, file="resultsCoexpr/cellTypeEnrich.csv", append=T)
close(fileConnSummary)


########################################################################################################

# consensus module preservation in individual networks
load("data/adjCoexprModules.RData")

multiData =vector("list",3)

multiData[[1]] =list(data= adjCoexpr)
multiData[[2]] =list(data= adj_HSCC_H)
multiData[[3]] =list(data= adj_HSCC_L)

names(multiData)=c("Consensus","HSCC_H", "HSCC_L")
checkSets(multiData, checkStructure = FALSE, useSets = NULL)

multiColor =vector("list",3)

multiColor[[1]] =as.vector(colorsCoexpr)
multiColor[[2]] =as.vector(colorsCoexpr)
multiColor[[3]] =as.vector(colorsCoexpr)

names(multiColor)=c("Consensus","HSCC_H", "HSCC_L")

modulePreservIndVsConsensus=modulePreservation(
  multiData,
  multiColor,
  dataIsExpr = F,
  networkType = "unsigned", 
  corFnc = "cor",
  corOptions = "use = 'p'",
  referenceNetworks = 1, 
  nPermutations = 200, 
  includekMEallInSummary = FALSE,
  restrictSummaryForGeneralNetworks = FALSE,
  calculateQvalue = FALSE,
  randomSeed = 12345, 
  maxGoldModuleSize = 1000, 
  maxModuleSize = 1000, 
  quickCor = 1, 
  ccTupletSize = 2, 
  calculateCor.kIMall = TRUE,
  useInterpolation = FALSE, 
  checkData = F, 
  greyName = "grey", 
  savePermutedStatistics = FALSE, 
  loadPermutedStatistics = FALSE, 
  permutedStatisticsFile = if (useInterpolation) "permutedStats-intrModules.RData" 
  else "permutedStats-actualModules.RData", 
  plotInterpolation = FALSE, 
  interpolationPlotFile = "modulePreservationInterpolationPlots.pdf", 
  discardInvalidOutput = TRUE,
  verbose = 3, indent = 0)


save(modulePreservIndVsConsensus, file="data/modulePreservation.RData")
load("data/modulePreservation.RData")

names(modulePreservIndVsConsensus)
names(modulePreservIndVsConsensus$preservation)
names(modulePreservIndVsConsensus$preservation$Z)
names(modulePreservIndVsConsensus$preservation$Z$ref.Consensus)
names(modulePreservIndVsConsensus$preservation$Z$ref.Consensus$inColumnsAlsoPresentIn.HSCC_H)


preservSummary=cbind(modulePreservIndVsConsensus$preservation$Z$ref.Consensus$inColumnsAlsoPresentIn.HSCC_H$moduleSize,
                     modulePreservIndVsConsensus$preservation$Z$ref.Consensus$inColumnsAlsoPresentIn.HSCC_H$Zsummary.pres,
                     modulePreservIndVsConsensus$preservation$Z$ref.Consensus$inColumnsAlsoPresentIn.HSCC_L$Zsummary.pres)

rownames(preservSummary)=rownames(modulePreservIndVsConsensus$preservation$Z$ref.Consensus$inColumnsAlsoPresentIn.HSCC_H)
colnames(preservSummary)=c("moduleSize", "Cons preserv in HSCC_H", "Cons preserv in HSCC_L")

write.csv(preservSummary, file="resultsCoexpr/coexprModulesPreserv.csv")
########################################################################################################
########################################################################################################
# save the results below for use with enrinchR
try(dir.create("resultsCoexpr/moduleGeneList"), silent = T)


coexprConnConsensus=intramodularConnectivity(adjCoexpr,  colorsCoexpr, scaleByMax=T)
#totalScaledConnectivity=coexprConnConsensus[,"kTotal"]/max(coexprConnConsensus[,"kTotal"])

coexprConnHigh=intramodularConnectivity(adj_HSCC_H,  colorsCoexpr, scaleByMax=T)
coexprConnLow=intramodularConnectivity(adj_HSCC_L,  colorsCoexpr, scaleByMax=T)

coexprResultsTable=cbind(colorsCoexpr, round(coexprConnConsensus[,"kWithin"],3), round(coexprConnHigh[,"kWithin"],3), round(coexprConnLow[,"kWithin"],3))
colnames(coexprResultsTable)=c("module", "consensus conn", "High Conn", "Low Conn")
rownames(coexprResultsTable)=geneNames

for (module in modulesCoexpr){
  print(module)
  currModuleInfo=cbind(rownames(coexprConnConsensus)[colorsCoexpr==module],round(coexprConnConsensus[coexprConnConsensus==module,"kWithin"],2))
  write.csv(currModuleInfo, file=paste("resultsCoexpr/moduleGeneList/module_", module, ".csv", sep=""), row.names=F, col.names=F)  
}
#############################################################################
# GO annotations

load("/home/dan/workDir/HDID2/data/transcriptInfoMouse.RData")
annotateMouseModulesGO(colorsCoexpr, transcriptInfoMouse, type="Coexpr")
##############################################################################
# record differential expression
resultsDEtotal=read.csv("resultsCoexpr/resultsDEtotal.csv", row.names=1)
summaryResultsShell=cbind(coexprResultsTable, resultsDEtotal[rownames(coexprResultsTable),1:3])

meanHSCC_H=rowMeans(HSCC_H)
meanHSCC_L=rowMeans(HSCC_L)

summaryResultsShell=cbind(summaryResultsShell, meanHSCC_H, meanHSCC_L)

##############################################################################
# find differentially variable genes

pvalVar=rep(1, length(geneNames))
names(pvalVar)=geneNames
pvalVar=pvalVar

for (gene in geneNames){
  pvalVar[gene]=var.test(x=HSCC_H[gene,], y=HSCC_L[gene,])$p.value
}

adjpOut=mt.rawp2adjp(pvalVar, proc=c( "BH"))
fdrVar=adjpOut$adjp[order(adjpOut$index),2]

sdHSCC_H=apply(HSCC_H, 1, sd)
sdHSCC_L=apply(HSCC_L, 1, sd)

summaryResultsShell=cbind(summaryResultsShell, pvalVar, fdrVar, sdHSCC_H, sdHSCC_L)

write.csv(summaryResultsShell, file="resultsCoexpr/summaryResultsShell.csv")
summaryResultsShell=read.csv("resultsCoexpr/summaryResultsShell.csv")

##################################################################################33

deGenes=geneNames[summaryResultsShell[,"PValue"] < 0.01]
dvGenes=geneNames[summaryResultsShell[,"pvalVar"] < 0.01]

modulesEnrichDEDV=moduleEnrichment(colorsCoexpr, union(deGenes, dvGenes))
affectedModulesDEDV=names(modulesEnrichDEDV)[modulesEnrichDEDV<(0.05/length(modulesCoexpr))]

fileConnSummary<-file("resultsCoexpr/SummaryResultsCoexpr.txt",  open="at")
writeLines(paste("\n "), fileConnSummary)

writeLines(paste("Modules affected by DE and/or DV expression changes ", affectedModulesDEDV, sep=','), fileConnSummary)
close(fileConnSummary)

############################################################################
# evaluate changes in edge strength
rawAdj1=adjacency(t(HSCC_H), power=1)
rawAdj2=adjacency(t(HSCC_L), power=1)

diffEdgesShell = diffEdges(rawAdj1, rawAdj2, n1=dim(HSCC_H)[2], n2=dim(HSCC_L)[2], pThreshold=0.01, adjThreshold=0.5, nCores=7)
  
save(diffEdgesShell, file="resultsCoexpr/diffEdgesShell.RData")
load("resultsCoexpr/diffEdgesShell.RData")

totalEdges=(length(geneNames))^2
affectedEdges=sum(diffEdgesShell)
edgeChangeRate=affectedEdges/totalEdges

geneChangeEdgeCount=rowSums(diffEdgesShell)
names(geneChangeEdgeCount)=geneNames

pValuesEdgeChange=rep(1,length(geneNames))
names(pValuesEdgeChange)=geneNames

for (gene in geneNames){
  pValuesEdgeChange[gene]=binom.test(x=geneChangeEdgeCount[gene], n=length(geneNames), p=edgeChangeRate, alternative  ="g")$p.value
}

adjpOut=mt.rawp2adjp(pValuesEdgeChange, proc=c( "BH"))
fdrEdgesChange=adjpOut$adjp[order(adjpOut$index),2]

summaryResultsShell=cbind(summaryResultsShell, pValuesEdgeChange, fdrEdgesChange,geneChangeEdgeCount )
write.csv(summaryResultsShell, file="resultsCoexpr/summaryResultsShell.csv")


genesChangedEdges=geneNames[pValuesEdgeChange < 0.01]
geneChangeEdgeCount[genesChangedEdges]
mean(geneChangeEdgeCount)
mean(geneChangeEdgeCount[genesChangedEdges])

median(geneChangeEdgeCount)
median(geneChangeEdgeCount[genesChangedEdges])


hist(geneChangeEdgeCount)

names(colorsCoexpr)=geneNames
edgesModuleEnrich=moduleEnrichment(colorsCoexpr, genesChangedEdges)
affectedModules=names(edgesModuleEnrich[edgesModuleEnrich<(0.05/length(modulesCoexpr))])

fileConnSummary<-file("resultsCoexpr/SummaryResultsCoexpr.txt",  open="at")
writeLines(paste("Modules affected by edge changes  ", affectedModules, sep=','), fileConnSummary)
close(fileConnSummary)


##############################################################################
summaryResultsShell=read.csv("resultsCoexpr/summaryResultsShell.csv")
# summaryResultsShell[,3:18]=round(summaryResultsShell[,3:18], 3)
# write.csv(summaryResultsShell, file="resultsCoexpr/summaryResultsShell.csv")


hubsResultsShell=summaryResultsShell[(summaryResultsShell[,"de.calls"] ==1 | summaryResultsShell[,"fdrVar"]<0.1 | summaryResultsShell[,"fdrEdgesChange"]<0.1) & summaryResultsShell[,"moduleScaledConn"] > 0.8, ]
write.csv(hubsResultsShell, file="resultsCoexpr/hubsResultsShell.csv")
##############################################################################
#calculate overlap with Mulligan 2006

MulliganTable=read.csv("data/SITable2_MulliganPonamarev_2006.csv", skip=1, header=T)
MulliganGenes=as.character(MulliganTable[,1])

genesShellPval=summaryResultsShell[(summaryResultsShell[,"PValue"] <0.01 | summaryResultsShell[,"pvalVar"]<0.01 | summaryResultsShell[,"pValuesEdgeChange"]<0.01),"Gene.Names"]
genesShellFdr=summaryResultsShell[(summaryResultsShell[,"fdrDE"] <0.01 | summaryResultsShell[,"fdrVar"]<0.01 | summaryResultsShell[,"fdrEdgesChange"]<0.01),"Gene.Names"]

intersect(MulliganGenes,genesShellPval)
length(intersect(MulliganGenes,genesShellFdr))
# length(intersect(MulliganGenes,genesShellFdr))
# [1] 252
fisher.test(summaryResultsShell[,"Gene.Names"] %in% genesShellPval, summaryResultsShell[,"Gene.Names"]  %in% MulliganGenes, alternative = "g")
fisher.test(summaryResultsShell[,"Gene.Names"] %in% genesShellFdr, summaryResultsShell[,"Gene.Names"]  %in% MulliganGenes, alternative = "g")

MulliganDE=summaryResultsShell[,"Gene.Names"]  %in% MulliganGenes

summaryResultsShell=cbind(summaryResultsShell, MulliganDE)
write.csv(summaryResultsShell, file="resultsCoexpr/summaryResultsShell.csv")

##############################################################################
hubsResultsShellCoexpr=read.csv("resultsCoexpr/hubsResultsShell.csv")[,"Gene.Names"]
hubsResultsShellCoSplicEx=read.csv("resultsCoSplicEx/hubsResultsShell.csv")[,"Gene.Names"]

intersect(hubsResultsShellCoexpr, hubsResultsShellCoSplicEx)

##############################################################################
# summaryResultsShell=read.csv("resultsCoexpr/summaryResultsCoexprShell.csv")
# 
# summaryResultsShell=cbind(summaryResultsShell[,2:4], coexprConnHigh[,"kWithin"], coexprConnLow[,"kWithin"], summaryResultsShell[,6:18])
# summaryResultsShell[,3:18]=round(summaryResultsShell[,3:18], 3)
# colnames(summaryResultsShell)[1]="Gene Names"
# colnames(summaryResultsShell)[4:5]=c("High Conn", "Low Conn")
# write.csv(summaryResultsShell, file="resultsCoexpr/summaryResultsShell.csv")

##########################################################################################################
hubsResultsShell=read.csv("resultsCoexpr/hubsResultsShell.csv")

colorsHubsAffected=rep("white", length(geneNames))
names(colorsHubsAffected)=geneNames
colorsHubsAffected[hubsResultsShell[,"Gene.Names"]]="black"

# GO annotations

load("/home/dan/workDir/HDID2/data/transcriptInfoMouse.RData")
annotateMouseModulesGO(colorsHubsAffected, transcriptInfoMouse, type="CoexprHubsAffected")

