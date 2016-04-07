library(foreach)
library(doMC)
registerDoMC()
library(WGCNA)
library(multtest)

getDoParWorkers()
options(cores=27)
getDoParWorkers()

setwd("/home/dan/workDir/alexShell")
source("/home/dan/workDir/functionDefinitions.R")

load("data/selectedData.RData")
geneNames=rownames(adjCoSplicEx)
############################################################################
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold.fromSimilarity(adjCoSplicEx, powerVector = powers, verbose = 5, moreNetworkConcepts=T)

plotNetConstruction(sft)
quartz.save("figures/netConstructionCoSplicEx.tif", type="tif", bg="white", dpi=300)
quartz.save("figures/netConstructionCoSplicEx.jpg", type="jpg", bg="white")

softPowerCoSplicEx=6
adjCoSplicEx=adjCoSplicEx^softPowerCoSplicEx
adjCoSplicEx[is.na(adjCoSplicEx)]=0
hierADJCoSplicEx = hclust(as.dist(1-adjCoSplicEx),method="average");

# code below might need modifications to make the number of modules between ~ 15 and 50, and to make number of grey genes less than ~ 2-3000
hybridCoSplicEx=cutreeHybrid(dendro = hierADJCoSplicEx, distM=1-adjCoSplicEx, cutHeight = 0.99995, maxCoreScatter=0.995, minClusterSize = 20, deepSplit = 4, minGap = NULL, maxAbsCoreScatter = NULL, minAbsGap = NULL, 
                             pamStage = T, maxPamDist=1,pamRespectsDendro = F, useMedoids = FALSE,  respectSmallClusters = F, verbose = 2, indent = 0)

colorsCoSplicEx = labels2colors(hybridCoSplicEx$labels)
names(colorsCoSplicEx)=geneNames
table(colorsCoSplicEx)
length(table(colorsCoSplicEx))
modulesCoSplicEx=names(table(colorsCoSplicEx))
sum(colorsCoSplicEx=="grey")
modulesCoSplicEx=names(table(colorsCoSplicEx))

save(adjCoSplicEx, colorsCoSplicEx,modulesCoSplicEx,geneNames, file="data/adjModulesCoSplicEx.RData")
load("data/adjModulesCoSplicEx.RData")
names(colorsCoSplicEx)=geneNames
modulesCoSplicEx=names(table(colorsCoSplicEx))

fileConnSummary<-file("resultsCoSplicEx/SummaryResultsCoSplicEx.txt",  open="at")

writeLines(paste("Number modules ", length(table(colorsCoSplicEx)), sep=','), fileConnSummary)
writeLines(paste("Number grey genes  ",  sum(colorsCoSplicEx=="grey"), sep=','), fileConnSummary)



########################################################################################################
# save the results below for use with enrinchR

try(dir.create("resultsCoSplicEx/moduleGeneList"), silent = T)

CoSplicExConn=intramodularConnectivity(adjCoSplicEx, colorsCoSplicEx, scaleByMax=T)
totalScaledConnectivity=CoSplicExConn[,"kTotal"]/max(CoSplicExConn[,"kTotal"])
CoSplicExConn=cbind(CoSplicExConn, totalScaledConnectivity)

CoSplicExConnHigh=intramodularConnectivity(adjCoSplicEx_High,  colorsCoSplicEx, scaleByMax=T)
CoSplicExConnLow=intramodularConnectivity(adjCoSplicEx_Low,  colorsCoSplicEx, scaleByMax=T)



for (module in modulesCoSplicEx){
  print(module)
   currModuleInfo=cbind(rownames(CoSplicExConn)[colorsCoSplicEx==module],CoSplicExConn[colorsCoSplicEx==module,"kWithin"])
  write.csv(currModuleInfo, file=paste("resultsCoSplicEx/moduleGeneList/module_", module, ".csv", sep=""), row.names=F, col.names=F)  
  }
#############################################################################
names(colorsCoSplicEx)=geneNames


neuronsList=read.csv("data/CahoyNeurons.csv", header=TRUE)
neuronsSymbols= neuronsList[,"Gene.Name"]

astrosList=read.csv("data/CahoyAstros.csv", header=TRUE)
astrosSymbols= astrosList[,"Gene.Name"]

oligosList=read.csv("data/CahoyOligos.csv", header=TRUE)
oligosSymbols= oligosList[,"Gene.Name"]

moduleEnrichmentNeurons = moduleEnrichment (colorsCoSplicEx, neuronsSymbols)
moduleEnrichmentAstros = moduleEnrichment (colorsCoSplicEx, astrosSymbols)
moduleEnrichmentOligos = moduleEnrichment (colorsCoSplicEx, oligosSymbols)

cellTypeEnrichment=round(cbind(moduleEnrichmentNeurons,moduleEnrichmentAstros, moduleEnrichmentOligos),4)
colnames(cellTypeEnrichment)=c("Neurons", "Astros", "Oligos")
rownames(cellTypeEnrichment)=modulesCoSplicEx


fileConnSummary<-file("resultsCoSplicEx/SummaryResultsCoSplicEx.txt",  open="at")
writeLines(paste("\n "), fileConnSummary)

writeLines(paste("Modules enriched in neuronal cell types ", cbind(moduleEnrichmentNeurons,moduleEnrichmentAstros, moduleEnrichmentOligos), sep=','), fileConnSummary)
close(fileConnSummary)


################################################################################3
# GO annotations

load("/home/dan/workDir/HDID2/data/transcriptInfoMouse.RData")
annotateMouseModulesGO(colorsCoSplicEx, transcriptInfoMouse, "CoSplicEx")

##############################################################################
########################################################################################################

# consensus module preservation in individual networks
load("data/adjModulesCoSplicEx.RData")

multiData =vector("list",3)

multiData[[1]] =list(data= adjCoSplicEx)
multiData[[2]] =list(data= adjCoSplicEx_High)
multiData[[3]] =list(data= adjCoSplicEx_Low)

names(multiData)=c("Consensus","HSCC_H", "HSCC_L")
checkSets(multiData, checkStructure = FALSE, useSets = NULL)

multiColor =vector("list",3)

multiColor[[1]] =as.vector(colorsCoSplicEx)
multiColor[[2]] =as.vector(colorsCoSplicEx)
multiColor[[3]] =as.vector(colorsCoSplicEx)

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
  randomSeed = 45, 
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


save(modulePreservIndVsConsensus, file="resultsCoSplicEx/modulePreservation.RData")
load("resultsCoSplicEx/modulePreservation.RData")

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

write.csv(preservSummary, file="resultsCoSplicEx/CoSplicExModulesPreserv.csv")




# evaluate splicing significance
sample_info=read.csv("data/sampleInfo.csv", header=T, row.names=1)
# divide the data in different groups
HSCC_H=selectedExonCounts[,sample_info[,"Line"]=="H" ]
HSCC_L=selectedExonCounts[,sample_info[,"Line"]=="L" ]

# evaluate differentially spliced genes
phenotypeVector=c(rep(1,dim(HSCC_H)[2]),rep(0,dim(HSCC_L)[2]))
names(phenotypeVector)=c(colnames(HSCC_H), colnames(HSCC_L))
geneSplicingSignif=splicingSignificance(cbind(HSCC_H, HSCC_L),exonGeneNameSelected, geneNames, phenotypeVector, nPerm=1000, nCores=28)

save(geneSplicingSignif, file="resultsCoSplicEx/splicingSignifResults.RData")
load("resultsCoSplicEx/splicingSignifResults.RData")
# # ##############################################################################

pValues=geneSplicingSignif[ ,"p value"]
names(pValues)=geneNames
dsplicedGenes=names(pValues)[pValues<0.01]

############################################################################################
# evaluate genes with high splicing variability

diffVarSplicing = diffVarSplicing(canberraListSelected, phenotypeVector, nPerm=500, nCores=28) 

save(diffVarSplicing,file="resultsCoSplicEx/diffVarSplicing.RData")
load("resultsCoSplicEx/diffVarSplicing.RData")

############################################################################################
# save summary differential splicing and differential splicing variability

summaryResults=cbind(colorsCoSplicEx, CoSplicExConn[,c("kWithin", "totalScaledConnectivity")], geneSplicingSignif, diffVarSplicing)
colnames(summaryResults)=c("module", "moduleScaledConn", "netScaledConn", "splicingSignif", "splicingPval", "splicingFDR", "splicingVarPval", "splicingVarFDR")
write.csv(summaryResults, file="resultsCoSplicEx/summaryResults.csv")

###########################################################################################

pValuesSplicing=summaryResults[ ,"splicingPval"]
names(pValuesSplicing)=geneNames
dSplicedGenes=names(pValuesSplicing)[pValuesSplicing<0.01]

pValuesSplicingVar=summaryResults[ ,"splicingVarPval"]
names(pValuesSplicingVar)=geneNames
dSplicedVarGenes=names(pValuesSplicingVar)[pValuesSplicingVar<0.01]


DSDVModuleEnrich=moduleEnrichment(colorsCoSplicEx, union(dSplicedGenes, dSplicedVarGenes))
affectedModulesDSDV=names(modulesCoSplicEx)[DSDVModuleEnrich<(0.05/length(modulesCoSplicEx))]

fileConnSummary<-file("resultsCoSplicEx/SummaryResultsCoSplicEx.txt",  open="at")
writeLines(paste("Modules overlap with DS/DSVAR genes ", DSDVModuleEnrich), fileConnSummary)
close(fileConnSummary)
# ##################################################################

#############################################################################################################333
samples1=colnames(HSCC_H)
samples2=colnames(HSCC_L)

diffEdgesMantel = diffEdgesMantel(canberraListSelected, samples1, samples2, adjThreshold=0.5, pThreshold=0.01, nCores=28, nboot=200)
save(diffEdgesMantel, file="resultsCoSplicEx/diffEdgesMantel.RData")

load("resultsCoSplicEx/diffEdgesMantel.RData")

totalEdges=(length(geneNames))^2
affectedEdges=sum(diffEdgesMantel)
edgeChangeRate=affectedEdges/totalEdges

geneChangeEdgeCount=rowSums(diffEdgesMantel)
names(geneChangeEdgeCount)=geneNames

pValuesEdgeChange=rep(1,length(geneNames))
names(pValuesEdgeChange)=geneNames

for (gene in geneNames){
  pValuesEdgeChange[gene]=binom.test(x=geneChangeEdgeCount[gene], n=length(geneNames), p=edgeChangeRate, alternative  ="g")$p.value
}

adjpOut=mt.rawp2adjp(pValuesEdgeChange, proc=c( "BH"))
fdrEdgesChange=adjpOut$adjp[order(adjpOut$index),2]

summaryResultsShell=read.csv("resultsCoSplicEx/summaryResultsShell.csv")

summaryResultsShell=cbind(summaryResultsShell, pValuesEdgeChange, fdrEdgesChange,geneChangeEdgeCount )
summaryResultsShell[,3:12]=round(summaryResultsShell[,3:12], 3)

write.csv(summaryResultsShell, file="resultsCoSplicEx/summaryResultsShell.csv")


genesChangedEdges=geneNames[pValuesEdgeChange < 0.01]
geneChangeEdgeCount[genesChangedEdges]
mean(geneChangeEdgeCount)
mean(geneChangeEdgeCount[genesChangedEdges])

median(geneChangeEdgeCount)
median(geneChangeEdgeCount[genesChangedEdges])


hist(geneChangeEdgeCount)

names(colorsCoSplicEx)=geneNames
edgesModuleEnrich=moduleEnrichment(colorsCoSplicEx, genesChangedEdges)
affectedModules=names(edgesModuleEnrich[edgesModuleEnrich<(0.05/length(modulesCoSplicEx))])

fileConnSummary<-file("resultsCoSplicEx/SummaryResultsCoSplicEx.txt",  open="at")
writeLines(paste("Modules affected by edge changes  ", affectedModules, sep=','), fileConnSummary)
close(fileConnSummary)
###########################################################################################################

summaryResultsShell=read.csv("resultsCoSplicEx/summaryResultsShell.csv")

hubsResultsShell=summaryResultsShell[(summaryResultsShell[,"splicingFDR"]<0.1 | summaryResultsShell[,"splicingVarFDR"]<0.1 | summaryResultsShell[,"fdrEdgesChange"]<0.1) & summaryResultsShell[,"moduleScaledConn"] > 0.8, ]
write.csv(hubsResultsShell, file="resultsCoSplicEx/hubsResultsShell.csv")
###########################################################################################################
# summaryResultsShell=read.csv("resultsCoSplicEx/summaryResultsShell.csv")
# 
# summaryResultsShell=cbind(summaryResultsShell[,1:3], CoSplicExConnHigh[,"kWithin"], CoSplicExConnLow[,"kWithin"], summaryResultsShell[,6:13])
# summaryResultsShell[,3:13]=round(summaryResultsShell[,3:13], 3)
# colnames(summaryResultsShell)[1]="Gene Names"
# colnames(summaryResultsShell)[4:5]=c("High Conn", "Low Conn")
# write.csv(summaryResultsShell, file="resultsCoSplicEx/summaryResultsShell.csv")

##################################################################################################################################################33

hubsResultsShell=read.csv("resultsCoSplicEx/hubsResultsShell.csv")

colorsHubsAffected=rep("white", length(geneNames))
names(colorsHubsAffected)=geneNames
colorsHubsAffected[hubsResultsShell[,"Gene.Names"]]="black"

# GO annotations

load("/home/dan/workDir/HDID2/data/transcriptInfoMouse.RData")
annotateMouseModulesGO(colorsHubsAffected, transcriptInfoMouse, type="CoSplicExHubsAffected")

