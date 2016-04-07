
###############################################################################################################################################

plotNetConstruction = function (sft)

{
sizeGrWindow(15, 9)
par(mfrow = c(2,3));
par(mar=c(5, 5, 3, 3))
par(lab=c(5, 7, 3))

cex1 = 1.5;

plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Fit",
     main = paste("Scale independence"), font.lab=2, cex.lab=1.2, xaxt = "n", yaxt = "n", col="red");

x.at <- powers
axis(1, at = x.at, labels = as.character(powers), cex.lab=1.5, font=2, lwd=2, lwd.ticks=2)

y.at <- c(seq(from = 0, to=1, by=.1))
axis(2, at = y.at, labels = format(y.at,digits=1), cex.lab=1.5, font=2, lwd=2, lwd.ticks=2)


abline(h=0.75,col="black")

mtext("A", NORTH<-3, at=-2, line=0.25, cex=2.5, font=2)
############################################################################

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity",type="l",
     main = paste("Mean Connectivity"),lwd=2, font.lab=2, cex.lab=1.2, xaxt = "n", yaxt = "n", col="red");


x.at <- powers
axis(1, at = x.at, labels = as.character(powers), cex.lab=1.5, font=2, lwd=2, lwd.ticks=2)

y_range=c(0,max(sft$fitIndices[,5]))

y_step=(y_range[2]-y_range[1])/7

y.at <- c(seq(from = y_range[1], to=y_range[2], by=y_step))
axis(2, at = y.at, labels = format(y.at,digits=1), cex.lab=1.5, font=2, lwd=2, lwd.ticks=2)

mtext("B", NORTH<-3, at=-2, line=0.25, cex=2.5, font=2)
############################################################################

# Density as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,8],
     xlab="Soft Threshold (power)",ylab="Density",type="l",
     main = paste("Density"),lwd=2, font.lab=2, cex.lab=1.2, xaxt = "n", yaxt = "n", col="red");


x.at <- powers
axis(1, at = x.at, labels = as.character(powers), cex.lab=1.5, font=2, lwd=2, lwd.ticks=2)

y_range=c(0,max(sft$fitIndices[,8]))

y_step=(y_range[2]-y_range[1])/5

y.at <- c(seq(from = y_range[1], to=y_range[2], by=y_step))
axis(2, at = y.at, labels = format(y.at,digits=1), cex.lab=1.5, font=2, lwd=2, lwd.ticks=2)

mtext("C", NORTH<-3, at=-2, line=0.25, cex=2.5, font=2)
############################################################################

# Centralization as a function of the soft-thresholding power
y_range=c(0,max(sft$fitIndices[,9]))

plot(sft$fitIndices[,1], sft$fitIndices[,9],
     xlab="Soft Threshold (power)",ylab="Centralization",type="l", ylim=y_range,
     main = paste("Centralization"),lwd=2, font.lab=2, cex.lab=1.2, xaxt = "n", yaxt = "n", col="red");


x.at <- powers
axis(1, at = x.at, labels = as.character(powers), cex.lab=1.5, font=2, lwd=2, lwd.ticks=2)

y_range=c(0,max(sft$fitIndices[,9]))

y_step=(y_range[2]-y_range[1])/5

y.at <- c(seq(from = y_range[1], to=y_range[2], by=y_step))
axis(2, at = y.at, labels = format(y.at,digits=1), cex.lab=1.5, font=2, lwd=2, lwd.ticks=2)

mtext("D", NORTH<-3, at=-2, line=0.25, cex=2.5, font=2)
############################################################################

# Heterogeneity as a function of the soft-thresholding power
y_range=c(0,max(sft$fitIndices[,10]))

plot(sft$fitIndices[,1], sft$fitIndices[,10],
     xlab="Soft Threshold (power)",ylab="Centralization",type="l", ylim=y_range,
     main = paste("Heterogeneity"),lwd=2, font.lab=2, cex.lab=1.2, xaxt = "n", yaxt = "n", col="red");

x.at <- powers
axis(1, at = x.at, labels = as.character(powers), cex.lab=1.5, font=2, lwd=2, lwd.ticks=2)

y_range=c(0,max(sft$fitIndices[,10]))

y_step=(y_range[2]-y_range[1])/5

y.at <- c(seq(from = y_range[1], to=y_range[2], by=y_step))
axis(2, at = y.at, labels = format(y.at,digits=1), cex.lab=1.5, font=2, lwd=2, lwd.ticks=2)

mtext("E", NORTH<-3, at=-2, line=0.25, cex=2.5, font=2)

}

####################################################################################

annotateMouseModulesGO = function (colorsCoexpr, transcriptInfoMouse, type)
  # this assumes the gene identifiers are gene symbols, whichr are the names of the colorsCoexpr
  
{
  library("GOstats")
  
  dirName=paste("results", type, sep="")
  try(dir.create(dirName), silent = F)
  
  
  fileName=paste("results", type,"/modules", type, "GO.csv", sep="")
  modules=names(table(colorsCoexpr))
  geneNames=names(colorsCoexpr)    
  networkTranscriptInfo= transcriptInfoMouse[which(transcriptInfoMouse[,"external_gene_name"] %in% geneNames),]
  
  pValue=0.01/length(modules)
  univ= unique(as.character(networkTranscriptInfo[,"entrezgene"]))
  
  univ=unique(univ)
  
  for (module in modules){
    print(module)
    geneNamesModule=geneNames[colorsCoexpr==module]
    
    moduleEntrez = unique(as.character(networkTranscriptInfo[which(networkTranscriptInfo[,"external_gene_name"] %in% geneNamesModule),"entrezgene"]))
    BPparams <- new("GOHyperGParams", geneIds = moduleEntrez, universeGeneIds = univ,  annotation = "org.Mm.eg.db",ontology = "BP", pvalueCutoff = pValue, conditional = T, testDirection = "over")
        
    BPres <- hyperGTest(BPparams)
    res=summary(BPres)
    print(res)
    res=cbind(rep(paste(module,"module BP"), dim(res)[1]),res)
    colnames(res)[1]="module and GO term"
    write.table(res, fileName, append=TRUE, sep =",", quote=F, row.names=F)
    
    
    CCparams <- new("GOHyperGParams", geneIds = moduleEntrez,
                    universeGeneIds = univ, annotation = "org.Mm.eg.db", 
                    ontology = "CC", pvalueCutoff = pValue, conditional = F,
                    testDirection = "over")
    
    CCres <- hyperGTest(CCparams)
    res=summary(CCres)
    res=cbind(rep.int(paste(module,"module CC"), dim(res)[1]),res)
    colnames(res)[1]="module and GO term"
    write.table(res, fileName, append=TRUE, sep =",", quote=F, row.names=F)
    
    print(res)
    
    MFparams <- new("GOHyperGParams", geneIds = moduleEntrez,
                    universeGeneIds = univ, annotation = "org.Mm.eg.db", 
                    ontology = "MF", pvalueCutoff = pValue, conditional = F,
                    testDirection = "over")
    
    MFres <- hyperGTest(MFparams)
    res=summary(MFres)
    res=cbind(rep.int(paste(module,"module MF"), dim(res)[1]),res)
    colnames(res)[1]="module and GO term"
    print(res)
    
    write.table(res, fileName, append=TRUE, sep =",", quote=F, row.names=F)
       
  }
    
}

##########################################################################################
differentialExpression = function (geneReads, groupFactor)
  
{
  d=DGEList(counts= geneReads, group= groupFactor)
  
  d <- estimateTagwiseDisp(d)
  de.tgw <- exactTest(d, dispersion="tagwise") 

  de.calls <- decideTestsDGE(de.tgw, p=0.05, adjust="BH")


  return(as.list(c(de.tgw, de.calls)))
}

##########################################################################################

differentialVariability = function (geneReads, groupFactor, nCores)
  
{
  library(car)

  leveneOut=apply(geneReads, 1, leveneTest, geneReads, groupFactor)
  
  
  
  return(as.list(c(de.tgw, de.calls)))
}


##########################################################################################
#cbind(HDID_M, HSNPT_M),exonGeneNameSelected, geneNames, phenotypeVectorM, nPerm=0
# exonReads=cbind(HDID_M, HSNPT_M)
# exonGeneNames=exonGeneNameSelected
# phenotypeVector=phenotypeVectorM
# nPerm=0

splicingSignificance = function (exonReads, exonGeneNames, geneNames, phenotypeVector, nPerm, nCores)
  
{
  options(cores=nCores)
  
  mantelGSS=rep(0,length(geneNames))
  names(mantelGSS)=geneNames
  mantelGSS_P= mantelGSS
  
  phenotypeDistances =stats::dist(phenotypeVector, method="manhattan")
  nSamples=dim(exonReads)[2] 
  sampleNames=colnames(exonReads)
  mantelList=foreach (gene_name = geneNames, .inorder=T, .verbose=T) %dopar% {
    
    currGeneExonCounts= t(exonReads[which(exonGeneNames==gene_name),]) # distance is computed between rows of the matrix
    
    
#     totalNumberReads=rep(0, length(geneNames))
#     names(totalNumberReads)=geneNames
#    for (gene_name in geneNames){
      
#     print(geneName)
#     print(paste("Numer exons is :", as.character(dim(currGeneExonCounts)[2]), sep=""))
#     totalNumberReads[gene_name]=sum(currGeneExonCounts)
#     print(sum(currGeneExonCounts)) 
#    }
        
    currGeneExonCounts= t(exonReads[which(exonGeneNames==gene_name),]) # distance is computed between rows of the matrix
    if (sum(currGeneExonCounts) > 0){
      exon_numbers =dim(currGeneExonCounts)[2]
      #print(exon_numbers)
      
      # compute significance of each exon individually
      exonSignif=0*(1:exon_numbers)
      for (exonIdx in 1:exon_numbers){
        currExonCounts= currGeneExonCounts[,exonIdx]
        if (sum(currExonCounts)>0){
          currExonDist=dist(currExonCounts, method="canberra")
          exonSignif[exonIdx]=vegan::mantel(phenotypeDistances, currExonDist, permutations=0, na.rm=T)$statistic
        }else{
          
          exonSignif[exonIdx]=0
        }
      }
      exonSignif[is.na(exonSignif)]=0
      
      # rank the exons in decreased order of significance
      rankedExons=length(exonSignif)-rank(exonSignif,ties.method="first")+1
      
      # record the cumulative significance of exon groups
      groupExonSignif=rep(0, exon_numbers)
      
      sequenceDistances=array(data=0,dim=c(exon_numbers,nSamples,nSamples),dimnames=list(1:exon_numbers, sampleNames, sampleNames))
      
      for (idx in 1:exon_numbers){
        exonIndexes=which(rankedExons <=idx)
        cumulativeDist= (1/idx)*as.matrix(dist(currGeneExonCounts[,exonIndexes], method="canberra"))
        cumulativeDist[is.na(cumulativeDist)]=0		
        groupExonSignif[idx]=vegan::mantel(phenotypeDistances, as.dist(cumulativeDist), permutations=0, na.rm = T)$statistic	
        sequenceDistances[idx,,]= cumulativeDist
      }
      
      mantelGSS[gene_name]=max(groupExonSignif)	
      max_idx=which(groupExonSignif==max(groupExonSignif, na.rm=T))
      max_idx= max_idx[1]
      
      mantelResults=vegan::mantel(phenotypeDistances, as.dist(sequenceDistances[max_idx,,]), permutations=nPerm, na.rm=TRUE)
      mantelResults
    } 
    # finished selecting the group of exons that maximizes splicing significance           
  }
  
  for (i in 1:length(geneNames)){
    if (!is.null(mantelList[[i]])){
    mantelGSS[i]=mantelList[[i]]$statistic
    mantelGSS_P[i]=mantelList[[i]]$signif	
    }
  }
  
  mantelGSS_P[is.na(mantelGSS_P)]=1
  mantelGSS[is.na(mantelGSS)]=0
  
  adjpOut=mt.rawp2adjp(mantelGSS_P, proc=c( "BH"))
  fdrValues=adjpOut$adjp[order(adjpOut$index),2]
  

  returnStruct=cbind(mantelGSS, mantelGSS_P, fdrValues)
  colnames(returnStruct)=c("splicing signif", "p value", "FDR value")
  return(returnStruct)

  #  return(mantelGSS)
  
}

##########################################################################################



#######################################################################
#########################################################################################
#######################################################################
# rawAdj1=adjacency(t(HSCC_H), power=1)
# rawAdj2=adjacency(t(HSCC_L), power=1)
# n1=dim(HSCC_H)[2]
# n2=dim(HSCC_L)[2]
# pThreshold=0.01
# adjThreshold=0.5
# nCores=7
# 
# diffEdgesShell = diffEdges(rawAdj1, rawAdj2, n1=dim(HSCC_H)[2], n2=dim(HSCC_L)[2], pThreshold=0.01, adjThreshold=0.5, nCores=14)



diffEdges = function (rawAdj1, rawAdj2, n1, n2, pThreshold=0.01, adjThreshold=0.5, nCores){
  library(psych)

  options(cores=nCores)
  geneNames=rownames(rawAdj1)
  nGenes=length(geneNames)

  diffEdgesMatrix <- foreach (gene=geneNames, .inorder=T, .verbose = T, .combine=cbind) %dopar% {
  
#   diffEdges=0*rawAdj1
#   for (gene in geneNames){ 
  
    vectorEdges=rep(1, nGenes)
    names(vectorEdges)=geneNames
    
    vector1=rawAdj1[gene,]
    vector2=rawAdj2[gene,]
    idxsVector=abs(vector1 -vector2) > adjThreshold

    for (gene2 in geneNames[idxsVector==T]){
      vectorEdges[gene2]=r.test(n=n1, r12=rawAdj1[gene,gene2], r34=rawAdj2[gene,gene2], n2=n2)$p
    }
    
    vectorEdges<pThreshold 
    #diffEdges[gene,]=vectorEdges<pThreshold
  }
  
  return(diffEdgesMatrix)
}
##########################################################################################

#rawMantelAdj1=, rawMantelADj2, canberraListSelected, samples1, samples2, adjThreshold=0.5, pThreshold=0.01, nCores=5, nboot=100


#######################################################################

diffEdgesMantel = function (canberraListSelected, samples1, samples2, adjThreshold=0.5, pThreshold=0.01, nCores=28 , nboot=200){
  library(psych)
  
  options(cores=nCores)
  geneNames=names(canberraListSelected)
  nGenes=length(geneNames)
  nSamples1=length(samples1)
  nSamples2=length(samples2)
  
  canberraListSelected1=vector(mode="list", length=nGenes)
  names(canberraListSelected1)=geneNames
  canberraListSelected2=canberraListSelected1
  
  for(gene in geneNames){
    canberraListSelected1[[gene]]=canberraListSelected[[gene]][samples1,samples1]
    canberraListSelected2[[gene]]=canberraListSelected[[gene]][samples2,samples2]
  }
  
  # reformat the data so one can use WGCNA adjacency function to construct CoSplicEx adjacency matrix
  lengthVector=length(as.vector(as.dist(canberraListSelected1[[1]])))
  distData=matrix(data=0, nrow=lengthVector, ncol=length(canberraListSelected1))
  colnames(distData)=names(canberraListSelected1)
  
  for(gene in names(canberraListSelected1)) {
    distData[,gene]=as.vector(as.dist(canberraListSelected1[[gene]]))
  } 
  adjCoSplicEx1=adjacency(distData,power=1)
  
  lengthVector=length(as.vector(as.dist(canberraListSelected2[[1]])))
  distData=matrix(data=0, nrow=lengthVector, ncol=length(canberraListSelected2))
  colnames(distData)=names(canberraListSelected2)
  
  for(gene in names(canberraListSelected2)) {
    distData[,gene]=as.vector(as.dist(canberraListSelected2[[gene]]))
  } 
  adjCoSplicEx2=adjacency(distData,power=1)
  
  
  #diffEdges=0*adjCoSplicEx1
  
  diffEdges=foreach (gene1=geneNames, .inorder=T,  .verbose = T, .combine=cbind) %dopar% {
    
    vector1=adjCoSplicEx1[gene1,]
    vector2=adjCoSplicEx2[gene1,]
    idxsVector=abs(vector1 -vector2) > adjThreshold
    
    diffEdgesCurr=rep(0,nGenes)
    names(diffEdgesCurr)=geneNames
    for (gene2 in geneNames[idxsVector]){

     mantelBoots1=rep(0, nboot)
     for (i in 1:nboot) {
        
        nBootResamples=round(.9*nSamples1)
        samples=sample(samples1, nBootResamples, replace = F)
        mantelBoots1[i]=vegan::mantel(canberraListSelected1[[gene1]][samples,samples], canberraListSelected1[[gene2]][samples,samples], permutations = 0 )$statistic
        
      }
      CI1=c( quantile(mantelBoots1, p=0.5*pThreshold), quantile(mantelBoots1, p=1-0.5*pThreshold))
 
      mantelBoots2=rep(0, nboot)
      for (i in 1:nboot) {
        nBootResamples=round(.9*nSamples2)
        samples=sample(samples2, nBootResamples, replace = F)
        mantelBoots1[i]=vegan::mantel(canberraListSelected2[[gene1]][samples,samples], canberraListSelected2[[gene2]][samples,samples], permutations = 0 )$statistic
      }
      CI2=c( quantile(mantelBoots2, p=0.5*pThreshold), quantile(mantelBoots2, p=1-0.5*pThreshold))
      
      lowerCIS=c(CI1[1], CI2[1])
      upperCIS=c(CI1[2], CI2[2])
      
      minUpper=min(upperCIS)
      maxLower=max(lowerCIS)
            
      diffEdgesCurr[gene2]=minUpper<maxLower
    }
    diffEdgesCurr 
  }
  return(diffEdges)
}
##########################################################################################


##########################################################################################
diffVarSplicing = function (canberraListSelected, phenotypeVector, nPerm, nCores)
{
  
  library(parallel)
  library(vegan)
  geneNames=names(canberraListSelected)
  samples=names(phenotypeVector)
  

  canberraListMatrices=foreach(gene=geneNames,.inorder=T,.verbose = F) %dopar% {
    
    as.matrix(canberraListSelected[[gene]])[samples, samples]

  } 
  names(canberraListMatrices)=geneNames
  
  mrppPvalues=rep(1, length(geneNames))
  names(mrppPvalues)=geneNames
  
  for (gene in geneNames) {
    print(which(geneNames==gene))
    if (sum(is.na(canberraListMatrices[[gene]]))==0){
    mrppPvalues[gene]=mrpp(canberraListMatrices[[gene]], grouping=phenotypeVector, permutations=nPerm, parallel=nCores)$Pvalue
    }
  }
    
  adjpOut=mt.rawp2adjp(mrppPvalues, proc=c( "BH"))
  mrppFdrValues=adjpOut$adjp[order(adjpOut$index),2]
  
  returnValue=cbind(mrppPvalues, mrppFdrValues)
  return(returnValue)
}

##########################################################################################


moduleEnrichment = function (colors, signifGenes)
  
{ 
  geneNames=names(colors)
  modules=names(table(colors))
  modulePvaluesEnrich=rep(1,length(modules))
  names(modulePvaluesEnrich)=modules
  for (module in modules){
    moduleGenes=geneNames[colors==module]
    modulePvaluesEnrich[module]=fisher.test(geneNames %in% moduleGenes, geneNames %in% signifGenes, alternative = "g")$p.value
  }
  return(modulePvaluesEnrich)
}


#############################################################################################
##########################################################################################

