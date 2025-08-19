rm(list=ls())
exprSetTPM=read.table('38002DEG_TPM.txt',header=T,sep="\t",row.names=1)
dim(exprSetTPM)

library(RColorBrewer)
library(WGCNA)
library(flashClust)
options(stringsAsFactors = FALSE)
allowWGCNAThreads()
dim(exprSetTPM)
expr=exprSetTPM
dim(expr)

srr57=read.table('sample.txt',header=T,sep="\t")
coldata1<-data.frame(srr=srr57$Sample,sample = srr57$Sample)
coldata=coldata1
head(coldata)
dim(coldata)

multiExpr = vector(mode = "list", length = 1)
multiExpr[[1]] = list(data = t(expr));
checkSets(multiExpr)

nSets=checkSets(multiExpr)$nSets

filterMultiExpr<-function(multiExpr,nSets)
{ gsg = goodSamplesGenesMS(multiExpr, verbose = 3);
print(gsg$allOK)
if (!gsg$allOK)
{
  # Print information about the removed genes:
  if (sum(!gsg$goodGenes) > 0)
    printFlush(paste("Removing genes:", paste(names(multiExpr[[1]]$data)[!gsg$goodGenes], collapse = ", ")))
  for (set in 1:nSets)
  {
    if (sum(!gsg$goodSamples[[set]]))
      printFlush(paste("In set", setLabels[set], "removing samples",paste(rownames(multiExpr[[set]]$data)[!gsg$goodSamples[[set]]], collapse = ", ")))
    # Remove the offending genes and samples
    multiExpr[[set]]$data = multiExpr[[set]]$data[gsg$goodSamples[[set]], gsg$goodGenes];
  } }
print(checkSets(multiExpr))
return(multiExpr)
}

multiExpr<-filterMultiExpr(multiExpr,nSets)
shortLabels="All"


powers = c(c(1:10), seq(from = 12, to=40, by=2))
for(type in c("signed","unsigned")){
  for(corU in c("cor","bicor")){
    # Initialize a list to hold the results of scale-free analysis
    powerTables = vector(mode = "list", length = nSets);
    # Call the network topology analysis function for each set in turn
    for(set in 1:nSets){
      powerTables[[set]] = list(data = pickSoftThreshold(multiExpr[[set]]$data, powerVector=powers, verbose = 2, corFnc = get(corU), networkType = type, blockSize=10000)[[2]])      }
    collectGarbage()
    # Plot the results:
    colors=brewer.pal(nSets,"Set1")
    # Will plot these columns of the returned scale free analysis tables
    plotCols = c(2,5,6,7)
    colNames = c("Scale Free Topology Model Fit", "Mean connectivity", "Median connectivity", "Max connectivity");
    # Get the minima and maxima of the plotted points
    ylim = matrix(NA, nrow = 2, ncol = 4);
    for (set in 1:nSets)
    {
      for (col in 1:length(plotCols))
      {
        ylim[1, col] = min(ylim[1, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE);
        ylim[2, col] = max(ylim[2, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE);
      }
    }
    # Plot the quantities in the chosen columns vs. the soft thresholding power
    # sizeGrWindow(8, 6)
    pdf(paste0("wgcna.choosePower.",gsub(" ","",type),".",corU,".pdf") )
    par(mfcol = c(2,2));
    par(mar = c(4.2, 4.2 , 2.2, 0.5))
    cex1 = 0.7;
    for (col in 1:length(plotCols)) for (set in 1:nSets)
    {
      if (set==1)
      {
        plot(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2], xlab="Soft Threshold (power)",ylab=colNames[col],type="n", ylim = ylim[, col], main = colNames[col]);
        addGrid()
      }
      if (col==1)
      {
        text(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
             labels=powers,cex=cex1,col=colors[set]);
      } else
        text(powerTables[[set]]$data[,1], powerTables[[set]]$data[,plotCols[col]], labels=powers,cex=cex1,col=colors[set]);
      if (col==1)
      {
        legend("bottomright", legend = shortLabels, col = colors, pch = 20) ;
      } else
        legend("topright", legend = shortLabels, col = colors, pch = 20) ;
    }
    dev.off()
    assign(paste0("powerTables.",gsub(" ","",type),".",corU) , powerTables)
  }
}
# examine plot and also print out power choice, plus default 12
powerTables=list()

for(type in c("unsigned", "signed", "signed hybrid"))
{
  for(corU in c("cor","bicor")){
    p=paste0("powerTables.",gsub(" ","",type),".",corU)
    print(p)
    print(lapply(get(p), function(x){ x<-x$data; x$Power[x$SFT.R.sq>0.8 & x$slope<0][1]} ))
    powerTables[[paste0(gsub(" ","",type),".",corU)]] <-get(p)[[1]]
  }
}

powers=c(12,14,24)
save(expr, multiExpr, powerTables, powers, file = "wgcna.salt.prep.Rdata")

#Server running
for(corM in c("bicor","pearson")){
  for(power in powers){
    # construct A2D5 network
    cgn =  blockwiseModules(
      # Input data
      multiExpr[[1]]$data,
      # Data checking options
      checkMissingData = TRUE,
      # Options for splitting data into blocks
      maxBlockSize =  30000,  # as individual network # 5000 for 4G memory, 20000 for 16G, 30000 for 32 G
      #randomSeed = 12345,
      # Network construction arguments: correlation options, use bicor instead of default pearson
      corType = corM,
      # Adjacency and topology overlap function options
      power = power, networkType = "signed", TOMType = "signed",
      # load previous TOMs
      saveTOMs = FALSE,
      # Basic tree cut options
      deepSplit = 2,  #default, known to reasonable
      minModuleSize = 300, #default 20, use 30 for transcriptome, or ncol(subDat)/2
      pamStage = TRUE, pamRespectsDendro = TRUE, #default, known to reasonable
      # Thredhold to merge modules: a height cut of 0.25 corresponding to correlation of 0.75
      mergeCutHeight = 0.25,
      # others
      reassignThreshold = 0,
      numericLabels = TRUE,
      verbose = 3)
    assign(paste0("cgnP",power,substring(corM,1,1)),cgn)
    save(list=grep("cgnP",ls(),value=TRUE), file = "wgcna.salt.Rdata")
  }
}
#Plotting
rm(list=ls())
library(RColorBrewer)
library(flashClust)
library(WGCNA)
options(stringsAsFactors = FALSE)
allowWGCNAThreads()
#Loading intermediate data
load("wgcna.salt.prep.Rdata")
load("wgcna.salt.Rdata")
ls()
i="cgnP24p"
net = get(i)
datExpr=as.data.frame(t(expr))
dim(datExpr)
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
table(moduleColors)
table(net$colors)


sizeGrWindow(12, 9)
mergedColors = labels2colors(net$colors)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

MEs1 = net$MEs;
geneTree = net$dendrograms[[1]];
MEs = moduleEigengenes(datExpr, moduleColors)$eigengenes
MET = orderMEs(MEs)
dim(MET)
save(MEs,MEs1 ,moduleLabels, moduleColors, geneTree,
           file = "02-networkConstruction-auto.RData")

#
pdf(file = "Eigengene adjacency heatmap.pdf", width = 12, height = 9)
plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2), plotDendrograms = FALSE, xLabelsAngle = 90)
dev.off()

#Importing phenotypic data
traitData=read.table('phenotype_data.txt',header=T,sep="\t")
allTraits = traitData
dim(allTraits)
#[1] 72  6
names(allTraits)
#[1] "Sample" "MDA"    "H2O2"   "Chla"   "Chlb"   "Chla.b"
datExpr=t(expr)
Samples = rownames(datExpr)
traitRows = match(Samples, allTraits$Sample);
datTraits = allTraits[traitRows, -1];
rownames(datTraits) = allTraits[traitRows, 1];
collectGarbage()
sampleTree2 = hclust(dist(datExpr), method = "average")
traitColors = numbers2colors(datTraits, signed = FALSE) #Use color to represent relevance
plotDendroAndColors(sampleTree2, traitColors,
                                         groupLabels = names(datTraits),
                                          main = "Sample dendrogram and trait heatmap")
save(datExpr, datTraits, file = "01-dataInput.RData")

load(file = "01-dataInput.RData");
load(file = "02-networkConstruction-auto.RData");
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
sizeGrWindow(10,6)
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation coefficient as a heat map
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 1.0,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
#Draw a picture for each module
pdf(file = "Greenyellowexpressionheatmap.pdf", width = 12, height = 9)
which.module="greenyellow"
ME=MEs[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(datExpr[,moduleColors==which.module ]) ),
        nrgcols=30,rlabels=F,rcols=which.module,
        main=which.module, cex.main=2)
par(mar=c(5, 4.2, 0, 0.7))
barplot(ME, col=which.module, main="", cex.main=2,
        ylab="eigengene expression",xlab="hub gene")
dev.off()


# WGCNA hub gene heatmap
MEs1=t(MEs)
MEs2=MEs1
c=data.frame(cbind(apply(MEs2[,1:3],1,mean),
                   apply(MEs2[,4:6],1,mean),
                   apply(MEs2[,7:9],1,mean),
                   apply(MEs2[,10:12],1,mean),
                   apply(MEs2[,13:15],1,mean),
                   apply(MEs2[,16:18],1,mean),
                   apply(MEs2[,19:21],1,mean),
                   apply(MEs2[,22:24],1,mean),
                   apply(MEs2[,25:27],1,mean),
                   apply(MEs2[,28:30],1,mean),
                   apply(MEs2[,31:33],1,mean),
                   apply(MEs2[,34:36],1,mean),
                   apply(MEs2[,37:39],1,mean),
                   apply(MEs2[,40:42],1,mean),
                   apply(MEs2[,43:45],1,mean),
                   apply(MEs2[,46:48],1,mean),
                   apply(MEs2[,49:51],1,mean),
                   apply(MEs2[,52:54],1,mean),
                   apply(MEs2[,55:57],1,mean),apply(MEs2[,58:60],1,mean),apply(MEs2[,61:63],1,mean),apply(MEs2[,64:66],1,mean),apply(MEs2[,67:69],1,mean),apply(MEs2[,70:72],1,mean) ))


class(c)
names(c) <- c("a", "b", "c","d","e","f","g","h","i","j",
              "k","l","m","n","o","p","q","r","s","t","u","v","w","x")


c$module=names(MEs)
#c$module=c("ME7(1624)","ME3(5241)","ME4(3946)","ME10(798)","ME2(5893)","ME9(1414)","ME5(3289)","ME6(2220)","ME1(6156)","ME11(422)","ME8(1590)","ME0(5409)" )
c$module=c("black(1624)","brown(5241)","yellow(3946)","purple(798)","blue(5893)","magenta(1414)","green(3289)","red(2220)","turquoise(6156)","greenyellow(422)","pink(1590)","grey(5409)" )
c <- c[, c("module","a", "b", "c","d","e","f","g","h","i","j","k","l","m","n","o","p","q","r","s","t","u","v","w","x")]
anno=read.table("annot0.txt",header=T,sep="\t")
Period=anno[,-1]
head(Period)

class(Period)
dim(Period)



per=c
dat2 <- per[,-1]*10
rownames(Period) = rownames(t(dat2))
row.names(dat2)=per$module
library(pheatmap)
Period$time<- factor(Period$time, levels=c("0.5h", "1h","3h","6h","9h","12h","15h","18h","24h","48h","72h","168h"), ordered=TRUE)
Period$treatment<- factor(Period$treatment, levels=c("CK","Salt"), ordered=TRUE)
ann_colors1 = list(time = c("0.5h"="#fffff5", "1h"= "#fab1ce","3h"="#EE7785","6h"="#ef5285","9h"="#C5E99B","12h"="#8CD790","15h"="#77AF9C","18h"="#379392","24h"="#285943","48h"="#0000FF","72h"="#8A2BE2","168h"="#4B0082"),
                   treatment = c("CK"="#47b8e0","Salt"="#dedcee"))

pheatmap(dat2[c(-12),],annotation_col = Period,
                   annotation_colors = ann_colors1,
                   scale="none",cluster_cols=F,cluster_rows =T,display_numbers = F,
                   fontsize = 15,show_colnames = F,
                   gaps_col = c(12))

table(moduleColors)
#moduleColors
#black        blue       brown       green greenyellow        grey     magenta        pink      purple         red
#1624        5893        5241        3289         422        5409        1414        1590         798        2220
#turquoise      yellow
#6156        3946
table(net$colors)
#0    1    2    3    4    5    6    7    8    9   10   11
#5409 6156 5893 5241 3946 3289 2220 1624 1590 1414  798  422
#Corresponding module genes
moduleGenes = names(net$colors)[net$colors%in%c("11")]
write.csv(moduleGenes,"greenyellow11.csv")

#Batch obtain the correlation between each module gene and each trait GS, gene and module feature vector MM
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))

names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep = "")

#names of those trait
traitNames = names(datTraits)

geneTraitSignificance = as.data.frame(cor(datExpr, datTraits, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))

names(geneTraitSignificance) = paste("GS.", traitNames, sep="")
names(GSPvalue) = paste("p.GS.", traitNames, sep = "")


dir.create("9_trait_module_Module membership vs gene significance", showWarnings = FALSE)
for (trait in traitNames){
  traitColumn=match(trait,traitNames)
  for (module in modNames){
    column = match(module, modNames)
    moduleGenes = moduleColors == module
    if (nrow(geneModuleMembership[moduleGenes,]) > 1){
      pdf(file = paste("9_trait_module_Module membership vs gene significance/9_", trait, "_", module,"_Module membership vs gene significance.pdf", sep = ""), width = 7, height = 7)
      par(mfrow = c(1,1))
      verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                         abs(geneTraitSignificance[moduleGenes, traitColumn]),
                         xlab = paste("Module Membership in", module, "module"),
                         ylab = paste("Gene significance for ",trait),
                         main = paste("Module membership vs. gene significance\n"),
                         cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
      dev.off()
    }
  }
}


names(EXPR)
probes = names(EXPR)

#################export GS and MM###############
geneInfo0 = data.frame(probes = probes,
                       moduleColor = moduleColors)

for (Tra in 1:ncol(geneTraitSignificance))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneTraitSignificance[,Tra],
                         GSPvalue[, Tra])
  names(geneInfo0) = c(oldNames,names(geneTraitSignificance)[Tra],
                       names(GSPvalue)[Tra])
}

for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[,mod],
                         MMPvalue[, mod])
  names(geneInfo0) = c(oldNames,names(geneModuleMembership)[mod],
                       names(MMPvalue)[mod])
}
geneOrder =order(geneInfo0$moduleColor)
geneInfo = geneInfo0[geneOrder, ]
write.table(geneInfo, file = "10_GS_and_MM.xls", sep = "\t", row.names = F)


# Recalculate the eigengenes of the module
MEs = moduleEigengenes(datExpr, moduleColors)$eigengenes
# Extracting phenotypic data
MDA = as.data.frame(datTraits$MDA)
H2O2=as.data.frame(datTraits$H2O2)
Chla=as.data.frame(datTraits$Chla)
Chlb=as.data.frame(datTraits$Chlb)
Chla.b=as.data.frame(datTraits$Chla.b)
names(MDA) = "MDA"
names(H2O2) = "H2O2"
names(Chla) = "Chla"
names(Chlb) = "Chlb"
names(Chla.b) = "Chla+b"
# Add to the corresponding module
MET = orderMEs(cbind(MEs, MDA,H2O2,Chla,Chlb,Chla.b))
#Plotting
sizeGrWindow(6,6);
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene dendrogram", marDendro = c(0,4,2,0),
                      plotHeatmaps = FALSE)
sizeGrWindow(6,6);
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2),
                      plotDendrograms = FALSE, xLabelsAngle = 90)
