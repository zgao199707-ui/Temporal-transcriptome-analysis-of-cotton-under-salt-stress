rm(list=ls())
countData <-read.table("rawdata/kallisto_estCount.txt", head=TRUE)
colData <-read.table("rawdata/sample.txt", head=TRUE)
all(colnames(countData)==rownames(colData))
colData$time <- factor(colData$time, levels = c("0.5","1","3","6","9","12","15","18","24","48","72","168"))
View(colData)
colData$condition <- factor(colData$condition, levels = c("Control","Salt"))

library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = colData,
                              design = ~ condition+time+condition:time)

# Filter genes with very low counts
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

dds = DESeq(dds)
resultsNames(dds)

mfResList = list()
library(tidyverse)

# Generate contrasts for each time point (Salt vs Control)
mfResList["SaltvsControl_time0.5h"] = results(dds,contrast = list(c("condition_Salt_vs_Control"))) %>%list()
mfResList["SaltvsControl_time1h"] = results(dds,contrast = list(c("condition_Salt_vs_Control","conditionSalt.time1"))) %>%list()
mfResList["SaltvsControl_time3h"] = results(dds,contrast = list(c("condition_Salt_vs_Control","conditionSalt.time3"))) %>%list()
mfResList["SaltvsControl_time6h"] = results(dds,contrast = list(c("condition_Salt_vs_Control","conditionSalt.time6"))) %>%list()
mfResList["SaltvsControl_time9h"] = results(dds,contrast = list(c("condition_Salt_vs_Control","conditionSalt.time9"))) %>%list()
mfResList["SaltvsControl_time12h"] = results(dds,contrast = list(c("condition_Salt_vs_Control","conditionSalt.time12"))) %>%list()
mfResList["SaltvsControl_time15h"] = results(dds,contrast = list(c("condition_Salt_vs_Control","conditionSalt.time15"))) %>%list()
mfResList["SaltvsControl_time18h"] = results(dds,contrast = list(c("condition_Salt_vs_Control","conditionSalt.time18"))) %>%list()
mfResList["SaltvsControl_time24h"] = results(dds,contrast = list(c("condition_Salt_vs_Control","conditionSalt.time24"))) %>%list()
mfResList["SaltvsControl_time48h"] = results(dds,contrast = list(c("condition_Salt_vs_Control","conditionSalt.time48"))) %>%list()
mfResList["SaltvsControl_time72h"] = results(dds,contrast = list(c("condition_Salt_vs_Control","conditionSalt.time72"))) %>%list()
mfResList["SaltvsControl_time168h"] = results(dds,contrast = list(c("condition_Salt_vs_Control","conditionSalt.time168"))) %>%list()

setwd("C:/Users/17294/Desktop/Time_RNAseq_analysis/Multi_factor_DEG2024.7.10/")
savePath = "C:/Users/17294/Desktop/Time_RNAseq_analysis/Multi_factor_DEG2024.7.10/"

# Save results for each time point
for (myname in names(mfResList)) {
  element = mfResList[[myname]]
  print(summary(element,alpha=.05)) # Print summary results
  write(paste0("#",capture.output(summary(element,alpha=.05))),
        paste(savePath,myname,".txt", sep=""))
  write.table(as.data.frame(element),
              file=paste(savePath,myname,".txt", sep=""),
              sep="\t",append = TRUE)
}

# Custom function: summarize DEGs across pairwise comparisons
myGetDiffSummary = function(dirpath,outcsvfile){
  pairwiselist=list.files(dirpath)
  pwdatalist = list()
  for(pwname in pairwiselist){
    pwdata = read.delim(paste0(dirpath,"/",pwname),comment.char = "#")
    pwdatalist[[pwname]] = pwdata
  }
  pwdsummary=data.frame()
  for(pwname in pairwiselist){
    genedfdata = pwdatalist[[pwname]]%>%filter(!is.na(padj),padj<=0.05)
    upgene = filter(genedfdata,log2FoldChange>0)
    downgene = filter(genedfdata,log2FoldChange<0)
    numbersum = data.frame(
      At.up   = filter(upgene,grepl("Gohir.A",rownames(upgene)))%>%nrow,
      At.down = filter(downgene,grepl("Gohir.A",rownames(downgene)))%>%nrow,
      At.total= filter(genedfdata,grepl("Gohir.A",rownames(genedfdata)))%>%nrow,
      Dt.up   = filter(upgene,grepl("Gohir.D",rownames(upgene)))%>%nrow,
      Dt.down = filter(downgene,grepl("Gohir.D",rownames(downgene)))%>%nrow,
      Dt.total= filter(genedfdata,grepl("Gohir.D",rownames(genedfdata)))%>%nrow,
      Up      = nrow(upgene),
      down    = nrow(downgene),
      total   = nrow(genedfdata),
      row.names = pwname)
    pwdsummary = bind_rows(pwdsummary,numbersum)
  }
  write.csv(pwdsummary,outcsvfile)
}

# Move pairwise comparison results into Multi_factor_pair_DEG folder and generate summary CSV
myGetDiffSummary("C:/Users/17294/Desktop/Time_RNAseq_analysis/Multi_factor_DEG2024.7.10/Multi_factor_pair_DEG/",
                 "C:/Users/17294/Desktop/Time_RNAseq_analysis/Multi_factor_DEG2024.7.10//Multifactor_diffsummary.csv")
