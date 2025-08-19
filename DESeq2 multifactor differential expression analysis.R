### Clear environment
rm(list=ls())
countData <-read.table("kallisto_estCount.txt", head=TRUE)
colData <-read.table("sample.txt", head=TRUE)
all(colnames(countData)==rownames(colData))
colData$time <- factor(colData$time, levels = c("0.5","1","3","6","9","12","15","18","24","48","72","168"))
View(colData)
colData$condition <- factor(colData$condition, levels = c("Control","Salt"))
library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData = countData,colData = colData,design = ~ condition+time+condition:time)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds.reduced.time <- DESeq(dds, test="LRT", reduced=~condition)
res.dds.reduced.time <- results(dds.reduced.time)
summary(res.dds.reduced.time)
dds.reduced.condition <- DESeq(dds, test="LRT", reduced=~time)
res.dds.reduced.condition<-results(dds.reduced.condition)
summary(res.dds.reduced.condition)
dds.reduced.ct <- DESeq(dds, test="LRT", reduced=~condition+time)
res.dds.reduced.ct<-results(dds.reduced.ct)
summary(res.dds.reduced.ct)
write.csv(res.dds.reduced.condition, file = "condition.LRT.csv")
write.csv(res.dds.reduced.ct, file = "condition_time.LRT.csv")
write.csv(res.dds.reduced.time, file = "time.LRT.csv")

####
# LRT analysis
rm(list=ls())
getwd()
library(DESeq2)
setwd("C:/Users/17294/Desktop/HEB_analysis/new_HEB/")
countData <-read.table("A+D.count.txt", head=TRUE)
colData <-read.table("sample_contiuous.time.txt", head=TRUE)
all(colnames(countData)==rownames(colData))
colData$time <- factor(colData$time, levels = c("0.5","1","3","6","9","12","15","18","24","48","72","168"))
colData$condition <- factor(colData$condition, levels = c("Control","Salt"))
colData$homoeolog <- factor(colData$homoeolog, levels = c("A","D"))
dds_full <- DESeqDataSetFromMatrix(countData = countData,colData = colData,design = ~time + homoeolog + condition + homoeolog:time + homoeolog:condition + time:condition)
keep <- rowSums(counts(dds_full)) >= 10
dds_full <- dds_full[keep, ]

# LRT: homoeolog × condition
dds.LRT.homoeolog_condition <- DESeq(dds_full, test="LRT", reduced=~time + homoeolog + condition + homoeolog:time + time:condition)
res.dds.LRT.homoeolog_condition<-results(dds.LRT.homoeolog_condition)
summary(res.dds.LRT.homoeolog_condition)
write.csv(res.dds.LRT.homoeolog_condition, file = "res.dds.LRT.homoeolog_condition.csv")

# LRT: homoeolog × time
dds.LRT.homoeolog_time <- DESeq(dds_full, test="LRT", reduced=~time + homoeolog + condition + homoeolog:condition + time:condition)
res.dds.LRT.homoeolog_time<-results(dds.LRT.homoeolog_time)
summary(res.dds.LRT.homoeolog_time)
write.csv(res.dds.LRT.homoeolog_time, file = "res.dds.LRT.homoeolog_time.csv")

# LRT: homoeolog × condition × time
dds_full2 <- DESeqDataSetFromMatrix(countData = countData,colData = colData,design = ~time + homoeolog + condition + homoeolog:time + homoeolog:condition + time:condition + homoeolog:time:condition)
keep <- rowSums(counts(dds_full2)) >= 10
dds_full2 <- dds_full2[keep, ]
dds.LRT.homoeolog_condition_time <- DESeq(dds_full2, test="LRT", reduced=~time + homoeolog + condition + homoeolog:time + homoeolog:condition + time:condition)
res.dds.LRT.homoeolog_condition_time<-results(dds.LRT.homoeolog_condition_time)
summary(res.dds.LRT.homoeolog_condition_time)
write.csv(res.dds.LRT.homoeolog_condition_time, file = "res.dds.LRT.homoeolog_condition_time.csv")

# Homoeolog effect
rm(list=ls())
getwd()
library(DESeq2)
setwd("C:/Users/17294/Desktop/HEB_analysis/new_HEB/")
countData <-read.table("A+D.count.txt", head=TRUE)
colData <-read.table("sample_contiuous.time.txt", head=TRUE)
all(colnames(countData)==rownames(colData))
colData$time <- factor(colData$time, levels = c("0.5","1","3","6","9","12","15","18","24","48","72","168"))
colData$condition <- factor(colData$condition, levels = c("Control","Salt"))
colData$homoeolog <- factor(colData$homoeolog, levels = c("A","D"))
dds_full <- DESeqDataSetFromMatrix(countData = countData,colData = colData,design = ~time + homoeolog + condition + homoeolog:time + homoeolog:condition + time:condition)
keep <- rowSums(counts(dds_full)) >= 10
dds_full <- dds_full[keep, ]
dds<-DESeq(dds_full)
resultsNames(dds)
res_12_homoeolog_full<- results(dds, name="homoeolog_D_vs_A")
summary(res_12_homoeolog_full)
write.csv(res_12_homoeolog_full, file = "res_12_homoeolog_full.csv")


# LRT using 7 selected time points
rm(list=ls())
getwd()
library(DESeq2)
setwd("C:/Users/17294/Desktop/HEB_analysis/new_HEB/LRT_7TIME/")
countData <-read.table("A+D.count.txt", head=TRUE)
colData <-read.table("sample_contiuous.time.txt", head=TRUE)
all(colnames(countData)==rownames(colData))
colData$time <- factor(colData$time, levels = c("0.5","3","6","9","12","15","18"))
colData$condition <- factor(colData$condition, levels = c("Control","Salt"))
colData$homoeolog <- factor(colData$homoeolog, levels = c("A","D"))
dds_full <- DESeqDataSetFromMatrix(countData = countData,colData = colData,design = ~time + homoeolog + condition + homoeolog:time + homoeolog:condition + time:condition)
keep <- rowSums(counts(dds_full)) >= 10
dds_full <- dds_full[keep, ]

# LRT: homoeolog × condition
dds.LRT.homoeolog_condition <- DESeq(dds_full, test="LRT", reduced=~time + homoeolog + condition + homoeolog:time + time:condition)
res.dds.LRT.homoeolog_condition<-results(dds.LRT.homoeolog_condition)
summary(res.dds.LRT.homoeolog_condition)
write.csv(res.dds.LRT.homoeolog_condition, file = "res.dds.LRT.homoeolog_condition.csv")

# LRT: homoeolog × time
dds.LRT.homoeolog_time <- DESeq(dds_full, test="LRT", reduced=~time + homoeolog + condition + homoeolog:condition + time:condition)
res.dds.LRT.homoeolog_time<-results(dds.LRT.homoeolog_time)
summary(res.dds.LRT.homoeolog_time)
write.csv(res.dds.LRT.homoeolog_time, file = "res.dds.LRT.homoeolog_time.csv")

# LRT: homoeolog × condition × time
dds_full2 <- DESeqDataSetFromMatrix(countData = countData,colData = colData,design = ~time + homoeolog + condition + homoeolog:time + homoeolog:condition + time:condition + homoeolog:time:condition)
keep <- rowSums(counts(dds_full2)) >= 10
dds_full2 <- dds_full2[keep, ]
dds.LRT.homoeolog_condition_time <- DESeq(dds_full2, test="LRT", reduced=~time + homoeolog + condition + homoeolog:time + homoeolog:condition + time:condition)
res.dds.LRT.homoeolog_condition_time<-results(dds.LRT.homoeolog_condition_time)
summary(res.dds.LRT.homoeolog_condition_time)
write.csv(res.dds.LRT.homoeolog_condition_time, file = "res.dds.LRT.homoeolog_condition_time.csv")

# Homoeolog effect (7 time points)
rm(list=ls())
getwd()
library(DESeq2)
setwd("C:/Users/17294/Desktop/HEB_analysis/new_HEB/LRT_7TIME/")
countData <-read.table("A+D.count.txt", head=TRUE)
colData <-read.table("sample_contiuous.time.txt", head=TRUE)
all(colnames(countData)==rownames(colData))
colData$time <- factor(colData$time, levels = c("0.5","3","6","9","12","15","18"))
colData$condition <- factor(colData$condition, levels = c("Control","Salt"))
colData$homoeolog <- factor(colData$homoeolog, levels = c("A","D"))
dds_full <- DESeqDataSetFromMatrix(countData = countData,colData = colData,design = ~time + homoeolog + condition + homoeolog:time + homoeolog:condition + time:condition)
keep <- rowSums(counts(dds_full)) >= 10
dds_full <- dds_full[keep, ]
dds<-DESeq(dds_full)
resultsNames(dds)
res_7_homoeolog_full<- results(dds, name="homoeolog_D_vs_A")
summary(res_7_homoeolog_full)
write.csv(res_7_homoeolog_full, file = "res_7_homoeolog_full.csv")
