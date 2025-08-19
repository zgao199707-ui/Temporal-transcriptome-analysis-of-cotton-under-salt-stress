getwd()
setwd("C:/Users/17294/Desktop/RNA_analysis/MetaCycle/")
rm(list=ls())
library("MetaCycle")

# Run MetaCycle analysis for CK and Salt datasets
meta2d(infile="38002CK_TPM_AVER(0-18_3H).csv", filestyle="csv",
       outdir="CK_new/", analysisStrategy = "auto",
       timepoints="Line1", outRawData=TRUE)

meta2d(infile="38002Salt_TPM_AVER(0-18_3H).csv", filestyle="csv",
       outdir="Salt_new/", analysisStrategy = "auto",
       timepoints="Line1", outRawData=TRUE)


### Plotting
rm(list=ls())
getwd()
setwd("C:/Users/17294/Desktop/Time2025/Metacycle/")

# Install and load required packages
library(ggplot2)
library(ggsignif)
library(rstatix)
library(tidyr)

# Read data
data <- read.table("three.rAMP.box.txt", sep = "\t", header = TRUE)
data <- log(data+1)

# Preview the data
head(data)

# Convert data from wide format to long format
long_data <- pivot_longer(data, 
                          cols = everything(), 
                          names_to = "group", 
                          values_to = "value")

# Preview the transformed data
head(long_data)

# Remove rows containing NA values  
long_data <- na.omit(long_data)  

# Perform pairwise t-tests and adjust p-values
pairwise_result <- long_data %>%
  t_test(value ~ group) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()

# Show t-test results
print(pairwise_result)


## Define comparison groups  
compare_groups <- list(  
  c("conserved.C", "conserved.S"),  
  c("gain.C", "gain.S"),  
  c("loss.C", "loss.S")  
)  

# Create boxplots with significance annotations  
ggplot(long_data, aes(x = group, y = value, fill = group)) +  
  geom_boxplot() +  
  geom_signif(comparisons = compare_groups,  # use defined comparison groups  
              map_signif_level = TRUE,  
              textsize = 6) +  
  theme_minimal() +  
  labs(title = "Pairwise Comparisons of rAMP",  
       y = "rAMP",   
       x = "") +  
  theme(  
    legend.position = "none",              
    axis.text.x = element_text(size = 18),  
    axis.title.x = element_text(size = 14),   
    axis.text.y = element_text(size = 14),  
    axis.title.y = element_text(size = 16),   
    plot.title = element_text(size = 18)    
  ) 

# Export options: PNG 1200x900 / PDF 12x9
