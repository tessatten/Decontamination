if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("decontam")

library("ape")
library("ggplot2")
library("plyr")
library("scales")
library("vegan")
library("knitr")
library("dplyr")
library("praise")
library("tidyverse")
library("phyloseq")
library("decontam")

#im going to add the ng info from Ronan's spreadsheet

ngData <- read.csv(file = "after_pcr_1.csv", stringsAsFactors = F) #import otu table
ngData <- subset(ngData, select = c(SampleID, PCR1_qubit_ngul, Sample_or_control))

colnames(ngData)[colnames(ngData)=="SampleID"] <- "Sample_ID"
ngData$Sample_ID <- paste("CVG", ngData$Sample_ID, sep="")

testngData <- merge(x = sample_data_table, y = ngData, by = "Sample_ID", all = TRUE)
write.csv(testngData, file = "testngData.csv")



#import into phyloseq object

ng_sample_matrix <- as.matrix(testngData) #make it a matrix

rownames(ng_sample_matrix) <- testngData$Sample_ID #change column names, it adds an X in front of the number
ng_sample_matrix <- ng_sample_matrix[,-1]  #remove null column
ng_sample_matrix <- as.data.frame(ng_sample_matrix) #make it a data frame again
#sam_data is depreciated so switched to sample_data
#required a data frame so switched it back to a data frame


ng_CVG_sample_table <- sample_data(ng_sample_matrix) #import to phyloseq
#it does import!!!
#but does not import ntc etc


ng_CVG_phyloseq <- merge_phyloseq(CVG_otu_table, CVG_tax_table, ng_CVG_sample_table) #merge phyloseq objects
########################################################

# so ng_CVG_phyloseq is the phyloseq object with the new info for decontam added

head(sample_data(ng_CVG_phyloseq))

df <- as.data.frame(sample_data(ng_CVG_phyloseq)) # Put sample_data into a ggplot-friendly data.frame
df$LibrarySize <- sample_sums(ng_CVG_phyloseq)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
ggplot(data=df, aes(x=Index, y=LibrarySize, color=Sample_or_Control)) + geom_point()