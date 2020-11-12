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

#im going to import the dataframes with OTUs of low confidence

otu_data_table <- read_csv(file = "CVG_otutable_updated.csv") #import otu table

tax_data_table <- read_tsv(file = "CVG_taxatable_updated.tsv") #import taxa table

otu_data_matrix <- data.matrix(otu_data_table) #make the otu table into a matrix (it has to be a matrix)
rownames(otu_data_matrix) <- otu_data_table$OTU #relabel the column names at the OTU numbers
otu_data_matrix <- otu_data_matrix[,-1] #remove null column

CVG_otu_table <- otu_table(otu_data_matrix, taxa_are_rows = T) #turn OTU table into a phyloseq object

tax_matrix <- as.matrix(tax_data_table) #make the taxa table into a matrix (it has to be a matrix)
rownames(tax_matrix) <- tax_data_table$OTU #relabel the column names to be the OTU numbers
tax_matrix <- tax_matrix[,-1]  #remove null column

CVG_tax_table <- tax_table(tax_matrix)



#dna_CVG_sample_table is the dataframe with ng of dna included

CVG_phyloseq <- merge_phyloseq(CVG_otu_table, CVG_tax_table, dna_CVG_sample_table, tree.fasttree) #merge phyloseq objects
########################################################
#above is the new phyloseq object incorportating the trimmed taxa

########################################################
#now i'm going to do to data cleaning checks

# Create table, number of features for each phyla
table(tax_table(CVG_phyloseq)[, "Phylum"], exclude = NULL)

# Compute prevalence of each feature, store as data.frame
prevdf = apply(X = otu_table(CVG_phyloseq),
               MARGIN = ifelse(taxa_are_rows(CVG_phyloseq), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(CVG_phyloseq),
                    tax_table(CVG_phyloseq))

plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})