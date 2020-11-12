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

#in order to further interrogate the contamination of the data, I will collapse the OTUs into their genera
plot_bar(dna_CVG_phyloseq, fill = "Phylum")

#the dataframe in use is dna_CVG_phyloseq
#the model to collapse by genus is testglom052 = tax_glom(CVG052, taxrank = "Genus")
dna_CVG_phyloseq <- merge_phyloseq(CVG_otu_table, CVG_tax_table, dna_CVG_sample_table) #merge phyloseq objects

dna_CVG_phyloseq_genus = tax_glom(dna_CVG_phyloseq, taxrank = "Genus")

#qucik access to taxonomic identifications
TEMP = tax_table(dna_CVG_phyloseq_genus)
head(TEMP)
TEMP["OTU_5"]

TEMP = otu_table(dna_CVG_phyloseq_genus)
head(TEMP)
TEMP["OTU_19"]
#this collapses 'taxa' from 599 to 210
#when dodgy classifications are removed, it's 593 to 195.
#i'm assuming i'm confident about my taxonoic classifications
#double check how the tax_glom works

########################################################

############################################################################

#first we can inspect the library size: see if theres a pattern between controls and samples

df <- as.data.frame(sample_data(dna_CVG_phyloseq_genus)) # Put sample_data into a ggplot-friendly data.frame
df$LibrarySize <- sample_sums(dna_CVG_phyloseq_genus)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
ggplot(data=df, aes(x=Index, y=LibrarySize, color=Sample_or_Control)) + geom_point()

############################################################################
#using the frequency meothd first
genuscontamdf.freq <- isContaminant(dna_CVG_phyloseq_genus, method="frequency", conc="DNA_conc")
head(genuscontamdf.freq)
hist(genuscontamdf.freq$p)

#using the prevalence method
sample_data(dna_CVG_phyloseq_genus)$is.neg <- sample_data(dna_CVG_phyloseq_genus)$Sample_or_Control == "Control"
contamdf.prev <- isContaminant(dna_CVG_phyloseq_genus, method="prevalence", neg="is.neg")
table(contamdf.prev$contaminant)

#using the combined method (with Fishers)
genuscontamdf.comb03 <- isContaminant(dna_CVG_phyloseq_genus, method="combined", conc="DNA_conc", threshold=0.3, neg = "is.neg", detailed = TRUE)
table(genuscontamdf.comb03$contaminant)
hist(genuscontamdf.comb03$p)

write.csv(contamdf.comb03, file = "genuscontamdf.comb03.csv")

############################################################################


#i'm going to look at the histogram of P values

hist(genuscontamdf.comb03$p)
ggsave("contamdfcombpspreas.pdf", width = 37, height = 25, units = c('cm'))

hist(genuscontamdf.freq$p)
hist(genuscontamdf.prev$p)
hist(contamdf.prev05$p)

ggplot(data=genuscontamdf.freq, aes(genuscontamdf.freq$p)) + 
  geom_histogram(breaks=seq(0, 1, by = 0.1), 
                 col="#208A5B", 
                 fill="#9C5ADA", 
                 alpha = .8) + 
  labs(title="Frequency-based Decontam: Genus") +
  labs(x="P value", y="Frequency") + theme(
    plot.title = element_text(color="black", size=28, face="bold"),
    axis.title.x = element_text(color="black", size=28, face="bold"),
    axis.title.y = element_text(color="black", size=28, face="bold"),
    axis.text=element_text(size=26),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 20))
ggsave("decontam_frequency_genus.pdf", width = 25, height = 25, units = c('cm'))

ggplot(data=genuscontamdf.prev, aes(genuscontamdf.prev$p)) + 
  geom_histogram(breaks=seq(0, 1, by = 0.1), 
                 col="#208A5B", 
                 fill="#9C5ADA", 
                 alpha = .8) + 
  labs(title="Prevalence-based Decontam: Genus") +
  labs(x="P value", y="Frequency") + theme(
    plot.title = element_text(color="black", size=28, face="bold"),
    axis.title.x = element_text(color="black", size=28, face="bold"),
    axis.title.y = element_text(color="black", size=28, face="bold"),
    axis.text=element_text(size=26),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 20))
ggsave("decontam_prevalence_genus.pdf", width = 25, height = 25, units = c('cm'))

ggplot(data=genuscontamdf.comb03, aes(genuscontamdf.comb03$p)) + 
  geom_histogram(breaks=seq(0, 1, by = 0.1), 
                 col="#208A5B", 
                 fill="#9C5ADA", 
                 alpha = .8) + 
  labs(title="Combined Decontam: Genus") +
  labs(x="P value", y="Frequency") + theme(
    plot.title = element_text(color="black", size=28, face="bold"),
    axis.title.x = element_text(color="black", size=28, face="bold"),
    axis.title.y = element_text(color="black", size=28, face="bold"),
    axis.text=element_text(size=26),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 20))
ggsave("decontam_combined_genus.pdf", width = 25, height = 25, units = c('cm'))
##############################################################

#using the combined method (with Fishers)
genuscontamdf.comb03 <- isContaminant(dna_CVG_phyloseq_genus, method="combined", conc="DNA_conc", threshold=0.3, neg = "is.neg", detailed = TRUE)
table(genuscontamdf.comb03$contaminant)

write.csv(genuscontamdf.comb03, file = "genuscontamdf.comb03.csv")

#using the frequency meothd first
genuscontamdf.freq02 <- isContaminant(dna_CVG_phyloseq_genus, method="frequency", threshold=0.2, conc="DNA_conc")
table(genuscontamdf.freq02$contaminant)
write.csv(genuscontamdf.freq02, file = "genuscontamdf.freq02.csv")

#using the prevalence method
genuscontamdf.prev02 <- isContaminant(dna_CVG_phyloseq_genus, method="prevalence", threshold=0.2, neg="is.neg")
table(genuscontamdf.prev02$contaminant)
write.csv(genuscontamdf.prev02, file = "genuscontamdf.prev02.csv")

##############################################################
#lets have a look
plot_frequency(dna_CVG_phyloseq_genus, taxa_names(dna_CVG_phyloseq_genus)[c(3, 5, 70)], conc="DNA_conc") + 
  xlab("DNA Concentration (Tapestation or Qubit)") +theme_bw()
ggsave("decontam_frdfgdfgeq_6.pdf", width = 37, height = 25, units = c('cm'))