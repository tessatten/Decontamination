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

otu_data_table <- read_csv(file = "CVG_otutable.csv") #import otu table
#otu_data_table <- read_csv(file = "CVG_otutable_foranosim.csv") #import otu table
otu_data_matrix <- data.matrix(otu_data_table) #make the otu table into a matrix (it has to be a matrix)
rownames(otu_data_matrix) <- otu_data_table$OTU #relabel the column names at the OTU numbers
otu_data_matrix <- otu_data_matrix[,-1] #remove null column

contam_CVG_otu_table <- otu_table(otu_data_matrix, taxa_are_rows = T) #turn OTU table into a phyloseq object


#first i will import the dna conc info from Ronan 

dnaData <- read.csv(file = "DNA_conc_all_samples.csv", stringsAsFactors = F) #import otu table
dnaData$Sample_ID <- paste("CVG", dnaData$Sample_ID, sep="")

#now i'll add these columns to the sample metadata file

testdnaData <- merge(x = sample_data_table, y = dnaData, by = "Sample_ID", all = TRUE)
write.csv(testdnaData, file = "testdnaData.csv")

#reimport checked metadata file
dnaData <- read.csv(file = "testdnaData.csv", stringsAsFactors = F) #import otu table
dnaData <- dnaData[,-1]  #remove null column
#import into phyloseq object

dna_sample_matrix <- as.matrix(dnaData) #make it a matrix

rownames(dna_sample_matrix) <- dnaData$Sample_ID #change column names, it adds an X in front of the number
dna_sample_matrix <- dna_sample_matrix[,-1]  #remove null column
dna_sample_matrix <- as.data.frame(dna_sample_matrix) #make it a data frame again
#sam_data is depreciated so switched to sample_data
#required a data frame so switched it back to a data frame

# adding in cat col
dna_sample_matrix$bioMass <- if_else(dna_sample_matrix$DNA_conc < 0.0014, "Low",
                                     if_else(dna_sample_matrix$DNA_conc >= 0.0014 & dna_sample_matrix$DNA_conc <= 0.0576, "Mid", "High"))
dna_sample_matrix$bioMass <- as.factor(dna_sample_matrix$bioMass)
class(dna_sample_matrix$bioMass)
dna_sample_matrix$bioMass
#i need to check here if the DNA concs are numerical
class(dna_sample_matrix$DNA_conc)
dna_sample_matrix$DNA_conc <- as.numeric(levels(dna_sample_matrix$DNA_conc))[dna_sample_matrix$DNA_conc]
class(dna_sample_matrix$DNA_conc)
test <- read.csv(file = "dna_CVG_sample_table.csv", stringsAsFactors = F) #import otu table
write.csv(dna_CVG_sample_table, file = "dna_CVG_sample_table.csv")
dna_CVG_sample_table <- sample_data(dna_sample_matrix) #import to phyloseq
#it does import!!!
#but does not import ntc etc

tax_data_table <- read_tsv(file = "CVG_taxatable_updated.tsv") #import taxa table

tax_matrix <- as.matrix(tax_data_table) #make the taxa table into a matrix (it has to be a matrix)
rownames(tax_matrix) <- tax_data_table$OTU #relabel the column names to be the OTU numbers
tax_matrix <- tax_matrix[,-1]  #remove null column

contam_CVG_tax_table <- tax_table(tax_matrix)

sample_names(all_OTUS_CVG_phyloseq)
all_OTUS_CVG_phyloseq <- merge_phyloseq(contam_CVG_otu_table, contam_CVG_tax_table, dna_CVG_sample_table) #merge phyloseq objects
genus_dna_CVG_phyloseq <-tax_glom(all_OTUS_CVG_phyloseq, taxrank = "Genus")
########################################################

#so now i have a phyloseq object that include the ng of dna present in each sample after pcr 1

head(sample_data(dna_CVG_phyloseq))

########################################################
#first we can inspect the library size: see if theres a pattern between controls and samples

df <- as.data.frame(sample_data(dna_CVG_phyloseq)) # Put sample_data into a ggplot-friendly data.frame
df$LibrarySize <- sample_sums(dna_CVG_phyloseq)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
ggplot(data=df, aes(x=Index, y=LibrarySize, color=Sample_or_Control)) + geom_point()

#this can be export to demo that the library size of many samples are smaller than those of controls...
########################################################
#working through more slowly
#see

########################################################

#now we are going to use the frequency method to look at contaminants

contamdf.freq <- isContaminant(all_OTUS_CVG_phyloseq, method="frequency", conc="DNA_conc", detailed = TRUE)
head(contamdf.freq)
write.csv(contamdf.freq, file = "OTU_decontam_freq0.1.csv")

sample_data(all_OTUS_CVG_phyloseq)$is.neg <- sample_data(all_OTUS_CVG_phyloseq)$Sample_or_Control == "Control"
contamdf.prev <- isContaminant(all_OTUS_CVG_phyloseq, method="prevalence", conc="DNA_conc", detailed = TRUE, neg="is.neg")
head(contamdf.prev)
write.csv(contamdf.prev, file = "OTU_decontam_prev0.1later.csv")
#how many are 'contaminants'?
table(contamdf.freq$contaminant)
table(contamdf.prev$contaminant)

#for interests sake i am looking at 0.2 as well
contamdf.freq02 <- isContaminant(dna_CVG_phyloseq, method="frequency", conc="DNA_conc", threshold=0.2)
table(contamdf.freq02$contaminant)


#which ones?
head(which(contamdf.freq$contaminant))
which(contamdf.freq$contaminant)

#let's have a look at them 
plot_frequency(dna_CVG_phyloseq, taxa_names(dna_CVG_phyloseq)[c(8,10)], conc="DNA_conc") + 
  xlab("DNA Concentration (Tapestation or Qubit)")

set.seed(101)
plot_frequency(dna_CVG_phyloseq, taxa_names(dna_CVG_phyloseq)[sample(which(contamdf.freq$contaminant),3)], conc="DNA_conc") +
  xlab("DNA Concentration (Tapestation or Qubit)")
decontam_freq_1
#im going to plot the "contaminant" OTUS in sets of 3s
plot_frequency(dna_CVG_phyloseq, taxa_names(dna_CVG_phyloseq)[c(9,21,80)], conc="DNA_conc") + 
  xlab("DNA Concentration (Tapestation or Qubit)") +theme_bw() 
#+ ylim(min=-0.005, max=0.06)
ggsave("decontam_freq_1.pdf", width = 37, height = 25, units = c('cm'))

plot_frequency(dna_CVG_phyloseq, taxa_names(dna_CVG_phyloseq)[c(101,119,129)], conc="DNA_conc") + 
  xlab("DNA Concentration (Tapestation or Qubit)") +theme_bw()
ggsave("decontam_freq_2.pdf", width = 37, height = 25, units = c('cm'))

plot_frequency(dna_CVG_phyloseq, taxa_names(dna_CVG_phyloseq)[c(138,168,175)], conc="DNA_conc") + 
  xlab("DNA Concentration (Tapestation or Qubit)") +theme_bw()
ggsave("decontam_freq_3.pdf", width = 37, height = 25, units = c('cm'))

plot_frequency(dna_CVG_phyloseq, taxa_names(dna_CVG_phyloseq)[c(224,230,284)], conc="DNA_conc") + 
  xlab("DNA Concentration (Tapestation or Qubit)") +theme_bw()
ggsave("decontam_freq_4.pdf", width = 37, height = 25, units = c('cm'))

plot_frequency(dna_CVG_phyloseq, taxa_names(dna_CVG_phyloseq)[c(331,479,542)], conc="DNA_conc") + 
  xlab("DNA Concentration (Tapestation or Qubit)") +theme_bw()
ggsave("decontam_freq_5.pdf", width = 37, height = 25, units = c('cm'))

plot_frequency(dna_CVG_phyloseq, taxa_names(dna_CVG_phyloseq)[c(562,573)], conc="DNA_conc") + 
  xlab("DNA Concentration (Tapestation or Qubit)") +theme_bw()
ggsave("decontam_freq_6.pdf", width = 37, height = 25, units = c('cm'))


########################################################

plot_frequency(dna_CVG_phyloseq, taxa_names(dna_CVG_phyloseq)[c(3, 202)], conc="DNA_conc") + 
  xlab("DNA Concentration (Tapestation or Qubit)") +theme_bw()
ggsave("decontam_frdfgdfgeq_6.pdf", width = 37, height = 25, units = c('cm'))

########################################################

#now we are going to use the prevalence method to look at contaminants

#using 0.1 threshold
sample_data(dna_CVG_phyloseq)$is.neg <- sample_data(dna_CVG_phyloseq)$Sample_or_Control == "Control"
contamdf.prev <- isContaminant(dna_CVG_phyloseq, method="prevalence", neg="is.neg")
table(contamdf.prev$contaminant)
head(which(contamdf.prev$contaminant))
#using 0.5 threshold
contamdf.prev05 <- isContaminant(dna_CVG_phyloseq, method="prevalence", neg="is.neg", threshold=0.5)
table(contamdf.prev05$contaminant)
write.csv(contamdf.prev05, file = "contamdf.prev05.csv")

#using 0.2 threshold
contamdf.prev02 <- isContaminant(dna_CVG_phyloseq, method="prevalence", neg="is.neg", threshold=0.2)
table(contamdf.prev02$contaminant)
write.csv(contamdf.prev02, file = "contamdf.prev02.csv")

# Make phyloseq object of presence-absence in negative controls and true samples
ps.pa <- transform_sample_counts(dna_CVG_phyloseq, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$Sample_or_Control == "Control", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$Sample_or_Control == "Sample", ps.pa)
# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                    contaminant=contamdf.prev$contaminant)
ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

?isContaminant

#using the combined method (with Fishers)
contamdf.comb <- isContaminant(dna_CVG_phyloseq, method="combined", conc="DNA_conc", neg="is.neg", detailed = TRUE)
table(contamdf.comb$contaminant)
contamdf.comb03 <- isContaminant(dna_CVG_phyloseq, method="combined", conc="DNA_conc", neg="is.neg", threshold=0.3, detailed = TRUE)
table(contamdf.comb03$contaminant)

write.csv(contamdf.comb03, file = "contamdf.comb03.csv")


#combined on OTU
OTUcontamdf.comb03 <- isContaminant(all_OTUS_CVG_phyloseq, method="combined", conc="DNA_conc", neg="is.neg", threshold=0.3, detailed = TRUE)
table(OTUcontamdf.comb03$contaminant)

write.csv(withTheTx, file = "OTUcontamdf.comb03.csv")
OTUcontamdf.comb03$OTU <- rownames(OTUcontamdf.comb03)
withTheTx = left_join(OTUcontamdf.comb03, tax_data_table, by = "OTU")

#using the both method (has to satisfy both criteria)
contamdf.both <- isContaminant(dna_CVG_phyloseq, method="both", conc="DNA_conc", neg="is.neg", detailed = TRUE)
table(contamdf.both$contaminant)
head(contamdf.both)


plot_frequency(dna_CVG_phyloseq, taxa_names(dna_CVG_phyloseq)[c(2,3,4)], conc="DNA_conc") + 
  xlab("DNA Concentration (Tapestation or Qubit)") +theme_bw()



##############################################################

#i'm going to look at the histogram of P values

hist(contamdf.comb$p)
ggsave("contamdfcombpspreas.pdf", width = 37, height = 25, units = c('cm'))

hist(contamdf.freq$p)
hist(contamdf.prev$p)
hist(contamdf.prev05$p)

ggplot(data=contamdf.freq, aes(contamdf.freq$p)) + 
  geom_histogram(breaks=seq(0, 1, by = 0.1), 
                 col="#CF3EA4", 
                 fill="#327EC1", 
                 alpha = .8) + 
  labs(title="Frequency-based Decontam: OTU") +
  labs(x="P value", y="Frequency") + theme(
    plot.title = element_text(color="black", size=28, face="bold"),
    axis.title.x = element_text(color="black", size=28, face="bold"),
    axis.title.y = element_text(color="black", size=28, face="bold"),
    axis.text=element_text(size=26),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 20))
ggsave("decontam_frequency_OTU.pdf", width = 25, height = 25, units = c('cm'))

ggplot(data=contamdf.prev, aes(contamdf.prev$p)) + 
  geom_histogram(breaks=seq(0, 1, by = 0.1), 
                 col="#CF3EA4", 
                 fill="#327EC1", 
                 alpha = .8) + 
  labs(title="Prevalence-based Decontam: OTU") +
  labs(x="P value", y="Frequency") + theme(
    plot.title = element_text(color="black", size=28, face="bold"),
    axis.title.x = element_text(color="black", size=28, face="bold"),
    axis.title.y = element_text(color="black", size=28, face="bold"),
    axis.text=element_text(size=26),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 20))
ggsave("decontam_prevalence_OTU.pdf", width = 25, height = 25, units = c('cm'))

ggplot(data=contamdf.comb, aes(contamdf.comb$p)) + 
  geom_histogram(breaks=seq(0, 1, by = 0.1), 
                 col="#CF3EA4", 
                 fill="#327EC1", 
                 alpha = .8) + 
  labs(title="Combined Decontam: OTU") +
  labs(x="P value", y="Frequency") + theme(
    plot.title = element_text(color="black", size=28, face="bold"),
    axis.title.x = element_text(color="black", size=28, face="bold"),
    axis.title.y = element_text(color="black", size=28, face="bold"),
    axis.text=element_text(size=26),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 20))
ggsave("decontam_combined_OTU.pdf", width = 25, height = 25, units = c('cm'))
##############################################################

##############################################################
#going to plot histograms of the DNA input for various samples and OTUs
Temperature <- airquality$Temp
hist(Temperature)


##############################################################
sample_data(No_BadTax_all_OTUS_CVG_phyloseq)$is.neg <- sample_data(No_BadTax_all_OTUS_CVG_phyloseq)$Sample_or_Control == "Control"

#test - separate into tape and qubit!


tape = subset_samples(No_BadTax_all_OTUS_CVG_phyloseq, tape_or_qubit=="tape")
qubit = subset_samples(No_BadTax_all_OTUS_CVG_phyloseq, tape_or_qubit=="qubit")

#yes, now each of those phyloseq object contain only the samples processed that way

#okay, lets look at tape first
#now we are going to use the frequency method to look at contaminants

tapecontamdf.freq <- isContaminant(tape, method="frequency", conc="DNA_conc", threshold=0.2)
head(tapecontamdf.freq)

#how many are 'contaminants'?
table(tapecontamdf.freq$contaminant)
hist(tapecontamdf.freq$p)

#now we are going to use the prevalence method to look at contaminants

tapecontamdf.prev <- isContaminant(tape, method="prevalence", neg="is.neg", threshold=0.2)
table(tapecontamdf.prev$contaminant)
hist(tapecontamdf.prev$p)
#using the combined method (with Fishers)

tapecontamdf.comb <- isContaminant(tape, method="combined", conc="DNA_conc", neg="is.neg", threshold=0.3, detailed = TRUE)
table(tapecontamdf.comb$contaminant)
hist(tapecontamdf.comb$p)
#okay, lets look at qubit first
#now we are going to use the frequency method to look at contaminants

qubitcontamdf.freq <- isContaminant(qubit, method="frequency", conc="DNA_conc", threshold=0.2)
head(qubitcontamdf.freq)

#how many are 'contaminants'?
table(qubitcontamdf.freq$contaminant)
hist(qubitcontamdf.freq$p)
#now we are going to use the prevalence method to look at contaminants

qubitcontamdf.prev <- isContaminant(qubit, method="prevalence", neg="is.neg", threshold=0.2)
table(qubitcontamdf.prev$contaminant)
hist(qubitcontamdf.prev$p)
#using the combined method (with Fishers)

qubitcontamdf.comb <- isContaminant(qubit, method="combined", conc="DNA_conc", neg="is.neg", threshold=0.3, detailed = TRUE)
table(qubitcontamdf.comb$contaminant)
hist(qubitcontamdf.comb$p)

##############################################################
#here i am checking if the 'contaminant' OTUs also appear invindividually in tape and qubit

tapecontamdf.comb["OTU_30",]
tapecontamdf.freq["OTU_30",]
tapecontamdf.prev["OTU_30",]

qubitcontamdf.comb["OTU_30",]
qubitcontamdf.freq["OTU_30",]
qubitcontamdf.prev["OTU_30",]

tapecontamdf.comb["OTU_117",]
tapecontamdf.freq["OTU_117",]
tapecontamdf.prev["OTU_117",]

qubitcontamdf.comb["OTU_117",]
qubitcontamdf.freq["OTU_117",]
qubitcontamdf.prev["OTU_117",]

tapecontamdf.comb["OTU_19",]
tapecontamdf.freq["OTU_19",]
tapecontamdf.prev["OTU_19",]

qubitcontamdf.comb["OTU_19",]
qubitcontamdf.freq["OTU_19",]
qubitcontamdf.prev["OTU_19",]

tapecontamdf.comb["OTU_101",]
tapecontamdf.freq["OTU_101",]
tapecontamdf.prev["OTU_101",]

qubitcontamdf.comb["OTU_101",]
qubitcontamdf.freq["OTU_101",]
qubitcontamdf.prev["OTU_101",]

tapecontamdf.comb["OTU_129",]
tapecontamdf.freq["OTU_129",]
tapecontamdf.prev["OTU_129",]

qubitcontamdf.comb["OTU_129",]
qubitcontamdf.freq["OTU_129",]
qubitcontamdf.prev["OTU_129",]

tapecontamdf.comb["OTU_25",]
tapecontamdf.freq["OTU_25",]
tapecontamdf.prev["OTU_25",]

qubitcontamdf.comb["OTU_25",]
qubitcontamdf.freq["OTU_25",]
qubitcontamdf.prev["OTU_25",]

tapecontamdf.comb["OTU_305",]
tapecontamdf.freq["OTU_305",]
tapecontamdf.prev["OTU_305",]

qubitcontamdf.comb["OTU_305",]
qubitcontamdf.freq["OTU_305",]
qubitcontamdf.prev["OTU_305",]

tapecontamdf.comb["OTU_41",]
tapecontamdf.freq["OTU_41",]
tapecontamdf.prev["OTU_41",]

qubitcontamdf.comb["OTU_41",]
qubitcontamdf.freq["OTU_41",]
qubitcontamdf.prev["OTU_41",]

tapecontamdf.comb["OTU_354",]
tapecontamdf.freq["OTU_354",]
tapecontamdf.prev["OTU_354",]

qubitcontamdf.comb["OTU_354",]
qubitcontamdf.freq["OTU_354",]
qubitcontamdf.prev["OTU_354",]

tapecontamdf.comb["OTU_75",]
tapecontamdf.freq["OTU_75",]
tapecontamdf.prev["OTU_75",]

qubitcontamdf.comb["OTU_75",]
qubitcontamdf.freq["OTU_75",]
qubitcontamdf.prev["OTU_75",]

tapecontamdf.comb["OTU_63",]
tapecontamdf.freq["OTU_63",]
tapecontamdf.prev["OTU_63",]

qubitcontamdf.comb["OTU_63",]
qubitcontamdf.freq["OTU_63",]
qubitcontamdf.prev["OTU_63",]

tapecontamdf.comb["OTU_106",]
tapecontamdf.freq["OTU_106",]
tapecontamdf.prev["OTU_106",]

qubitcontamdf.comb["OTU_106",]
qubitcontamdf.freq["OTU_106",]
qubitcontamdf.prev["OTU_106",]

tapecontamdf.comb["OTU_587",]
tapecontamdf.freq["OTU_587",]
tapecontamdf.prev["OTU_587",]

qubitcontamdf.comb["OTU_587",]
qubitcontamdf.freq["OTU_587",]
qubitcontamdf.prev["OTU_587",]

tapecontamdf.comb["OTU_61",]
tapecontamdf.freq["OTU_61",]
tapecontamdf.prev["OTU_61",]

qubitcontamdf.comb["OTU_61",]
qubitcontamdf.freq["OTU_61",]
qubitcontamdf.prev["OTU_61",]

tapecontamdf.comb["OTU_256",]
tapecontamdf.freq["OTU_256",]
tapecontamdf.prev["OTU_256",]

qubitcontamdf.comb["OTU_256",]
qubitcontamdf.freq["OTU_256",]
qubitcontamdf.prev["OTU_256",]

tapecontamdf.comb["OTU_214",]
tapecontamdf.freq["OTU_214",]
tapecontamdf.prev["OTU_214",]

qubitcontamdf.comb["OTU_214",]
qubitcontamdf.freq["OTU_214",]
qubitcontamdf.prev["OTU_214",]

tapecontamdf.comb["OTU_279",]
tapecontamdf.freq["OTU_279",]
tapecontamdf.prev["OTU_279",]

qubitcontamdf.comb["OTU_279",]
qubitcontamdf.freq["OTU_279",]
qubitcontamdf.prev["OTU_279",]

tapecontamdf.comb["OTU_53",]
tapecontamdf.freq["OTU_53",]
tapecontamdf.prev["OTU_53",]

qubitcontamdf.comb["OTU_53",]
qubitcontamdf.freq["OTU_53",]
qubitcontamdf.prev["OTU_53",]

tapecontamdf.comb["OTU_397",]
tapecontamdf.freq["OTU_397",]
tapecontamdf.prev["OTU_397",]

qubitcontamdf.comb["OTU_397",]
qubitcontamdf.freq["OTU_397",]
qubitcontamdf.prev["OTU_397",]

tapecontamdf.comb["OTU_222",]
tapecontamdf.freq["OTU_222",]
tapecontamdf.prev["OTU_222",]

qubitcontamdf.comb["OTU_222",]
qubitcontamdf.freq["OTU_222",]
qubitcontamdf.prev["OTU_222",]

tapecontamdf.comb["OTU_265",]
tapecontamdf.freq["OTU_265",]
tapecontamdf.prev["OTU_265",]

qubitcontamdf.comb["OTU_265",]
qubitcontamdf.freq["OTU_265",]
qubitcontamdf.prev["OTU_265",]

tapecontamdf.comb["OTU_23",]
tapecontamdf.freq["OTU_23",]
tapecontamdf.prev["OTU_23",]

qubitcontamdf.comb["OTU_23",]
qubitcontamdf.freq["OTU_23",]
qubitcontamdf.prev["OTU_23",]

tapecontamdf.comb["OTU_38",]
tapecontamdf.freq["OTU_38",]
tapecontamdf.prev["OTU_38",]

qubitcontamdf.comb["OTU_38",]
qubitcontamdf.freq["OTU_38",]
qubitcontamdf.prev["OTU_38",]

tapecontamdf.comb["OTU_56",]
tapecontamdf.freq["OTU_56",]
tapecontamdf.prev["OTU_56",]

qubitcontamdf.comb["OTU_56",]
qubitcontamdf.freq["OTU_56",]
qubitcontamdf.prev["OTU_56",]

tapecontamdf.comb["OTU_35",]
tapecontamdf.freq["OTU_35",]
tapecontamdf.prev["OTU_35",]

qubitcontamdf.comb["OTU_35",]
qubitcontamdf.freq["OTU_35",]
qubitcontamdf.prev["OTU_35",]

tapecontamdf.comb["OTU_280",]
tapecontamdf.freq["OTU_280",]
tapecontamdf.prev["OTU_280",]

qubitcontamdf.comb["OTU_280",]
qubitcontamdf.freq["OTU_280",]
qubitcontamdf.prev["OTU_280",]

tapecontamdf.comb["OTU_171",]
tapecontamdf.freq["OTU_171",]
tapecontamdf.prev["OTU_171",]

qubitcontamdf.comb["OTU_171",]
qubitcontamdf.freq["OTU_171",]
qubitcontamdf.prev["OTU_171",]

tapecontamdf.comb["OTU_548",]
tapecontamdf.freq["OTU_548",]
tapecontamdf.prev["OTU_548",]

qubitcontamdf.comb["OTU_548",]
qubitcontamdf.freq["OTU_548",]
qubitcontamdf.prev["OTU_548",]

tapecontamdf.comb["OTU_65",]
tapecontamdf.freq["OTU_65",]
tapecontamdf.prev["OTU_65",]

qubitcontamdf.comb["OTU_65",]
qubitcontamdf.freq["OTU_65",]
qubitcontamdf.prev["OTU_65",]

tapecontamdf.comb["OTU_40",]
tapecontamdf.freq["OTU_40",]
tapecontamdf.prev["OTU_40",]

qubitcontamdf.comb["OTU_40",]
qubitcontamdf.freq["OTU_40",]
qubitcontamdf.prev["OTU_40",]

tapecontamdf.comb["OTU_170",]
tapecontamdf.freq["OTU_170",]
tapecontamdf.prev["OTU_170",]

qubitcontamdf.comb["OTU_170",]
qubitcontamdf.freq["OTU_170",]
qubitcontamdf.prev["OTU_170",]

tapecontamdf.comb["OTU_111",]
tapecontamdf.freq["OTU_111",]
tapecontamdf.prev["OTU_111",]

qubitcontamdf.comb["OTU_111",]
qubitcontamdf.freq["OTU_111",]
qubitcontamdf.prev["OTU_111",]

tapecontamdf.comb["OTU_147",]
tapecontamdf.freq["OTU_147",]
tapecontamdf.prev["OTU_147",]

qubitcontamdf.comb["OTU_147",]
qubitcontamdf.freq["OTU_147",]
qubitcontamdf.prev["OTU_147",]

tapecontamdf.comb["OTU_205",]
tapecontamdf.freq["OTU_205",]
tapecontamdf.prev["OTU_205",]

qubitcontamdf.comb["OTU_205",]
qubitcontamdf.freq["OTU_205",]
qubitcontamdf.prev["OTU_205",]

tapecontamdf.comb["OTU_208",]
tapecontamdf.freq["OTU_208",]
tapecontamdf.prev["OTU_208",]

qubitcontamdf.comb["OTU_208",]
qubitcontamdf.freq["OTU_208",]
qubitcontamdf.prev["OTU_208",]

tapecontamdf.comb["OTU_160",]
tapecontamdf.freq["OTU_160",]
tapecontamdf.prev["OTU_160",]

qubitcontamdf.comb["OTU_160",]
qubitcontamdf.freq["OTU_160",]
qubitcontamdf.prev["OTU_160",]

tapecontamdf.comb["OTU_115",]
tapecontamdf.freq["OTU_115",]
tapecontamdf.prev["OTU_115",]

qubitcontamdf.comb["OTU_115",]
qubitcontamdf.freq["OTU_115",]
qubitcontamdf.prev["OTU_115",]

tapecontamdf.comb["OTU_133",]
tapecontamdf.freq["OTU_133",]
tapecontamdf.prev["OTU_133",]

qubitcontamdf.comb["OTU_133",]
qubitcontamdf.freq["OTU_133",]
qubitcontamdf.prev["OTU_133",]

tapecontamdf.comb["OTU_181",]
tapecontamdf.freq["OTU_181",]
tapecontamdf.prev["OTU_181",]

qubitcontamdf.comb["OTU_181",]
qubitcontamdf.freq["OTU_181",]
qubitcontamdf.prev["OTU_181",]

tapecontamdf.comb["OTU_17",]
tapecontamdf.freq["OTU_17",]
tapecontamdf.prev["OTU_17",]

qubitcontamdf.comb["OTU_17",]
qubitcontamdf.freq["OTU_17",]
qubitcontamdf.prev["OTU_17",]

tapecontamdf.comb["OTU_530",]
tapecontamdf.freq["OTU_530",]
tapecontamdf.prev["OTU_530",]

qubitcontamdf.comb["OTU_530",]
qubitcontamdf.freq["OTU_530",]
qubitcontamdf.prev["OTU_530",]

tapecontamdf.comb["OTU_427",]
tapecontamdf.freq["OTU_427",]
tapecontamdf.prev["OTU_427",]

qubitcontamdf.comb["OTU_427",]
qubitcontamdf.freq["OTU_427",]
qubitcontamdf.prev["OTU_427",]

tapecontamdf.comb["OTU_220",]
tapecontamdf.freq["OTU_220",]
tapecontamdf.prev["OTU_220",]

qubitcontamdf.comb["OTU_220",]
qubitcontamdf.freq["OTU_220",]
qubitcontamdf.prev["OTU_220",]

tapecontamdf.comb["OTU_340",]
tapecontamdf.freq["OTU_340",]
tapecontamdf.prev["OTU_340",]

qubitcontamdf.comb["OTU_340",]
qubitcontamdf.freq["OTU_340",]
qubitcontamdf.prev["OTU_340",]

tapecontamdf.comb["OTU_185",]
tapecontamdf.freq["OTU_185",]
tapecontamdf.prev["OTU_185",]

qubitcontamdf.comb["OTU_185",]
qubitcontamdf.freq["OTU_185",]
qubitcontamdf.prev["OTU_185",]

tapecontamdf.comb["OTU_195",]
tapecontamdf.freq["OTU_195",]
tapecontamdf.prev["OTU_195",]

qubitcontamdf.comb["OTU_195",]
qubitcontamdf.freq["OTU_195",]
qubitcontamdf.prev["OTU_195",]

tapecontamdf.comb["OTU_340",]
tapecontamdf.freq["OTU_340",]
tapecontamdf.prev["OTU_340",]

qubitcontamdf.comb["OTU_340",]
qubitcontamdf.freq["OTU_340",]
qubitcontamdf.prev["OTU_340",]

tapecontamdf.comb["OTU_455",]
tapecontamdf.freq["OTU_455",]
tapecontamdf.prev["OTU_455",]

qubitcontamdf.comb["OTU_455",]
qubitcontamdf.freq["OTU_455",]
qubitcontamdf.prev["OTU_455",]

tapecontamdf.comb["OTU_86",]
tapecontamdf.freq["OTU_86",]
tapecontamdf.prev["OTU_86",]

qubitcontamdf.comb["OTU_86",]
qubitcontamdf.freq["OTU_86",]
qubitcontamdf.prev["OTU_86",]

tapecontamdf.comb["OTU_176",]
tapecontamdf.freq["OTU_176",]
tapecontamdf.prev["OTU_176",]

qubitcontamdf.comb["OTU_176",]
qubitcontamdf.freq["OTU_176",]
qubitcontamdf.prev["OTU_176",]

tapecontamdf.comb["OTU_91",]
tapecontamdf.freq["OTU_91",]
tapecontamdf.prev["OTU_91",]

qubitcontamdf.comb["OTU_91",]
qubitcontamdf.freq["OTU_91",]
qubitcontamdf.prev["OTU_91",]

tapecontamdf.comb["OTU_18",]
tapecontamdf.freq["OTU_18",]
tapecontamdf.prev["OTU_18",]

qubitcontamdf.comb["OTU_18",]
qubitcontamdf.freq["OTU_18",]
qubitcontamdf.prev["OTU_18",]

tapecontamdf.comb["OTU_4",]
tapecontamdf.freq["OTU_4",]
tapecontamdf.prev["OTU_4",]

qubitcontamdf.comb["OTU_4",]
qubitcontamdf.freq["OTU_4",]
qubitcontamdf.prev["OTU_4",]

#########################################
#viewing the contminants coloured by sample biomass

GP20 = prune_taxa("OTU_30", all_OTUS_CVG_phyloseq)
GP.chl = prune_samples(sample_sums(GP20)>=1, GP20)
plot_bar(GP.chl, fill="bioMass")

GP20 = prune_taxa("OTU_117", all_OTUS_CVG_phyloseq)
GP.chl = prune_samples(sample_sums(GP20)>=1, GP20)
plot_bar(GP.chl, fill="bioMass")

GP20 = prune_taxa("OTU_19", dna_CVG_phyloseq)
GP.chl = prune_samples(sample_sums(GP20)>=1, GP20)
plot_bar(GP.chl, fill="bioMass")

GP20 = prune_taxa("OTU_101", all_OTUS_CVG_phyloseq)
GP.chl = prune_samples(sample_sums(GP20)>=1, GP20)
plot_bar(GP.chl, fill="bioMass")

GP20 = prune_taxa("OTU_129", all_OTUS_CVG_phyloseq)
GP.chl = prune_samples(sample_sums(GP20)>=1, GP20)
plot_bar(GP.chl, fill="bioMass")

GP20 = prune_taxa("OTU_25", dna_CVG_phyloseq)
GP.chl = prune_samples(sample_sums(GP20)>=1, GP20)
plot_bar(GP.chl, fill="bioMass")

GP20 = prune_taxa("OTU_305", dna_CVG_phyloseq)
GP.chl = prune_samples(sample_sums(GP20)>=1, GP20)
plot_bar(GP.chl, fill="bioMass")

GP20 = prune_taxa("OTU_41", all_OTUS_CVG_phyloseq)
GP.chl = prune_samples(sample_sums(GP20)>=1, GP20)
plot_bar(GP.chl, fill="bioMass")

GP20 = prune_taxa("OTU_354", dna_CVG_phyloseq)
GP.chl = prune_samples(sample_sums(GP20)>=1, GP20)
plot_bar(GP.chl, fill="bioMass")

GP20 = prune_taxa("OTU_75", dna_CVG_phyloseq)
GP.chl = prune_samples(sample_sums(GP20)>=1, GP20)
plot_bar(GP.chl, fill="bioMass")

GP20 = prune_taxa("OTU_63", dna_CVG_phyloseq)
GP.chl = prune_samples(sample_sums(GP20)>=1, GP20)
plot_bar(GP.chl, fill="bioMass")

GP20 = prune_taxa("OTU_106", dna_CVG_phyloseq)
GP.chl = prune_samples(sample_sums(GP20)>=1, GP20)
plot_bar(GP.chl, fill="bioMass")

GP20 = prune_taxa("OTU_587", dna_CVG_phyloseq)
GP.chl = prune_samples(sample_sums(GP20)>=1, GP20)
plot_bar(GP.chl, fill="bioMass")

GP20 = prune_taxa("OTU_61", dna_CVG_phyloseq)
GP.chl = prune_samples(sample_sums(GP20)>=1, GP20)
plot_bar(GP.chl, fill="bioMass")

GP20 = prune_taxa("OTU_256", dna_CVG_phyloseq)
GP.chl = prune_samples(sample_sums(GP20)>=1, GP20)
plot_bar(GP.chl, fill="bioMass")

GP20 = prune_taxa("OTU_214", dna_CVG_phyloseq)
GP.chl = prune_samples(sample_sums(GP20)>=1, GP20)
plot_bar(GP.chl, fill="bioMass")

GP20 = prune_taxa("OTU_279", dna_CVG_phyloseq)
GP.chl = prune_samples(sample_sums(GP20)>=1, GP20)
plot_bar(GP.chl, fill="bioMass")

GP20 = prune_taxa("OTU_53", dna_CVG_phyloseq)
GP.chl = prune_samples(sample_sums(GP20)>=1, GP20)
plot_bar(GP.chl, fill="bioMass")

GP20 = prune_taxa("OTU_397", dna_CVG_phyloseq)
GP.chl = prune_samples(sample_sums(GP20)>=1, GP20)
plot_bar(GP.chl, fill="bioMass")

GP20 = prune_taxa("OTU_222", dna_CVG_phyloseq)
GP.chl = prune_samples(sample_sums(GP20)>=1, GP20)
plot_bar(GP.chl, fill="bioMass")

GP20 = prune_taxa("OTU_265", all_OTUS_CVG_phyloseq)
GP.chl = prune_samples(sample_sums(GP20)>=1, GP20)
plot_bar(GP.chl, fill="bioMass")

GP20 = prune_taxa("OTU_23", dna_CVG_phyloseq)
GP.chl = prune_samples(sample_sums(GP20)>=1, GP20)
plot_bar(GP.chl, fill="bioMass")

GP20 = prune_taxa("OTU_38", all_OTUS_CVG_phyloseq)
GP.chl = prune_samples(sample_sums(GP20)>=1, GP20)
plot_bar(GP.chl, fill="bioMass")

GP20 = prune_taxa("OTU_56", dna_CVG_phyloseq)
GP.chl = prune_samples(sample_sums(GP20)>=1, GP20)
plot_bar(GP.chl, fill="bioMass")

GP20 = prune_taxa("OTU_35", dna_CVG_phyloseq)
GP.chl = prune_samples(sample_sums(GP20)>=1, GP20)
plot_bar(GP.chl, fill="bioMass")

GP20 = prune_taxa("OTU_280", dna_CVG_phyloseq)
GP.chl = prune_samples(sample_sums(GP20)>=1, GP20)
plot_bar(GP.chl, fill="bioMass")

GP20 = prune_taxa("OTU_171", dna_CVG_phyloseq)
GP.chl = prune_samples(sample_sums(GP20)>=1, GP20)
plot_bar(GP.chl, fill="bioMass")

GP20 = prune_taxa("OTU_548", dna_CVG_phyloseq)
GP.chl = prune_samples(sample_sums(GP20)>=1, GP20)
plot_bar(GP.chl, fill="bioMass")

GP20 = prune_taxa("OTU_65", dna_CVG_phyloseq)
GP.chl = prune_samples(sample_sums(GP20)>=1, GP20)
plot_bar(GP.chl, fill="bioMass")

GP20 = prune_taxa("OTU_40", dna_CVG_phyloseq)
GP.chl = prune_samples(sample_sums(GP20)>=1, GP20)
plot_bar(GP.chl, fill="bioMass")

GP20 = prune_taxa("OTU_170", dna_CVG_phyloseq)
GP.chl = prune_samples(sample_sums(GP20)>=1, GP20)
plot_bar(GP.chl, fill="bioMass")

GP20 = prune_taxa("OTU_111", dna_CVG_phyloseq)
GP.chl = prune_samples(sample_sums(GP20)>=1, GP20)
plot_bar(GP.chl, fill="bioMass")

GP20 = prune_taxa("OTU_147", dna_CVG_phyloseq)
GP.chl = prune_samples(sample_sums(GP20)>=1, GP20)
plot_bar(GP.chl, fill="bioMass")

GP20 = prune_taxa("OTU_205", dna_CVG_phyloseq)
GP.chl = prune_samples(sample_sums(GP20)>=1, GP20)
plot_bar(GP.chl, fill="bioMass")

GP20 = prune_taxa("OTU_208", dna_CVG_phyloseq)
GP.chl = prune_samples(sample_sums(GP20)>=1, GP20)
plot_bar(GP.chl, fill="bioMass")

GP20 = prune_taxa("OTU_160", dna_CVG_phyloseq)
GP.chl = prune_samples(sample_sums(GP20)>=1, GP20)
plot_bar(GP.chl, fill="bioMass")

GP20 = prune_taxa("OTU_115", dna_CVG_phyloseq)
GP.chl = prune_samples(sample_sums(GP20)>=1, GP20)
plot_bar(GP.chl, fill="bioMass")

GP20 = prune_taxa("OTU_133", dna_CVG_phyloseq)
GP.chl = prune_samples(sample_sums(GP20)>=1, GP20)
plot_bar(GP.chl, fill="bioMass")

GP20 = prune_taxa("OTU_181", dna_CVG_phyloseq)
GP.chl = prune_samples(sample_sums(GP20)>=1, GP20)
plot_bar(GP.chl, fill="bioMass")

GP20 = prune_taxa("OTU_17", dna_CVG_phyloseq)
GP.chl = prune_samples(sample_sums(GP20)>=1, GP20)
plot_bar(GP.chl, fill="bioMass")

GP20 = prune_taxa("OTU_530", dna_CVG_phyloseq)
GP.chl = prune_samples(sample_sums(GP20)>=1, GP20)
plot_bar(GP.chl, fill="bioMass")

GP20 = prune_taxa("OTU_427", dna_CVG_phyloseq)
GP.chl = prune_samples(sample_sums(GP20)>=1, GP20)
plot_bar(GP.chl, fill="bioMass")

GP20 = prune_taxa("OTU_220", all_OTUS_CVG_phyloseq)
GP.chl = prune_samples(sample_sums(GP20)>=1, GP20)
plot_bar(GP.chl, fill="bioMass")

GP20 = prune_taxa("OTU_340", dna_CVG_phyloseq)
GP.chl = prune_samples(sample_sums(GP20)>=1, GP20)
plot_bar(GP.chl, fill="bioMass")

GP20 = prune_taxa("OTU_185", dna_CVG_phyloseq)
GP.chl = prune_samples(sample_sums(GP20)>=1, GP20)
plot_bar(GP.chl, fill="bioMass")

GP20 = prune_taxa("OTU_195", dna_CVG_phyloseq)
GP.chl = prune_samples(sample_sums(GP20)>=1, GP20)
plot_bar(GP.chl, fill="bioMass")

GP20 = prune_taxa("OTU_340", dna_CVG_phyloseq)
GP.chl = prune_samples(sample_sums(GP20)>=1, GP20)
plot_bar(GP.chl, fill="bioMass")

GP20 = prune_taxa("OTU_455", dna_CVG_phyloseq)
GP.chl = prune_samples(sample_sums(GP20)>=1, GP20)
plot_bar(GP.chl, fill="bioMass")

GP20 = prune_taxa("OTU_86", dna_CVG_phyloseq)
GP.chl = prune_samples(sample_sums(GP20)>=1, GP20)
plot_bar(GP.chl, fill="bioMass")

GP20 = prune_taxa("OTU_176", dna_CVG_phyloseq)
GP.chl = prune_samples(sample_sums(GP20)>=1, GP20)
plot_bar(GP.chl, fill="bioMass")

GP20 = prune_taxa("OTU_91", dna_CVG_phyloseq)
GP.chl = prune_samples(sample_sums(GP20)>=1, GP20)
plot_bar(GP.chl, fill="bioMass")

GP20 = prune_taxa("OTU_18", dna_CVG_phyloseq)
GP.chl = prune_samples(sample_sums(GP20)>=1, GP20)
plot_bar(GP.chl, fill="bioMass")

GP20 = prune_taxa("OTU_4", dna_CVG_phyloseq)
GP.chl = prune_samples(sample_sums(GP20)>=1, GP20)
plot_bar(GP.chl, fill="bioMass")

GP20 = prune_taxa("OTU_106", genus_dna_CVG_phyloseq)
GP.chl = prune_samples(sample_sums(GP20)>=1, GP20)
plot_bar(GP.chl)

GP20 = subset_taxa(genus_dna_CVG_phyloseq, Genus=="Micrococcus")

GP20 = subset_taxa(CVG_genus_phyloseq, Genus=="Citrobacter")
GP.chl = prune_samples(sample_sums(GP20)>=1, GP20)
plot_bar(GP.chl)

GP20 = subset_taxa(CVG_genus_phyloseq, Genus=="Fusobacterium")
GP.chl = prune_samples(sample_sums(GP20)>=1, GP20)
plot_bar(GP.chl)

GP20 = subset_taxa(CVG_genus_phyloseq, Genus=="Granulicatella")
GP.chl = prune_samples(sample_sums(GP20)>=1, GP20)
plot_bar(GP.chl)

GP20 = subset_taxa(CVG_genus_phyloseq, Genus=="Zymomonas")
GP.chl = prune_samples(sample_sums(GP20)>=1, GP20)
plot_bar(GP.chl)