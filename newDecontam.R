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
library("patchwork")


pop_taxa = function(physeq, badTaxa){
  allTaxa = taxa_names(physeq)
  allTaxa <- allTaxa[!(allTaxa %in% badTaxa)]
  return(prune_taxa(allTaxa, physeq))
}
#Which you would then use as

badTaxa = c("OTU_180", "OTU_268", "OTU_344" , "OTU_354", "OTU_360","OTU_372")
No_BadTax_all_OTUS_CVG_phyloseq = pop_taxa(all_OTUS_CVG_phyloseq, badTaxa)
sample_names(No_BadTax_all_OTUS_CVG_phyloseq)


all_OTUS_CVG_phyloseq <- merge_phyloseq(contam_CVG_otu_table, contam_CVG_tax_table, dna_CVG_sample_table) #allll OTUs
genus_dna_CVG_phyloseq <-tax_glom(all_OTUS_CVG_phyloseq, taxrank = "Genus") #all OTUs by genus

No_BadTax_all_OTUS_genus_phyloseq  <-tax_glom(No_BadTax_all_OTUS_CVG_phyloseq, taxrank = "Genus") #all OTUs by genus no badtax

############## freq 0.1 in OTUs

#now we are going to use the frequency method to look at contaminants

OTU0.1contamdf.freq <- isContaminant(No_BadTax_all_OTUS_CVG_phyloseq, method="frequency", conc="DNA_conc", detailed = TRUE)
head(OTU0.1contamdf.freq)
table(OTU0.1contamdf.freq$contaminant)
OTU0.1contamdf.freq$OTU <- rownames(OTU0.1contamdf.freq)

hist(OTU0.1contamdf.freq$p)

withTheTx = left_join(OTU0.1contamdf.freq, tax_data_table, by = "OTU")
write.csv(withTheTx, file = "finalOTU_decontam_freq0.1.csv")

ggplot(data=OTU0.1contamdf.freq, aes(OTU0.1contamdf.freq$p)) + 
  geom_histogram(breaks=seq(0, 1, by = 0.1), 
                 col="#CF3EA4", 
                 fill="#327EC1", 
                 alpha = .8) + 
  labs(title="Frequency-based Decontam: OTU") +
  labs(x="P value", y="Frequency")+theme_bw() + theme(
    plot.title = element_text(color="black", size=22, face="bold"),
    axis.title.x = element_text(color="black", size=22, face="bold"),
    axis.title.y = element_text(color="black", size=22, face="bold"),
    axis.text=element_text(size=20),
    legend.title = element_text(color = "black", size = 20),
    legend.text = element_text(color = "black", size = 20)) 
ggsave("decontam_frequency_OTU_new.pdf", width = 25, height = 25, units = c('cm'))


############## prev 0.1 in OTUs
 

OTU0.1contamdf.prev <- isContaminant(No_BadTax_all_OTUS_CVG_phyloseq, method="prevalence", neg="is.neg", detailed = TRUE)
head(OTU0.1contamdf.prev)
table(OTU0.1contamdf.prev$contaminant)
OTU0.1contamdf.prev$OTU <- rownames(OTU0.1contamdf.prev)

hist(OTU0.1contamdf.prev$p)

withTheTx = left_join(OTU0.1contamdf.prev, tax_data_table, by = "OTU")
write.csv(withTheTx, file = "finalOTU_decontam_prev0.1.csv")

ggplot(data=OTU0.1contamdf.prev, aes(OTU0.1contamdf.prev$p)) + 
  geom_histogram(breaks=seq(0, 1, by = 0.1), 
                 col="#CF3EA4", 
                 fill="#327EC1", 
                 alpha = .8) + 
  labs(title="Prevalence-based Decontam: OTU") +
  labs(x="P value", y="Frequency")+theme_bw() + theme(
    plot.title = element_text(color="black", size=22, face="bold"),
    axis.title.x = element_text(color="black", size=22, face="bold"),
    axis.title.y = element_text(color="black", size=22, face="bold"),
    axis.text=element_text(size=20),
    legend.title = element_text(color = "black", size = 20),
    legend.text = element_text(color = "black", size = 20)) 
ggsave("decontam_prevalence_OTU_new.pdf", width = 25, height = 25, units = c('cm'))

############## comb 0.3 in OTUs
OTU0.3contamdf.comb <- isContaminant(No_BadTax_all_OTUS_CVG_phyloseq, method="comb", neg="is.neg", conc="DNA_conc", detailed = TRUE, threshold=0.3)
head(OTU0.3contamdf.comb)
table(OTU0.3contamdf.comb$contaminant)
OTU0.3contamdf.comb$OTU <- rownames(OTU0.3contamdf.comb)

hist(OTU0.3contamdf.comb$p)

withTheTx = left_join(OTU0.3contamdf.comb, tax_data_table, by = "OTU")
write.csv(withTheTx, file = "finalOTU_decontam_comb0.3.csv")

ggplot(data=OTU0.3contamdf.comb, aes(OTU0.3contamdf.comb$p)) + 
  geom_histogram(breaks=seq(0, 1, by = 0.1), 
                 col="#CF3EA4", 
                 fill="#327EC1", 
                 alpha = .8) + 
  labs(title="Combined Decontam: OTU") +
  labs(x="P value", y="Frequency")+theme_bw() + theme(
    plot.title = element_text(color="black", size=22, face="bold"),
    axis.title.x = element_text(color="black", size=22, face="bold"),
    axis.title.y = element_text(color="black", size=22, face="bold"),
    axis.text=element_text(size=20),
    legend.title = element_text(color = "black", size = 20),
    legend.text = element_text(color = "black", size = 20)) 
ggsave("decontam_combined_OTU_new.pdf", width = 25, height = 25, units = c('cm'))

######### individ

GP20 = prune_taxa("OTU_598", No_BadTax_all_OTUS_CVG_phyloseq)
GP.chl = prune_samples(sample_sums(GP20)>=1, GP20)
plot_bar(GP.chl, fill="bioMass")

tapecontamdf.comb["OTU_598",]
tapecontamdf.freq["OTU_598",]
tapecontamdf.prev["OTU_598",]

qubitcontamdf.comb["OTU_598",]
qubitcontamdf.freq["OTU_598",]
qubitcontamdf.prev["OTU_598",]

GP20 = prune_taxa("OTU_42", No_BadTax_all_OTUS_CVG_phyloseq)
GP.chl = prune_samples(sample_sums(GP20)>=1, GP20)
plot_bar(GP.chl, fill="bioMass")

tapecontamdf.comb["OTU_42",]
tapecontamdf.freq["OTU_42",]
tapecontamdf.prev["OTU_42",]

qubitcontamdf.comb["OTU_42",]
qubitcontamdf.freq["OTU_42",]
qubitcontamdf.prev["OTU_42",]

GP20 = prune_taxa("OTU_439", No_BadTax_all_OTUS_CVG_phyloseq)
GP.chl = prune_samples(sample_sums(GP20)>=1, GP20)
plot_bar(GP.chl, fill="bioMass")

tapecontamdf.comb["OTU_439",]
tapecontamdf.freq["OTU_439",]
tapecontamdf.prev["OTU_439",]

qubitcontamdf.comb["OTU_439",]
qubitcontamdf.freq["OTU_439",]
qubitcontamdf.prev["OTU_439",]

GP20 = prune_taxa("OTU_137", No_BadTax_all_OTUS_CVG_phyloseq)
GP.chl = prune_samples(sample_sums(GP20)>=1, GP20)
plot_bar(GP.chl, fill="bioMass")

tapecontamdf.comb["OTU_137",]
tapecontamdf.freq["OTU_137",]
tapecontamdf.prev["OTU_137",]

qubitcontamdf.comb["OTU_137",]
qubitcontamdf.freq["OTU_137",]
qubitcontamdf.prev["OTU_137",]

GP20 = prune_taxa("OTU_156", No_BadTax_all_OTUS_CVG_phyloseq)
GP.chl = prune_samples(sample_sums(GP20)>=1, GP20)
plot_bar(GP.chl, fill="bioMass")

tapecontamdf.comb["OTU_156",]
tapecontamdf.freq["OTU_156",]
tapecontamdf.prev["OTU_156",]

qubitcontamdf.comb["OTU_156",]
qubitcontamdf.freq["OTU_156",]
qubitcontamdf.prev["OTU_156",]

GP20 = prune_taxa("OTU_188", No_BadTax_all_OTUS_CVG_phyloseq)
GP.chl = prune_samples(sample_sums(GP20)>=1, GP20)
plot_bar(GP.chl, fill="bioMass")

tapecontamdf.comb["OTU_188",]
tapecontamdf.freq["OTU_188",]
tapecontamdf.prev["OTU_188",]

qubitcontamdf.comb["OTU_188",]
qubitcontamdf.freq["OTU_188",]
qubitcontamdf.prev["OTU_188",]

GP20 = prune_taxa("OTU_260", No_BadTax_all_OTUS_CVG_phyloseq)
GP.chl = prune_samples(sample_sums(GP20)>=1, GP20)
plot_bar(GP.chl, fill="bioMass")

tapecontamdf.comb["OTU_260",]
tapecontamdf.freq["OTU_260",]
tapecontamdf.prev["OTU_260",]

qubitcontamdf.comb["OTU_260",]
qubitcontamdf.freq["OTU_260",]
qubitcontamdf.prev["OTU_260",]

GP20 = prune_taxa("OTU_303", No_BadTax_all_OTUS_CVG_phyloseq)
GP.chl = prune_samples(sample_sums(GP20)>=1, GP20)
plot_bar(GP.chl, fill="bioMass")

tapecontamdf.comb["OTU_303",]
tapecontamdf.freq["OTU_303",]
tapecontamdf.prev["OTU_303",]

qubitcontamdf.comb["OTU_303",]
qubitcontamdf.freq["OTU_303",]
qubitcontamdf.prev["OTU_303",]

GP20 = prune_taxa("OTU_413", No_BadTax_all_OTUS_CVG_phyloseq)
GP.chl = prune_samples(sample_sums(GP20)>=1, GP20)
plot_bar(GP.chl, fill="bioMass")

tapecontamdf.comb["OTU_413",]
tapecontamdf.freq["OTU_413",]
tapecontamdf.prev["OTU_413",]

qubitcontamdf.comb["OTU_413",]
qubitcontamdf.freq["OTU_413",]
qubitcontamdf.prev["OTU_413",]

GP20 = prune_taxa("OTU_76", No_BadTax_all_OTUS_CVG_phyloseq)
GP.chl = prune_samples(sample_sums(GP20)>=1, GP20)
plot_bar(GP.chl, fill="bioMass")

tapecontamdf.comb["OTU_76",]
tapecontamdf.freq["OTU_76",]
tapecontamdf.prev["OTU_76",]

qubitcontamdf.comb["OTU_76",]
qubitcontamdf.freq["OTU_76",]
qubitcontamdf.prev["OTU_76",]

GP20 = prune_taxa("OTU_223", No_BadTax_all_OTUS_CVG_phyloseq)
GP.chl = prune_samples(sample_sums(GP20)>=1, GP20)
plot_bar(GP.chl, fill="bioMass")

tapecontamdf.comb["OTU_223",]
tapecontamdf.freq["OTU_223",]
tapecontamdf.prev["OTU_223",]

qubitcontamdf.comb["OTU_223",]
qubitcontamdf.freq["OTU_223",]
qubitcontamdf.prev["OTU_223",]
##################### switch it up!!!! now for genus

############## freq 0.2 in genus

#now we are going to use the frequency method to look at contaminants

genus0.2contamdf.freq <- isContaminant(No_BadTax_all_OTUS_genus_phyloseq, method="frequency", conc="DNA_conc", detailed = TRUE, threshold=0.2)
head(genus0.2contamdf.freq)
table(genus0.2contamdf.freq$contaminant)
genus0.2contamdf.freq$OTU <- rownames(genus0.2contamdf.freq)

hist(genus0.2contamdf.freq$p)

withTheTx = left_join(genus0.2contamdf.freq, tax_data_table, by = "OTU")
write.csv(withTheTx, file = "finalgenus_decontam_freq0.2.csv")

p2 = ggplot(data=genus0.2contamdf.freq, aes(genus0.2contamdf.freq$p)) + 
  geom_histogram(breaks=seq(0, 1, by = 0.1), 
                 col="#CF3EA4", 
                 fill="#327EC1", 
                 alpha = .8) + 
  labs(title="Frequency-based Decontam: Genus") +
  labs(x="P value", y="Frequency")+theme_bw() + theme(
    plot.title = element_text(color="black", size=20, face="bold"),
    axis.title.x = element_text(color="black", size=20, face="bold"),
    axis.title.y = element_text(color="black", size=20, face="bold"),
    axis.text=element_text(size=20),
    legend.title = element_text(color = "black", size = 20),
    legend.text = element_text(color = "black", size = 20)) 
ggsave("decontam_frequency_genus_new.pdf", width = 25, height = 25, units = c('cm'))


############## prev 0.2 in OTUs
sample_data(No_BadTax_all_OTUS_genus_phyloseq)$is.neg <- sample_data(No_BadTax_all_OTUS_genus_phyloseq)$Sample_or_Control == "Control"


genus0.2contamdf.prev <- isContaminant(No_BadTax_all_OTUS_genus_phyloseq, method="prevalence", neg="is.neg", detailed = TRUE, threshold=0.2)
head(genus0.2contamdf.prev)
table(genus0.2contamdf.prev$contaminant)
genus0.2contamdf.prev$OTU <- rownames(genus0.2contamdf.prev)

hist(genus0.2contamdf.prev$p)

withTheTx = left_join(genus0.2contamdf.prev, tax_data_table, by = "OTU")
write.csv(withTheTx, file = "finalgenus_decontam_prev0.2.csv")

p3 = ggplot(data=genus0.2contamdf.prev, aes(genus0.2contamdf.prev$p)) + 
  geom_histogram(breaks=seq(0, 1, by = 0.1), 
                 col="#CF3EA4", 
                 fill="#327EC1", 
                 alpha = .8) + 
  labs(title="Prevalence-based Decontam: Genus") +
  labs(x="P value", y="Frequency")+theme_bw() + theme(
    plot.title = element_text(color="black", size=20, face="bold"),
    axis.title.x = element_text(color="black", size=20, face="bold"),
    axis.title.y = element_text(color="black", size=20, face="bold"),
    axis.text=element_text(size=20),
    legend.title = element_text(color = "black", size = 20),
    legend.text = element_text(color = "black", size = 20)) 
ggsave("decontam_prevalence_genus_new.pdf", width = 25, height = 25, units = c('cm'))

############## comb 0.3 in OTUs
genus0.3contamdf.comb <- isContaminant(No_BadTax_all_OTUS_genus_phyloseq, method="prevalence", neg="is.neg", conc="DNA_conc", detailed = TRUE, threshold=0.3)
head(genus0.3contamdf.comb)
table(genus0.3contamdf.comb$contaminant)
genus0.3contamdf.comb$OTU <- rownames(genus0.3contamdf.comb)

hist(genus0.3contamdf.comb$p)

withTheTx = left_join(genus0.3contamdf.comb, tax_data_table, by = "OTU")
write.csv(withTheTx, file = "finalgenus_decontam_comb0.3.csv")

p1 = ggplot(data=genus0.3contamdf.comb, aes(genus0.3contamdf.comb$p)) + 
  geom_histogram(breaks=seq(0, 1, by = 0.1), 
                 col="#CF3EA4", 
                 fill="#327EC1", 
                 alpha = .8) + 
  labs(title="Combined Decontam: Genus") +
  labs(x="P value", y="Frequency")+theme_bw() + theme(
    plot.title = element_text(color="black", size=20, face="bold"),
    axis.title.x = element_text(color="black", size=20, face="bold"),
    axis.title.y = element_text(color="black", size=20, face="bold"),
    axis.text=element_text(size=20),
    legend.title = element_text(color = "black", size = 20),
    legend.text = element_text(color = "black", size = 20)) 
ggsave("decontam_combined_genus_new.pdf", width = 25, height = 25, units = c('cm'))


decontamPLot = p1 / (p2 | p3)
ggsave("decontam_p_plots.pdf", width = 30, height = 25, units = c('cm'))


