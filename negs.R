negs <- read.csv(file = "negs.csv", stringsAsFactors = F) #import otu table


library(tidyverse)
negs$counts <- NA

negs$counts <- rowSums(negs[, 2:8])   
negs <- dplyr::filter(negs, negs$counts > 0)


write.csv(negs, file = "negsandblank.csv")

#subset just the otus in neg 1
ntc1freqs <- subset(negs, select = -c(BLANK, NTC3, NTC4, NTC5, NTC6, NTC7))
ntc1freqs <- ntc1freqs[!(ntc1freqs$NTC1==0),]

#subset just the OTUs in the surrounding samples
s4159 <- subset(otu_data_table, select = c(OTU, CVG14V54159))
s4159 <- s4159[!(s4159$CVG14V54159==0),]
