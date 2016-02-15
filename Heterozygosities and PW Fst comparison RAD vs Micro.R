
install.packages("hierfstat")

##### COMPARISON OF RAD AND MICROSATELLITE HETEROZYGOSITIES ##########

library(adegenet)
library(hierfstat)
library(ggplot2)
library(reshape2)


### DATASETS ---------------------------------------------------------------------------------------------

RADhets <- read.delim("/media//dan//34D5D1CE642D7E36//2013076_Hanfling_Bernd/Stacks/Stacks_analyses_TRIMMED/Pure_cru_M_2/populations_r07_p17_Ho_filetered/Obs_hets_Variant_only.txt") ## RAD heterozygosities

M2_temp <- read.fstat("~/../Dropbox/PhD/Dans_PhD_Shared/Data/RAD/RADsample_microsat_genotypes_1.dat") ## M2 hets
M2_pops <- read.delim("~/../Dropbox/PhD/Dans_PhD_Shared/Data/RAD/RADsample_microsat_POPS.txt", header = F)

pop(M2_temp) <- M2_pops$V1 ## attach names

## subset microsat data to only shared RAD pops
SepM2 <- seppop(M2_temp)
names(SepM2)

M2 <- repool(SepM2$BF, SepM2$CAKE,SepM2$"CA-LK", SepM2$COP,SepM2$MOAT, SepM2$OBY, SepM2$OU, SepM2$PED , SepM2$POLEN, SepM2$PRO,  SepM2$SD, SepM2$SK, SepM2$STEC, SepM2$STYV, SepM2$TROM, SepM2$TU, SepM2$WEN )
pop(M2)


## Read in heterozygosities for pops shared between RAD and Microsats. (Doesn't include BOR or V populations in RAD)
both <- read.delim("~/../Dropbox/PhD/Dans_PhD_Shared/Papers/Phylogeography\ paper/Final_data_files/Micro_vs_RAD_Hobs.txt", header = T, sep = " ")



# Obseved Heterozygosities comparisons -------------------------------------------------------------

## First clalculate M2 heterozygosities

summary(M2)
M2_hier <- genind2hierfstat(M2) ## convert to hierfstat format

M2_hier_stats <- basic.stats(M2_hier) 
M2_Hobs <-  M2_hier_stats$Ho  ## Hobs are in here
M2_Hobs_df <- as.data.frame(M2_Hobs) # Add names
names(M2_Hobs_df) <- levels(pop(M2))

M2_Hobs_means<- data.frame(M2 = colMeans(M2_Hobs_df, na.rm = T)) ## Get mean of each population 

both$M2 <- M2_Hobs_means$M2 ## add to the other data
M2_M3_RAD  <- both
M2_M3_RAD

## RAD / M2

###  Pearson's correlation
RAD_M2_Het_corr <- cor.test(M2_M3_RAD$RAD, M2_M3_RAD$M2)
RAD_M2_Het_corr$estimate

## Paired t-test
RAD_M2_Het_T_test <- t.test(M2_M3_RAD$RAD, M2_M3_RAD$M2, paired = T)
RAD_M2_Het_T_test


## RAD / M3 

###  Pearson's correlation
RAD_M3_Het_corr <- cor.test(M2_M3_RAD$RAD, M2_M3_RAD$M3)
RAD_M3_Het_corr

## Paired t-test
RAD_M3_Het_T_test <- t.test(M2_M3_RAD$RAD, M2_M3_RAD$M3, paired = T)
RAD_M3_Het_T_test


## M2 / M3 


###  Pearson's correlation
M2_M3_Het_corr <- cor.test(both$M2, both$M3)
M2_M3_Het_corr

## Paired t-test
M2_M3_Het_T_test <- t.test(both$M2, both$M3, paired = T)
M2_M3_Het_T_test
mean(both$M3)



# FST ------------------------------------------------------------------------------------------------

### First calculate PW fsts for M2 dataset

M2_ppFsts <- pp.fst(M2_hier) ## calculate PW fsts (Weir & Cockerham 1984)
M2names <- levels(pop(M2)) ## get pop names
M2names

## add pop names to matirx
M2_fsts <- M2_ppFsts$fst.pp
M2_fsts.df <- as.data.frame(M2_fsts)
names(M2_fsts.df) = M2names
row.names(M2_fsts.df) <- M2names

M2_fsts_vect <- as.vector(t(M2_fsts.df)) ## make as vector

## make vector lables

names1 <- rep(M2names,16) 
names2 <- paste(c(rep("BF", 17),rep("CAKE", 17),rep("CA.LK", 17),rep("COP", 17),rep("MOAT", 17),rep("OBY", 17),rep("OU", 17),rep("PED", 17),rep("POLEN", 17),rep("PRO", 17),rep("SD", 17),rep("SK", 17),rep("STEC", 17),rep("STYV", 17),rep("TROM", 17),rep("TU", 17), rep("WEN", 17)))
pasted_names <- paste(names2, names1, sep ="_")

lablled_vector <- data.frame(pasted_names, M2_fsts_vect)

write.table(lablled_vector, "~/../Desktop/M2_fsts_vect", sep = "\t") ## reformatted by hand so all PW comparisons match up

M2_M3_RAD_PWfsts <- read.delim("~/../Dropbox/PhD/Dans_PhD_Shared/Papers/Phylogeography paper/Final_data_files/M2_M3_RAD_PWfsts.csv", sep = "\t")

M2_M3_RAD_PWfsts


# FST comparisons ----------------------------------------------------------------------------------

#####T-tests 

?t.test
## M2 - M3 ##

M2_M3 <- t.test(M2_M3_RAD_PWfsts$M2_fsts_vect, M2_M3_RAD_PWfsts$M3, paired = T)
M2_M3

## M2 - RAD ##
M2_RAD <- t.test(M2_M3_RAD_PWfsts$M2_fsts_vect, M2_M3_RAD_PWfsts$RAD, paired = T)
M2_RAD
## M3 - RAD ##
M3_RAD <- t.test(M2_M3_RAD_PWfsts$M3, M2_M3_RAD_PWfsts$RAD, paired = T)
M3_RAD

## Pearsons correlation tests

RAD_M2_Test <- cor.test(M2_M3_RAD_PWfsts$RAD, M2_M3_RAD_PWfsts$M2_fsts_vect)
RAD_M3_Test <- cor.test(M2_M3_RAD_PWfsts$RAD, M2_M3_RAD_PWfsts$M3)
M2_M3_Test <- cor.test(M2_M3_RAD_PWfsts$M2_fsts_vect, M2_M3_RAD_PWfsts$M3)

RAD_M2_Test
RAD_M3_Test
M2_M3_Test


## PLOTS --------------------------------------------------------------------------------------

# ---- Heterozygosities ------

## RAD / M3 

plot(both$RAD, both$Micro, pch = 16, xlab = "RADseq H obs", ylab = "Microsatellites H obs")
text(both$RAD, both$Micro, both$Pop, pos = 4, cex = .7)




M2_M3_RAD_PWfsts.melted <- melt(M2_M3_RAD_PWfsts)

gplot <- ggplot(M2_M3_RAD_PWfsts.melted, aes(pasted_names, value, fill = variable)) + geom_bar(position="dodge", stat = "identity")
gplot + theme(axis.text.x=element_text(angle=-90))



## -----Allele Counts ---------------------------------------------------------------------------------------

# Allele counts per pop
setwd("/media/dan/34D5D1CE642D7E36/2013076_Hanfling_Bernd/Stacks/Stacks_analyses_TRIMMED/Pure_cru_M_2/populations_r07_p17_Ho_filetered")


RADdata <- read.PLINK('plink.raw', , chunkSize=1000, parallel = TRUE, n.cores=3, saveNbAlleles=T)
pops <- read.delim("../populations.txt", header = F)

pop(RADdata) <- pops$V1

RADdata_sep <- seppop(RADdata)
names(RADdata_sep)

levels(pop(RADdata))

allele_SUMs <- data.frame(BF = glSum(RADdata_sep$BF), BOR = glSum(RADdata_sep$BOR), CAKE = glSum(RADdata_sep$CAKE), CALK = glSum(RADdata_sep$CALK), COP = glSum(RADdata_sep$COP), MOAT = glSum(RADdata_sep$MOAT), OBY = glSum(RADdata_sep$OBY), OU = glSum(RADdata_sep$OU), PED = glSum(RADdata_sep$PED), POLEN = glSum(RADdata_sep$POLEN), PRO = glSum(RADdata_sep$PRO), SD = glSum(RADdata_sep$SD), SK = glSum(RADdata_sep$SK), STEC = glSum(RADdata_sep$STEC), STYV = glSum(RADdata_sep$STYV), TROM = glSum(RADdata_sep$TROM), TU = glSum(RADdata_sep$TU), V = glSum(RADdata_sep$V), WEN = glSum(RADdata_sep$WEN))
sum(allele_SUM$Allele_counts)
Allele_counts <- colSums(allele_SUMs)




