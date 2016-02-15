
##### COMPARISON OF RAD AND MICROSATELLITE HETEROZYGOSITIES ##########

############# RAD obs. hets ###################

### All RAD hets ###

RADhets <- read.delim("/media//dan//34D5D1CE642D7E36//2013076_Hanfling_Bernd/Stacks/Stacks_analyses_TRIMMED/Pure_cru_M_2/populations_r07_p17_Ho_filetered/Obs_hets_Variant_only.txt")

## quick check
hetsbar <- barplot(RADhets$Obs.Het)
axis(1, RADhets$X..Pop.ID, at = hetsbar, las = 2) ## all looks as expected

## Read in heterozygosities for pops shared between RAD and Microsats. (Doesn't include BOR or V populations in RAD)
both <- read.delim("~/Dropbox/PhD/Dans_PhD_Shared/Papers/Phylogeography\ paper/Final_data_files/Micro_vs_RAD_Hobs.txt", header = T, sep = " ")
both
plot(both$RAD, both$Micro, pch = 16, xlab = "RADseq H obs", ylab = "Microsatellites H obs")
text(both$RAD, both$Micro, both$Pop, pos = 4, cex = .7)


# -----------------------------------------------------------------------------------------------------


########### Misc Microsat summary stats #####################

setwd("~/Dropbox/PhD/Dans_PhD_Shared/Data/Microsatellites/DAPC/Complete dataset outputs")

micro_all_pure <- read.fstat("PurecruplusManc_stripped_and_checked_temp.DAT")

micro_all_pure

sumstats <- summary(micro_all_pure)

par(mfrow=c(2,2))
plot(sumstats$pop.eff,sumstats$pop.nall,xlab="Colonies sample size",ylab="Number of alleles",main="Alleles numbers and sample sizes")
regres <- lm(sumstats$pop.nall~sumstats$pop.eff)
abline(regres)
text(sumstats$pop.eff,sumstats$pop.nall,lab=names(sumstats$pop.eff))

barplot(sumstats$loc.nall,ylab="Number of alleles", main="Number of alleles per locus")
barplot(sumstats$Hexp-sumstats$Hobs,main="Heterozygosity: expected-observed",ylab="Hexp - Hobs")
barplot(sumstats$pop.eff,main="Sample sizes per population",ylab="Number of genotypes",las=3)

sumstats




###################################################################################################







