
library(adegenet)

setwd("/media/dan//34D5D1CE642D7E36//2013076_Hanfling_Bernd/Stacks/Stacks_analyses_TRIMMED/Pure_cru_M_2/populations_r07_p17/")

RAD <- read.genepop('batch_1.gen') ## same data as in the DAPC but in genpop format
RADpops <- read.delim("populations.txt", header = F)
pop(RAD) <- RADpops$V1
pop(RAD) ## OK!


RADsub <- RAD[c(1:9,15:117,128:155,162:170),]  ## remove BOR, STEC and V
#length(indNames(RADsub))

RADsubgen <- genind2genpop(RADsub) ## need to convert to genpop object

MicnRADcoords <- RADcoordsall[c(1,3:5,15,17:21,24,25,27:29,32),] ## coordinates for these pops only


RADDgen <- dist.genpop(RADsubgen,method=2) ##Generates the genetic distance semi-matrix

RADDgeo <- dist(MicnRADcoords) ## Generates the Euclidian distance semi matrix. Note- vignette says to use "dist(body$other$xy),I dont know why but this wont work for me!

RADibd <- mantel.randtest(RADDgen,RADDgeo)

RADibd  ## gives the output for the test  NOTE that the "observation" is the R value. Sqaure to get R^2 (obviously)

plot(RADibd) ##plots histogram fgor test


################## MICROSATS ##########################

M1 <- read.fstat("~/Dropbox/PhD/Dan's PhD (Shared)/Data/Microsatellites/Purecruchecked.DAT")

MicroCoords <- read.csv("~/Dropbox/PhD/Dan's PhD (Shared)/Data/Microsatellites/R stats/Adegenet IBD/Mantelcoordinates.csv",header = TRUE, row.names = 1) ##table containing xy coordinates

############# All populations (M1) #####################

M1@other <- (MicroCoords) ## Attach coordinates to the GENIND file

M1gen <- genind2genpop(M1)

M1Dgen <- dist.genpop(M1gen,method=2) ##Generates the genetic distance semi-matrix

M1Dgeo <- dist(M1gen$other) ## Generates the Euclidian distance semi matrix. Note- vignette says to use "dist(body$other$xy),I dont know why but this wont work for me!

M1ibd <- mantel.randtest(M1Dgen,M1Dgeo)

M1ibd  ## gives the output for the test  NOTE that the "observation" is the R value. Sqaure to get R^2 (obviously)

plot(ibd) ##plots histogram fgor test


popnames <- levels(pop(M1)) ## get names from full microsat
M1_Hap1_noSTECnames <- popnames[c(1:9,11,14:49)] ## get only M1 hap 1 and non Danubian names
M3names <- levels(pop(RADsub)) ## get names for RAD and micro only


##### M1 hap 1 no STEC ################

M1sep <- seppop(M1)
M1_hap1_noSTEC <- repool(M1sep$CCS, M1sep$FFF, M1sep$CAKE, M1sep$BOK, M1sep$MVW, M1sep$MVWZ, M1sep$NLP, M1sep$MY20, M1sep$VIIKCA, M1sep$FFG, M1sep$GR1, M1sep$TU, M1sep$POLEN, M1sep$MOAT, M1sep$GFP, M1sep$BF, M1sep$OTOM, M1sep$RAIL, M1sep$PRIM, M1sep$RM, M1sep$HOLT, M1sep$UMCA, M1sep$OST, M1sep$LMCA, M1sep$AL, M1sep$EK, M1sep$SD, M1sep$GD, M1sep$STEC, M1sep$STYV, M1sep$KAP, M1sep$SK, M1sep$SA, M1sep$CALK, M1sep$OU, M1sep$EST, M1sep$EST2, M1sep$UKR, M1sep$PRO, M1sep$COP, M1sep$OBY, M1sep$PED, M1sep$TROM, M1sep$WEN, M1sep$GAM)
pop(M1_hap1_noSTEC)
M1_hap1_noSTEC

M1H1coords <- MicroCoords[c(1:9,11,14:32,34:49),] ## subset coordinates file

M1H1gen <- genind2genpop(M1_hap1_noSTEC)


M1H1Dgen <- dist.genpop(M1H1gen,method=2) ##Generates the genetic distance semi-matrix

M1H1Dgeo <- dist(M1H1coords) ## Generates the Euclidian distance semi matrix. Note- vignette says to use "dist(body$other$xy),I dont know why but this wont work for me!

M1H1ibd <- mantel.randtest(M1H1Dgen,M1H1Dgeo)

M1H1ibd  ## gives the output for the test  NOTE that the "observation" is the R value. Sqaure to get R^2 (obviously)

plot(M1H1ibd) ##plots histogram fgor test

par(mfrow = c(1,1))

smoothScatter(M1H1Dgeo, M1H1Dgen, nbin = 250, main = "a) M1 Hap1, no STEC data IBD", xlab = "Dgeo", ylab = "Ggen")
points(M1H1Dgeo, M1H1Dgen, pch = ".", cex = 3)

abline(lm(M1H1Dgen~M1H1Dgeo),  col="darkred",lty=2, lwd = 1.5)  ## Adds regression line
M1H1lm <- summary(lm(M1H1Dgen~M1H1Dgeo))
textxy(1,0.2, paste0("adj. R squared = ",round(M1H1lm$adj.r.squared,3),"***"), cex = 1)



############ M3 RAD & Micro pops only ###################

M1sep <- seppop(M1) ## separate into a separate genind object for each pop

M3 <- repool(M1sep$BF, M1sep$CAKE, M1sep$CALK, M1sep$COP, M1sep$MOAT, M1sep$OBY, M1sep$OU, M1sep$PED, M1sep$POLEN, M1sep$PRO, M1sep$SD, M1sep$SK, M1sep$STYV, M1sep$TROM, M1sep$TU, M1sep$WEN)
## repool only the pops than we want in M3

M3gen <- genind2genpop(M3)


M3Dgen <- dist.genpop(M3gen,method=2) ##Generates the genetic distance semi-matrix

M3Dgeo <- dist(MicnRADcoords) ## Generates the Euclidian distance semi matrix. Note- vignette says to use "dist(body$other$xy),I dont know why but this wont work for me!

M3ibd <- mantel.randtest(M3Dgen,M3Dgeo)

M3ibd  ## gives the output for the test  NOTE that the "observation" is the R value. Sqaure to get R^2 (obviously)

plot(M3ibd) ##plots histogram fgor test



############# M2 Microsat subset ###############


M2 <- read.fstat("/home/dan/Dropbox/PhD/Dan's PhD (Shared)/Data/RAD/RADsample_microsat_genotypes_1.dat") ## Load file - converts to GENIND object.
### this file contains only the RAD individuals! This should match up exactly with the RAD data file - ie: 170 individuals, all from the same pops. 
## note the amount of missing data is quite bad in some pops. 

M2pops <- read.delim("/home/dan/Dropbox/PhD/Dan's PhD (Shared)/Data/RAD/RADsample_microsat_POPS.txt", header = F)
pop(M2) <- M2pops$V1


M2sep <- seppop(M2)
M2sub <- repool(M2sep$BF, M2sep$CAKE, M2sep$"CA-LK", M2sep$COP, M2sep$MOAT, M2sep$OBY, M2sep$OU, M2sep$PED, M2sep$POLEN, M2sep$PRO, M2sep$SD, M2sep$SK, M2sep$STYV, M2sep$TROM, M2sep$TU, M2sep$WEN)

pop(M2sub)

M2subgen <- genind2genpop(M2sub)

M2Dgen <- dist.genpop(M2subgen,method=2) ##Generates the genetic distance semi-matrix

M2Dgeo <- dist(MicnRADcoords) ## Generates the Euclidian distance semi matrix. Note- vignette says to use "dist(body$other$xy),I dont know why but this wont work for me!

M2ibd <- mantel.randtest(M2Dgen,M2Dgeo)

M2ibd  ## gives the output for the test  NOTE that the "observation" is the R value. Sqaure to get R^2 (obviously)
M2_R2 = (M2ibd$obs)^2

plot(M2ibd) ##plots histogram fgor test

####### Density Plots  ##############
par(mfrow = c(2,2))
par(mar = c(2,2,1,2))
par(oma = c(1,1,2,1))

######## M1 ###########

smoothScatter(M1H1Dgeo, M1H1Dgen, nbin = 250, main = "a) M1 Hap1, no STEC data IBD", xlab = "Dgeo", ylab = "Ggen")
points(M1H1Dgeo, M1H1Dgen, pch = ".", cex = 3)

abline(lm(M1H1Dgen~M1H1Dgeo),  col="darkred",lty=2, lwd = 1.5)  ## Adds regression line
M1H1lm <- summary(lm(M1H1Dgen~M1H1Dgeo))
textxy(1,0.2, paste0("adj. R squared = ",round(M1H1lm$adj.r.squared,3),"***"), cex = 1)


######### RAD ##########

smoothScatter(RADDgeo, RADDgen, nbin = 250, main = "b) RADseq data IBD", xlab = "Dgeo", ylab = "Ggen")
points(RADDgeo, RADDgen, pch = ".", cex = 3)

abline(lm(RADDgen~RADDgeo), col="darkred",lty=2, lwd = 1.5)  ## Adds regression line

RADlm <- summary(lm(RADDgen~RADDgeo))
textxy(1,0.34, paste0("adj. R squared = ",round(RADlm$adj.r.squared,3),"***"), cex = 1)
#textxy(1,0.34,expression(paste("R"^"2", "= 0.7185***")), cex = 1)

######## M2 ###########

smoothScatter(M2Dgeo, M2Dgen, nbin = 250, main = "c) M2 data IBD", xlab = "Dgeo", ylab = "Ggen")
points(M2Dgeo, M2Dgen, pch = ".", cex = 3)

abline(lm(M2Dgen~M2Dgeo),  col="darkred",lty=2, lwd = 1.5)  ## Adds regression line
M2lm <- summary(lm(M2Dgen~M2Dgeo))
textxy(1,0.2, paste0("adj. R squared = ",round(M2lm$adj.r.squared,3),"***"), cex = 1)


######## M3 ###########

smoothScatter(M3Dgeo, M3Dgen, nbin = 250, main = "d) M3 data IBD", xlab = "Dgeo", ylab = "Ggen")
points(M3Dgeo, M3Dgen, pch = ".", cex = 3)

abline(lm(M3Dgen~M3Dgeo),  col="darkred",lty=2, lwd = 1.5)  ## Adds regression line
M3lm <- summary(lm(M3Dgen~M3Dgeo))
textxy(1,0.2, paste0("adj. R squared = ",round(M3lm$adj.r.squared,3),"***"), cex = 1)


title(main = "Isolation by distance signal compared between 18,908 SNPs and 13 microsatellites in 149 \n pure crucian carp in Europe", outer = T)


