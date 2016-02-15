library(adegenet)

setwd("/media/dan//34D5D1CE642D7E36//2013076_Hanfling_Bernd/Stacks/Stacks_analyses_V3/Pure_crucian_analysis_V2/populations_r_0.7_p_17_m8_singleSNP/")

RAD <- read.genepop('batch_1.gen')

pops <- read.delim('../populations.txt', header = F)
## made this file using cut from the pop_codes.txt file used in populations. 

pop(RAD) <- pops$V1

## Subset to only shared Mircrosats and RAD sammples and get rid of STEC
RADsub <- RAD[c(1:9,15:117,128:155,162:170),]

pop(RADsub)

RADsubgen <- genind2genpop(RADsub)

#### Get coordinates from file ####
RADcoordsall <- read.csv("/home/dan/Dropbox/PhD/Dan's PhD (Shared)/Data/Microsatellites/DAPC/Complete dataset outputs/Microsat_vs_RAD/RAD sample coordinates.csv", header =T)
RADcoordsall

#### Subset coordinates data ######
MicnRADcoords <- RADcoordsall[c(1,3:5,15,17:21,24,25,27:29,32),]
MicnRADcoords


Dgen <- dist.genpop(RADsubgen,method=2) ##Generates the genetic distance semi-matrix

Dgeo <- dist(MicnRADcoords) ## Generates the Euclidian distance semi matrix. Note- vignette says to use "dist(body$other$xy),I dont know why but this wont work for me!

ibd <- mantel.randtest(Dgen,Dgeo)

ibd  ## gives the output for the test  NOTE that the "observation" is the R value. Sqaure to get R^2 (obviously)

plot(ibd) ##plots histogram fgor test

smoothScatter(Dgen, Dgeo)

abline(lm(Dgeo~Dgen), col="red",lty=2)  ## Adds regression line


