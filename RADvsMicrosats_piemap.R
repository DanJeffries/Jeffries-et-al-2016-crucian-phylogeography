## Microsat Vs RAD DAPC
library(adegenet)
library(maps)
library(mapdata)
library(mapplots)
par(mar = c(1,1,1,1)) ## set plot window margins
par(pin = c(7,7)) ## set plot window size


setwd("/media/dan//34D5D1CE642D7E36//2013076_Hanfling_Bernd/Stacks/Stacks_analyses_TRIMMED/Pure_cru_only/populations_r07_p17_m8/")

## made the PLINK input files using vcftools from the batch_1.vcf file outputted by stacks

RAD <- read.PLINK('batch_1_plink_recoded_single_snp.raw',map = "batch_1_single_snp.plink.map", chunkSize=1000, parallel = TRUE, n.cores=3, saveNbAlleles=T)

tail(RAD$loc.names)




##check data

pop(RAD) ## contains individual names instead - change this:

pops <- read.delim('../populations.txt', header = F)
## made this file using cut from the pop_codes.txt file used in populations. 

pop(RAD) <- pops$V1

pop(RAD) ## looks good
nInd(RAD)
indNames(RAD)
other(RAD) ## nothing in here
ploidy(RAD)
head(alleles(RAD)) ## nothing in here

## all looks good!

## Subset for the populations for which I have RAD and 13 Microsats (basically just exclude V)
## Also excluded STEC, as I had to exclude this in microsats as well

RADsub <- RAD[c(1:9,15:117,128:155,162:170),]
indNames(RADsub)


## check missing data ##

plot(RADsub, posi = 'topleft', cex = 0.5) ## can't get pop names on here yet - the yaxis is unsuppressable. 
axis(2, at = seq(1,164,3), labels = rev(pop(RADsub)[seq(1,length(pop(RADsub)), 3)]), las = 2, cex = 0.01) ## complicated axis plotting - need to reverse the populations, cos the function plots down to up. And also need to only plot every 3rd pop name because the axis gets too cramped

## SD has a lot of dropout due to poor DNA quality, but it may still be ok. 

par(mar= c(2,2,3,2))
par(mfrow = c(3,2))
grp <- find.clusters(RADsub, max.n.clust=50, n.start = 100)
## chose 200 (all) PCs


## Do quick test DAPC and perform spline analysis to find optimal number of PCs to retain.
RADdapc2 <- dapc(RADsub, n.da=100, n.pca=50, grp$grp)
temp <- optim.a.score(RADdapc2) ## suggests 1 but I chose 12

##So now do the real thing, chosing 12 PCs

RADdapc <- dapc(RADsub, grp$grp)

RADdapc

par(mfrow = c(1,2))
scatter(RADdapc, scree.da=T, posi.da="topleft", bg = "gray80", pch= 20, solid=.6, cex=3)
table.value(table(pop(RADsub), grp$grp), col.lab=paste("inf", 1:10), row.lab=paste("Pop", levels(pop(RADsub)))) ## take a look at the cluster membership table


contrib <- loadingplot(RADdapc$var.contr, axis=2, thres=.00114292, lab.jitter=1, main = "Vertical axis (LD2) loading plot") ##Looks at the contribution of each allele to the variation in the data  (Vert)
## E to west in N. Europe axis
contrib <- loadingplot(RADdapc$var.contr, axis=1, thres=.07, lab.jitter=1, main = "Horizontal axis (LD1) loading plot") ##  (Horiz) 
## North to south variation axis

summary(contributions[,2])

contributions <- RADdapc$var.contr ## the loadings of each locus in the DAPC for retained DA's
quantile(contributions[,2], .99) ## finding the 99th percentile of these values - i.e. value of "informativeness" above which one percent of loci are left


contrib <- loadingplot(RADdapc$var.contr, axis=2, thres=.00114292, lab.jitter=1, main = "Vertical axis (LD2) loading plot") # makes the contrib variable which contains $var.names - the names of the loci above the set threshold.

highly_informative_loci <- contrib$var.names
## The thing is these are the internal names used by DAPC, not the real loc.names. 

All_loci <- RAD$loc.names


write.csv(All_loci, "/media/dan//34D5D1CE642D7E36//2013076_Hanfling_Bernd/Stacks/Stacks_analyses_TRIMMED/Pure_cru_only/populations_r07_p17_m8/All_loci_DAPC_indx.csv")



write.csv(highly_informative_loci, "/media/dan//34D5D1CE642D7E36//2013076_Hanfling_Bernd/Stacks/Stacks_analyses_TRIMMED/Pure_cru_only/populations_r07_p17_m8/highly_informative_loci_NEUvariation.csv")

length(contributions[,2])
length(highly_informative)


summary(RADdapc$var.contr)



###### Export results #############

RADclustermembs <- RADdapc$posterior

row.names(RADclustermembs) <- pop(RADsub)

write.csv(RADclustermembs, "/home/dan/Dropbox/PhD/Dan's PhD (Shared)/Data/Microsatellites/DAPC/Complete dataset outputs/Microsat_vs_RAD/RAD_allClusters_clusteremembs.csv")

RADclustmembs_read <- read.csv("/home/dan/Dropbox/PhD/Dan's PhD (Shared)/Data/Microsatellites/DAPC/Complete dataset outputs/Microsat_vs_RAD/RAD_allClusters_clusteremembs.csv", header = T)
names(RADclustmembs_read) <- c("Population", "Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5", "Cluster 6", "Cluster 7", "Cluster 8", "Cluster 9")
head(RADclustmembs_read)

RADpies <- sapply(split(RADclustmembs_read[2:10], RADclustmembs_read$Population), colMeans) ## This just makes population cluster averages for each cluster
RADpies


write.csv(RADpies, "/home/dan/Dropbox/PhD/Dan's PhD (Shared)/Data/Microsatellites/DAPC/Complete dataset outputs/Microsat_vs_RAD/RAD_allCluster_pies.csv")

#### Get coordinates from file ####
RADcoordsall <- read.csv("/home/dan/Dropbox/PhD/Dan's PhD (Shared)/Data/Microsatellites/DAPC/Complete dataset outputs/Microsat_vs_RAD/RAD sample coordinates.csv", header =T)
RADcoordsall

#### Subset coordinates data ######
MicnRADcoords <- RADcoordsall[c(1,3:5,15,17:21,24,25,27:29,32),]
MicnRADcoords

length(sort(levels(pop(RADsub)))) ## check lengths match up
length(sort(MicnRADcoords$population))

## 16 populations (without STEC)! Looks good 



################################################################
###################### Microsat DAPC- RAD samples only ###########################
################################################################



microsats <- read.fstat("/home/dan/Dropbox/PhD/Dan's PhD (Shared)/Data/RAD/RADsample_microsat_genotypes_1.dat") ## Load file - converts to GENIND object.
### this file contains only the RAD individuals! This should match up exactly with the RAD data file - ie: 170 individuals, all from the same pops. 
## note the amount of missing data is quite bad in some pops. 
pop(microsats)

pop_microsats <- read.delim("/home/dan/Dropbox/PhD/Dan's PhD (Shared)/Data/RAD/RADsample_microsat_POPS.txt", header = F)
pop_microsats
pop(microsats) <- pop_microsats$V1
pop(microsats)


obj<- seppop(RAD_only_microsats)
names(obj)

MicrosubnoSTEC <- repool(obj$BF, obj$CAKE , obj$TU , obj$POLEN , obj$MOAT, obj$SD, obj$STYV, obj$SK, obj$'CA-LK', obj$OU, obj$PRO, obj$COP, obj$OBY, obj$PED, obj$TROM, obj$WEN)

pop(MicrosubnoSTEC)


par(mfrow = c(3,2))
grp <- find.clusters (MicrosubnoSTEC, max.n.clust=100)  
## kept all PCs
## chose 9 clusters


micronoSTECdapc2 <- dapc(MicrosubnoSTEC, n.da=100, n.pca=50)
temp <- optim.a.score(micronoSTECdapc2) ## These two lines calculate the a-score showing how many pCs to retain in the DAPC analyses to follow
## Suggests 7 PCs

micronoSTECdapc1 <- dapc(MicrosubnoSTEC, grp$grp) ## RUN DAPC - Chose 7 PCs (from spline analysis) and chose 2 LDs

par(mfrow = c(1,2))
scatter(micronoSTECdapc1, scree.da=T, posi.da="topleft", bg="gray80", pch=20, solid=.4, cex=3) ## good combo
table.value(table(pop(MicrosubnoSTEC), grp$grp), col.lab=paste("inf", 1:30), row.lab=paste("ori", MicrosubnoSTEC$pop.names))  ##Table showing cluster assignment


## Export cluster memberships for the map ####

Microclustermembs <- (micronoSTECdapc1$posterior) ##For mapplots


row.names(Microclustermembs) <- pop(MicrosubnoSTEC)

write.csv(Microclustermembs, "/home/dan/Dropbox/PhD/Dan's PhD (Shared)/Data/Microsatellites/DAPC/Complete dataset outputs/Microsat_vs_RAD/Micro_AllCluster_clusteremembs.csv")


Microclustmembs_read <- read.csv("/home/dan/Dropbox/PhD/Dan's PhD (Shared)/Data/Microsatellites/DAPC/Complete dataset outputs/Microsat_vs_RAD/Micro_AllCluster_clusteremembs.csv", header = T)
names(Microclustmembs_read) <- c("Population", "Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5", "Cluster 6", "Cluster 7", "Cluster 8", "Cluster 9")
head(Microclustmembs_read)
Micropies <- sapply(split(Microclustmembs_read[2:10], Microclustmembs_read$Population), colMeans) ## This just makes population cluster averages for each cluster

write.csv(Micropies, "/home/dan/Dropbox/PhD/Dan's PhD (Shared)/Data/Microsatellites/DAPC/Complete dataset outputs/Microsat_vs_RAD/All_cluster_Micropies.csv")

Micropies

## Already have coords read in.
## But check the lengths anyway
length(sort(levels(pop(RADsub)))) ## check
length(sort(MicnRADcoords$population)) ## check 


################################################################
###################### Microsat DAPC- Full RAD populations ###########################
################################################################



microsats <- read.fstat("/home/dan/Dropbox/PhD/Dan's PhD (Shared)/Data/Microsatellites/DAPC/Complete dataset outputs/PurecruplusManc_stripped_and_checked.DAT") ## Load file - converts to GENIND object.
pop(microsats)

### this file contains all the microsat data!

### Subset! ##### (NO STEC)
MicrsatFULL_RAD_POPSub <- microsats[c(342:378, 14:40, 281:334, 537:557,614:642, 654:672, 683:726, 744:825)]
pop(MicrsatFULL_RAD_POPSub)

## Check sample sizes ##
par(mfrow = c(2,2))
length(pop(RADsub))
length(pop(MicrosubnoSTEC))
length(pop(MicrsatFULL_RAD_POPSub))

length(levels(pop(RADsub)))
length(levels(pop(MicrosubnoSTEC)))
length(levels(pop(MicrsatFULL_RAD_POPSub)))



table(sort(pop(RADsub))) ## RAD data
table(sort(pop(Microsub))) ## Microsats for individual RAD samples only
table(sort(pop(MicrsatFULL_RAD_POPSub))) ## Microsats for same pops as RAD data, but all individuals in those pops.


par(mfrow = c(3,2))
grp <- find.clusters (MicrsatFULL_RAD_POPSub, max.n.clust=100, n.start = 100)  
## kept all PCs
## chose 9 clusters


micrsatFULLRADPOPSubdapc2 <- dapc(MicrsatFULL_RAD_POPSub, n.da=100, n.pca=50)
temp <- optim.a.score(micrsatFULLRADPOPSubdapc2) ## These two lines calculate the a-score showing how many pCs to retain in the DAPC analyses to follow
## Suggests 8 PCs

micrsatFULLRADPOPSubdapc1 <- dapc(MicrsatFULL_RAD_POPSub, grp$grp) ## RUN DAPC - Chose 7 PCs (from spline analysis) and chose 2 LDs

mycol = rainbow(9)
mycols = transp(mycol, 0.6)

par(mfrow = c(1,2))
scatter(micrsatFULLRADPOPSubdapc1, scree.da=T, posi.da="bottomleft", bg="gray80", pch=20, solid=.4, cex=3, col = mycol) ## good combo
table.value(table(pop(MicrsatFULL_RAD_POPSub), grp$grp), col.lab=paste("inf", 1:30), row.lab=paste("ori", MicrsatFULL_RAD_POPSub$pop.names))  ##Table showing cluster assignment


#par(pin = c(7,4)) ## set different plot boundaries
#compoplot(microdapc1, posi=list(x=900,y=0.8), txt.leg=paste("Cluster", 1:4), lab="", ncol=1, xlab="individuals", col = mycol)  ## structure-like plot!
## not saved

#contrib <- loadingplot(microdapc1$var.contr, axis=2, thres=.07, lab.jitter=1, main = "Vertical axis (LD1) loading plot") ##Looks at the contribution of each allele to the variation in the data  (Vert)

#contrib <- loadingplot(microdapc1$var.contr, axis=1, thres=.07, lab.jitter=1, main = "Horizontal axis (LD2) loading plot") ##  (Horiz) 




## Export cluster memberships for the map ####

MicroFULLclustermembs <- (micrsatFULLRADPOPSubdapc1$posterior) ##For mapplots


row.names(MicroFULLclustermembs) <- pop(MicrsatFULL_RAD_POPSub)

write.csv(MicroFULLclustermembs, "/home/dan/Dropbox/PhD/Dan's PhD (Shared)/Data/Microsatellites/DAPC/Complete dataset outputs/Microsat_vs_RAD/Micro_RAD_FULL_POP_clusteremembs.csv")


MicroFULLclustmembs_read <- read.csv("/home/dan/Dropbox/PhD/Dan's PhD (Shared)/Data/Microsatellites/DAPC/Complete dataset outputs/Microsat_vs_RAD/Micro_RAD_FULL_POP_clusteremembs.csv", header = T)
names(MicroFULLclustmembs_read) <- c("Population", "Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5", "Cluster 6", "Cluster 7", "Cluster 8", "Cluster 9")
head(MicroFULLclustmembs_read)
MicropiesFULL <- sapply(split(MicroFULLclustmembs_read[2:10], MicroFULLclustmembs_read$Population), colMeans) ## This just makes population cluster averages for each cluster

write.csv(MicropiesFULL, "/home/dan/Dropbox/PhD/Dan's PhD (Shared)/Data/Microsatellites/DAPC/Complete dataset outputs/Microsat_vs_RAD/RAD_FULL_POP_Micropies.csv")

MicropiesFULL



#############################
## Looking at the post probs of individual cluster assignment across datasets ##

### RAD indiv. micros only ##
micronoSTEC_posteriors <- micronoSTECdapc1$posterior
micronoSTEC_posteriors_maxes <- apply(micronoSTEC_posteriors, 1, function(x) max(x))

## All indivs in RAD pops ##
micrsatFULLRADPOPSub_posteriors <- micrsatFULLRADPOPSubdapc1$posterior
micrsatFULLRADPOPSub_maxes <-  apply(micrsatFULLRADPOPSub_posteriors, 1, function(x) max(x))

## RAD data ###
RAD_posteriors <- RADdapc$posterior
RAD_maxes <-  apply(RAD_posteriors, 1, function(x) max(x))

## means ## 

mean(RAD_maxes)
mean(micronoSTEC_posteriors_maxes) 
mean(micrsatFULLRADPOPSub_maxes)

## graphs ##

par(mfrow=c(3,1))
barplot(sort(RAD_maxes))
barplot(sort(micronoSTEC_posteriors_maxes)) ## 
barplot(sort(micrsatFULLRADPOPSub_maxes))




#######################################################
######### FULL COMPARATIVE FIGURE ######################
######################################################
library(maps)
library(mapdata)
library(mapplots)
par(mfrow = c(3,2))

#### Get coordinates from file ####
RADcoordsall <- read.csv("/home/dan/Dropbox/PhD/Dan's PhD (Shared)/Data/Microsatellites/DAPC/Complete dataset outputs/Microsat_vs_RAD/RAD sample coordinates.csv", header =T)
RADcoordsall

#### Subset coordinates data ######
MicnRADcoords <- RADcoordsall[c(1,3:5,15,17:21,24,25,27:29,32),]
MicnRADcoords

######################  RAD ##############################

RADpies <- read.csv("/home/dan/Dropbox/PhD/Dan's PhD (Shared)/Data/Microsatellites/DAPC/Complete dataset outputs/Microsat_vs_RAD/RAD_allCluster_pies.csv")

mycols = rainbow(8)
mycol <- transp(mycols, 0.6)

map("worldHires", xlim=c(-10, 45), ylim=c(43,72), col="gray90",  fill=TRUE)##plots the area of the map that I am interested in (just Europe, leave out x/ylim arguments for whole world view)
points(MicnRADcoords$lon, MicnRADcoords$lat, pch = 16, cex = 0.7)


## Pies ## 

add.pie(RADpies$BF  ,x=  -2  ,y=	50	,labels="",radius=1,edges=200,clockwise=T,	col	=	mycol)
add.pie(RADpies$CAKE	,x=	-2	,y=	54	,labels="",radius=1,edges=200,clockwise=T,	col	=	mycol)
add.pie(RADpies$CA.LK	,x=	25.76	,y=	62.26	,labels="",radius=1,edges=200,clockwise=T,	col	=	mycol)
add.pie(RADpies$COP	,x=	12.55	,y=	53	,labels="",radius=1,edges=200,clockwise=T,	col	=	mycol)
add.pie(RADpies$MOAT	,x=	2	,y=	54	,labels="",radius=1,edges=200,clockwise=T,	col	=	mycol)
add.pie(RADpies$OBY	,x=	17.79	,y=	62.21	,labels="",radius=1,edges=200,clockwise=T,	col	=	mycol)
add.pie(RADpies$OU	,x=	25.47	,y=	65.01	,labels="",radius=1,edges=200,clockwise=T,	col	=	mycol)
add.pie(RADpies$PED	,x=	8.34	,y=	55.73	,labels="",radius=1,edges=200,clockwise=T,	col	=	mycol)
add.pie(RADpies$POLEN	,x=	22.02	,y=	53.83	,labels="",radius=1,edges=200,clockwise=T,	col	=	mycol)
add.pie(RADpies$PRO	,x=	42	,y=	47.29	,labels="",radius=1,edges=200,clockwise=T,	col	=	mycol)
add.pie(RADpies$SD	,x=	14	,y=	61.5	,labels="",radius=1,edges=200,clockwise=T,	col	=	mycol)
add.pie(RADpies$SK	,x=	16.15	,y=	55.55	,labels="",radius=1,edges=200,clockwise=T,	col	=	mycol)
add.pie(RADpies$STYV	,x=	16.27	,y=	58.56	,labels="",radius=1,edges=200,clockwise=T,	col	=	mycol)
add.pie(RADpies$TROM	,x=	18.95	,y=	69.65	,labels="",radius=1,edges=200,clockwise=T,	col	=	mycol)
add.pie(RADpies$TU	,x=	19.3	,y=	52.74	,labels="",radius=1,edges=200,clockwise=T,	col	=	mycol)
add.pie(RADpies$WEN	,x=	13.31	,y=	59.66	,labels="",radius=1,edges=200,clockwise=T,	col	=	mycol)
# n = 16 = Ok!


scatter(RADdapc, scree.da=T, posi.da="topleft", bg = "gray80", pch= 20, solid=.6, cex=3, col = mycol)


################### MICROSATS - ONLY RAD SAMPLES ##################################

Micropies <- read.csv("/home/dan/Dropbox/PhD/Dan's PhD (Shared)/Data/Microsatellites/DAPC/Complete dataset outputs/Microsat_vs_RAD/All_cluster_Micropies.csv", header = T)

map("worldHires", xlim=c(-10, 45), ylim=c(43,72), col="gray90",  fill=TRUE)##plots the area of the map that I am interested in (just Europe, leave out x/ylim arguments for whole world view)
points(MicnRADcoords$lon, MicnRADcoords$lat, pch = 16, cex = 0.7)


##### Set Micro colours #########
mycols = rainbow(12)
mycol = transp(mycols, 0.6)

add.pie(Micropies$BF  ,x=  -2	,y=	50	,labels="",radius=1,edges=200,clockwise=T,	col	=	mycol)
add.pie(Micropies$CAKE	,x=  -2	,y=	54	,labels="",radius=1,edges=200,clockwise=T,	col	=	mycol)
add.pie(Micropies$CA.LK	,x=	25.76	,y=	62.26	,labels="",radius=1,edges=200,clockwise=T,	col	=	mycol)
add.pie(Micropies$COP	,x=  12.55	,y=	53	,labels="",radius=1,edges=200,clockwise=T,	col	=	mycol)
add.pie(Micropies$MOAT	,x=  2	,y=	54	,labels="",radius=1,edges=200,clockwise=T,	col	=	mycol)
add.pie(Micropies$OBY	,x=	17.79	,y=	62.21	,labels="",radius=1,edges=200,clockwise=T,	col	=	mycol)
add.pie(Micropies$OU	,x=	25.47	,y=	65.01	,labels="",radius=1,edges=200,clockwise=T,	col	=	mycol)
add.pie(Micropies$PED	,x=  8.34	,y=	55.73	,labels="",radius=1,edges=200,clockwise=T,	col	=	mycol)
add.pie(Micropies$POLEN	,x=	22.02	,y=	53.83	,labels="",radius=1,edges=200,clockwise=T,	col	=	mycol)
add.pie(Micropies$PRO	,x=	42	,y=	47.29	,labels="",radius=1,edges=200,clockwise=T,	col	=	mycol)
add.pie(Micropies$SD	,x=	14	,y=	61.5	,labels="",radius=1,edges=200,clockwise=T,	col	=	mycol)
add.pie(Micropies$SK	,x=	16.15	,y=	55.55	,labels="",radius=1,edges=200,clockwise=T,	col	=	mycol)
add.pie(Micropies$STYV	,x=	16.27	,y=	58.56	,labels="",radius=1,edges=200,clockwise=T,	col	=	mycol)
add.pie(Micropies$TROM	,x=	18.95	,y=	69.65	,labels="",radius=1,edges=200,clockwise=T,	col	=	mycol)
add.pie(Micropies$TU	,x=	19.3	,y=	52.74	,labels="",radius=1,edges=200,clockwise=T,	col	=	mycol)
add.pie(Micropies$WEN	,x=	13.31	,y=	59.66	,labels="",radius=1,edges=200,clockwise=T,	col	=	mycol)

scatter(micronoSTECdapc1, scree.da=T, posi.da="bottomright", bg="gray80", pch=20, solid=.4, cex=3, col = mycol) ## good combo


################ MICROSATS - FULL RAD POPULATION SAMPLES ##############################

MicropiesFULL <- read.csv("/home/dan/Dropbox/PhD/Dan's PhD (Shared)/Data/Microsatellites/DAPC/Complete dataset outputs/Microsat_vs_RAD/RAD_FULL_POP_Micropies.csv", header = T)

mycols = rainbow(9)
mycol = transp(mycols, 0.6)


map("worldHires", xlim=c(-10, 45), ylim=c(43,72), col="gray90",  fill=TRUE)##plots the area of the map that I am interested in (just Europe, leave out x/ylim arguments for whole world view)
points(MicnRADcoords$lon, MicnRADcoords$lat, pch = 16, cex = 0.7)


add.pie(MicropiesFULL$BF  ,x=  -2  ,y=	50	,labels="",radius=1,edges=200,clockwise=T,	col	=	mycol)
add.pie(MicropiesFULL$CAKE	,x=  -2	,y=	54	,labels="",radius=1,edges=200,clockwise=T,	col	=	mycol)
add.pie(MicropiesFULL$CALK	,x=	25.76	,y=	62.26	,labels="",radius=1,edges=200,clockwise=T,	col	=	mycol)
add.pie(MicropiesFULL$COP	,x=  12.55	,y=	53	,labels="",radius=1,edges=200,clockwise=T,	col	=	mycol)
add.pie(MicropiesFULL$MOAT	,x=  2	,y=	54	,labels="",radius=1,edges=200,clockwise=T,	col	=	mycol)
add.pie(MicropiesFULL$OBY	,x=	17.79	,y=	62.21	,labels="",radius=1,edges=200,clockwise=T,	col	=	mycol)
add.pie(MicropiesFULL$OU	,x=	25.47	,y=	65.01	,labels="",radius=1,edges=200,clockwise=T,	col	=	mycol)
add.pie(MicropiesFULL$PED	,x=  8.34	,y=	55.73	,labels="",radius=1,edges=200,clockwise=T,	col	=	mycol)
add.pie(MicropiesFULL$POLEN	,x=	22.02	,y=	53.83	,labels="",radius=1,edges=200,clockwise=T,	col	=	mycol)
add.pie(MicropiesFULL$PRO	,x=	42	,y=	47.29	,labels="",radius=1,edges=200,clockwise=T,	col	=	mycol)
add.pie(MicropiesFULL$SD	,x=	14	,y=	61.5	,labels="",radius=1,edges=200,clockwise=T,	col	=	mycol)
add.pie(MicropiesFULL$SK	,x=	16.15	,y=	55.55	,labels="",radius=1,edges=200,clockwise=T,	col	=	mycol)
add.pie(MicropiesFULL$STYV	,x=	16.27	,y=	58.56	,labels="",radius=1,edges=200,clockwise=T,	col	=	mycol)
add.pie(MicropiesFULL$TROM	,x=	18.95	,y=	69.65	,labels="",radius=1,edges=200,clockwise=T,	col	=	mycol)
add.pie(MicropiesFULL$TU	,x=	19.3	,y=	52.74	,labels="",radius=1,edges=200,clockwise=T,	col	=	mycol)
add.pie(MicropiesFULL$WEN	,x=	13.31	,y=	59.66	,labels="",radius=1,edges=200,clockwise=T,	col	=	mycol)


scatter(micrsatFULLRADPOPSubdapc1, scree.da=T, posi.da="bottomleft", bg="gray80", pch=20, solid=.4, cex=3, col = mycol) ## good combo
