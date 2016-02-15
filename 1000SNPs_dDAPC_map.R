#install.packages("mapdata")
#install.packages("mapplots")


library(adegenet)
library(maps)
library(mapdata)
library(mapplots)
par(mar = c(3,3,3,3)) ## set plot window margins
par(pin = c(4,4))
par(mfrow = c(1,1)) ## set plot window size
mycolsolid = c("black", "red", "darkgreen", "blue") # set colours
mycol <- transp(mycolsolid, alpha = .8)

setwd("~/Dropbox/PhD/Dans_PhD_Shared/Papers/Phylogeography paper/Final_data_files/")

rndm_snps<- read.fstat("Pure_cru_1000_rndm_FSTAT2.dat")
length(indNames(rndm_snps))

pops <- c(rep("SWE2", 10), rep("GBR7", 10), rep("SWE9", 10), rep("SWE8", 10), rep("DEN1", 10), rep("FIN4", 8), rep("DEN4", 5), rep("GBR8", 9), rep("HUN2", 6), rep("GBR4", 9), rep("POL4", 10), rep("DEN2", 8), rep("RUS1", 9), rep("SWE12", 9), rep("POL3", 10), rep("SWE10", 9), rep("FIN3", 10), rep("SWE14", 9), rep("NOR2", 9))

pop(rndm_snps) <- pops

cords <- read.csv("Mantelcoordinatels_newcodes.txt", header = T)  ## load my microsat coordinates file. 

length(cords$population)

RADcords <- cords[c(20,4,3,38,44,18,45,39,46,17,43,31,36,33,34,47,16,51,48),]

grp <- find.clusters(rndm_snps, max.n.clust=50, n.start = 100)

RADdapc2 <- dapc(rndm_snps, n.da=100, n.pca=50, grp$grp)
temp <- optim.a.score(RADdapc2) ## suggests 1 but I chose 12

##So now do the real thing, chosing 12 PCs

RADdapc <- dapc(rndm_snps, grp$grp)


RADclustermembs <- RADdapc$posterior

row.names(RADclustermembs) <- pop(rndm_snps)

write.csv(RADclustermembs, "1000snp_cluster_membs.csv")

RADpies <- read.csv("1000snp_cluster_membs.csv", header = T)
names(RADpies) <- c("Population", "Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5", "Cluster 6")

RADpie_means <- sapply(split(RADpies[2:7], RADpies$Population), colMeans) ## This just makes population cluster averages for each cluster

write.csv(RADpie_means, "1000snp_mean_cluster_membs.csv")

##################### Map ###########################

RADpies <- read.csv("1000snp_mean_cluster_membs.csv")
mycols = rainbow(8)
mycol <- transp(mycols, 0.6)


svg("~/Dropbox/PhD/Dans_PhD_Shared/Papers//Phylogeography paper//MANUSCRIPT//SOM_FIGS/SVGs/1000_snp_map.svg")

map("worldHires", xlim=c(-10, 45), ylim=c(43,72), col="gray90",  fill=TRUE)##plots the area of the map that I am interested in (just Europe, leave out x/ylim arguments for whole world view)
points(RADcords$lon, RADcords$lat, pch = 16, cex = 0.7)


## Pies ## 

add.pie(RADpies$GBR8  ,x=  -2  ,y=  50	,labels="",radius=1,edges=200,clockwise=T,	col	=	mycol)
add.pie(RADpies$GBR4	,x=	-2	,y=	54	,labels="",radius=1,edges=200,clockwise=T,	col	=	mycol)
add.pie(RADpies$FIN3	,x=	25.76	,y=	62.26	,labels="",radius=1,edges=200,clockwise=T,	col	=	mycol)
add.pie(RADpies$DEN1	,x=	12.55	,y=	53	,labels="",radius=1,edges=200,clockwise=T,	col	=	mycol)
add.pie(RADpies$DEN4  ,x=	13	,y=	57	,labels="",radius=1,edges=200,clockwise=T,	col	=	mycol)
add.pie(RADpies$GBR7	,x=	2	,y=	54	,labels="",radius=1,edges=200,clockwise=T,	col	=	mycol)
add.pie(RADpies$SWE12	,x=	17.79	,y=	62.21	,labels="",radius=1,edges=200,clockwise=T,	col	=	mycol)
add.pie(RADpies$FIN4	,x=	25.47	,y=	65.01	,labels="",radius=1,edges=200,clockwise=T,	col	=	mycol)
add.pie(RADpies$DEN2	,x=	8.34	,y=	55.73	,labels="",radius=1,edges=200,clockwise=T,	col	=	mycol)
add.pie(RADpies$POL4	,x=	22.02	,y=	53.83	,labels="",radius=1,edges=200,clockwise=T,	col	=	mycol)
add.pie(RADpies$RUS1	,x=	42	,y=	47.29	,labels="",radius=1,edges=200,clockwise=T,	col	=	mycol)
add.pie(RADpies$SWE2	,x=	14	,y=	61.5	,labels="",radius=1,edges=200,clockwise=T,	col	=	mycol)
add.pie(RADpies$SWE8	,x=	16.15	,y=	55.55	,labels="",radius=1,edges=200,clockwise=T,	col	=	mycol)
add.pie(RADpies$SWE10	,x=	16.27	,y=	58.56	,labels="",radius=1,edges=200,clockwise=T,	col	=	mycol)
add.pie(RADpies$NOR2	,x=	18.95	,y=	69.65	,labels="",radius=1,edges=200,clockwise=T,	col	=	mycol)
add.pie(RADpies$POL3	,x=	19.3	,y=	52.74	,labels="",radius=1,edges=200,clockwise=T,	col	=	mycol)
add.pie(RADpies$SWE14	,x=	13.31	,y=	59.66	,labels="",radius=1,edges=200,clockwise=T,	col	=	mycol)
add.pie(RADpies$HUN2  ,x=	19.17	,y=	46.49	,labels="",radius=1,edges=200,clockwise=T,	col	=	mycol)
add.pie(RADpies$SWE9  ,x=	14.31	,y=	59.66	,labels="",radius=1,edges=200,clockwise=T,	col	=	mycol)
# n = 17 = Ok!

dev.off()

scatter(RADdapc, scree.da=T, posi.da="bottomright", bg = "gray80", pch= 20, solid=.6, cex=3, col = mycol)
table.value(table(pop(rndm_snps), grp$grp), col.lab=paste("inf", 1:10), row.lab=paste("Pop", levels(pop(rndm_snps)))) ## take a look at the cluster membership table



