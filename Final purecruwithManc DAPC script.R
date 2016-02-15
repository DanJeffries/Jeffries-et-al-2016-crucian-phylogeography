## This is all the pure_cru microsat data , inclusnig that from RAD samples analysed in Manchester
## Data and Outputs saved in 'C:\Users\Dan\Dropbox\PhD\Dan's PhD (Shared)\Data\Microsatellites\DAPC\Complete dataset outputs'

library(adegenet) ### Load package

mycol = c("darkgreen", "lightblue", "goldenrod","orange",  "green", "purple", "red", "blue", "lightgreen", "darkred", "turquoise", "brown" )


body <- read.fstat("/home/dan/Dropbox/PhD/Dan's PhD (Shared)/Data/Microsatellites/DAPC/Complete dataset outputs/PurecruplusManc_stripped_and_checked.DAT") ## Load file - converts to GENIND object.


body$loc.names     ## Check loci names

body$pop.names     ## Check population names


par(mar = c(0,0,0,0)) ## set plot window margins
par(pin = c(1,1)) ## set plot window size
par(mfrow = c(1,1))
WHOLE_EU_grp <- find.clusters (body, max.n.clust=100)  ## Finds clusters. The BIC scores are not very informative here...## chose 120 PCs and 4 clusters, which captures the broad structure in whole of EU


table.value(table(pop(body), WHOLE_EU_grp$grp), col.lab=paste("inf", 1:30), row.lab=paste("ori", body$pop.names))  ##Table showing cluster assignment
## Saved in 'C:\Users\Dan\Dropbox\PhD\Dan's PhD (Shared)\Data\Microsatellites\DAPC\Complete dataset outputs\WholeEU_cluster_assignments_table.svg'



WHOLE_EU_dapc2 <- dapc(body, n.da=100, n.pca=50)
temp <- optim.a.score(WHOLE_EU_dapc2) ## These two lines calculate the a-score showing how many pCs to retain in the DAPC analyses to follow
## Saved in 'C:\Users\Dan\Dropbox\PhD\Dan's PhD (Shared)\Data\Microsatellites\DAPC\Complete dataset outputs\WholeEU_Spline_interpolation.svg'



WHOLE_EU_dapc1 <- dapc(body, WHOLE_EU_grp$grp) ## RUN DAPC - Chose 24 PCs (from spline analysis) and chose 3 LDs

head(WHOLE_EU_dapc1$tab)
PC1 <- WHOLE_EU_dapc1$tab$"PCA-pc.1"

PC2 <- WHOLE_EU_dapc1$tab$"PCA-pc.2"

row.names(PC1) <- pop(body)

plot(PC1, PC2)

str(PC1)
scatter(PC1, PC2)
?scatter
scatter(WHOLE_EU_dapc1, scree.da=T, posi.da="topleft", bg="gray80", pch=20, solid=.4, cex=3, col = mycol) ## good combo
## Saved in 'C:\Users\Dan\Dropbox\PhD\Dan's PhD (Shared)\Data\Microsatellites\DAPC\Complete dataset outputs\DAPC scatter.svg'

par(mfrow = c(1,2))
scatter(WHOLE_EU_dapc1, scree.da=T, posi.da="topleft", bg="gray80", pch=20, solid=.4, cex=3) ## good combo
table.value(table(pop(body), WHOLE_EU_grp$grp), col.lab=paste("inf", 1:30), row.lab=paste("ori", body$pop.names))  ##Table showing cluster assignment


par(pin = c(7,4)) ## set different plot boundaries
compoplot(WHOLE_EU_dapc1, posi=list(x=900,y=0.8), txt.leg=paste("Cluster", 1:4), lab="", ncol=1, xlab="individuals", col = mycol)  ## structure-like plot!
## not saved

contrib <- loadingplot(WHOLE_EU_dapc1$var.contr, axis=2, thres=.07, lab.jitter=1, main = "Vertical axis (LD1) loading plot") ##Looks at the contribution of each allele to the variation in the data  (Vert)
## Saved in 'C:\Users\Dan\Dropbox\PhD\Dan's PhD (Shared)\Data\Microsatellites\DAPC\Complete dataset outputs\Vertical axis (LD1) loading plot

contrib <- loadingplot(WHOLE_EU_dapc1$var.contr, axis=1, thres=.07, lab.jitter=1, main = "Horizontal axis (LD2) loading plot") ##  (Horiz) 
## Saved in 'C:\Users\Dan\Dropbox\PhD\Dan's PhD (Shared)\Data\Microsatellites\DAPC\Complete dataset outputs\Horizontal axis (LD2) loading plot'

## Write the cluster membership data to a file for plotting on a map

Clustermemberships <- (WHOLE_EU_dapc1$posterior) ##For mapplots

row.names(Clustermemberships)= pop(body) ## add them to the cluster memberships
head(Clustermemberships) ##check
tail(Clustermemberships)
Clustermemberships

write.csv(Clustermemberships, file = "/home/dan/Dropbox/PhD/Dan's PhD (Shared)/Data/Microsatellites/DAPC/Complete dataset outputs/Whole EU 22 clusters/WholeEU_Clustermemberships.csv")  ## Wrote cluster memberships to csv file.

testcols = rainbow(15)
pie(rep(1,15), col = testcols)

testcols

###########################################################################
################### Northern Europe Without STEC ##########################
###########################################################################

mycols = c("darkgreen", "#CCFF00FF",  "blue", "purple", "darkred",  "red", "orange","green","goldenrod", "lightgreen", "turquoise", "brown" )


mycol = rainbow(20)
setwd("/home/dan/Dropbox/PhD/Dan's PhD (Shared)/Data/Microsatellites/DAPC/Complete dataset outputs/Whole EU 4 clusters")

table(pop(body))


hap1noSTEC <- body[c(1:166, 176:196,214:582,614:743,754:839),] ## Only populations in N.EU (mtDNA haplotype 1 excluding PRO, GEW's and CA-CR)
table(pop(hap1noSTEC)) ##check


hap1noSTECgrp <- find.clusters (hap1noSTEC, max.n.clust=100)  ## Finds clusters, here i chose 120 PC's and 12 clusters according to the BIC scores.
## chose 100 PCs and 3 clusters (the BIC score was ambiguous)!
## Saved in 'C:\Users\Dan\Dropbox\PhD\Dan's PhD (Shared)\Data\Microsatellites\DAPC\Complete dataset outputs\NorthEU only\find_clusters BIC scrores.svg'
## and 'C:\Users\Dan\Dropbox\PhD\Dan's PhD (Shared)\Data\Microsatellites\DAPC\Complete dataset outputs\NorthEU only\find_clusters NEU PC_variance.svg'


table.value(table(pop(hap1noSTEC), hap1noSTECgrp$grp), col.lab=paste("inf", 1:30), row.lab=paste("ori", hap1noSTEC$pop.names))  ##Table showing cluster assignment

hap1noSTECdapc2 <- dapc(hap1noSTEC, n.da=100, n.pca=50)
temp <- optim.a.score(dapc2) ## These two lines calculate the a-score showing how many pCs to retain in the DAPC analyses to follow
## Saved as 'C:\Users\Dan\Dropbox\PhD\Dan's PhD (Shared)\Data\Microsatellites\DAPC\Complete dataset outputs\NorthEU only\Spline.svg'


hap1noSTECdapc1 <- dapc(hap1noSTEC, hap1noSTECgrp$grp) ## RUN DAPC - Chose 24 PCs (from spline analysis) and chose 2 LDs
## Saved as 'C:\Users\Dan\Dropbox\PhD\Dan's PhD (Shared)\Data\Microsatellites\DAPC\Complete dataset outputs\NorthEU only\DAPC NEU PC_variance.svg'
## Saved as 'C:\Users\Dan\Dropbox\PhD\Dan's PhD (Shared)\Data\Microsatellites\DAPC\Complete dataset outputs\NorthEU only\DAPC LD F_statistic'


scatter(hap1noSTECdapc1, scree.da=T, posi.da="topleft", bg="white", pch=20, solid=.6, cex=3, col = mycols) ## good combo
## Saved in' C:\Users\Dan\Dropbox\PhD\Dan's PhD (Shared)\Data\Microsatellites\DAPC\Complete dataset outputs\NorthEU only\NEU_DAPC_scatter.svg'


par(pin = c(7,4)) ## set different plot boundaries
compoplot(hap1noSTECdapc1, posi=list(x=900,y=0.8), txt.leg=paste("Cluster", 1:3), lab="", ncol=1, xlab="individuals", col = c('blue', 'darkgreen', 'red'))  ## structure-like plot!
##not saved

contrib <- loadingplot(hap1noSTECdapc1$var.contr, axis=2, thres=.07, lab.jitter=1, main = 'Loading plot Axis2 (LD2)') ##Looks at the contribution of each allele to the variation in the data Axis 2 (Vert)
## Saved as 'C:\Users\Dan\Dropbox\PhD\Dan's PhD (Shared)\Data\Microsatellites\DAPC\Complete dataset outputs\NorthEU only\Loading plot axis2.svg'

contrib <- loadingplot(hap1noSTECdapc1$var.contr, axis=1, thres=.07, lab.jitter=1) ## Axis 1 (Horiz) 
## Saved as 'C:\Users\Dan\Dropbox\PhD\Dan's PhD (Shared)\Data\Microsatellites\DAPC\Complete dataset outputs\NorthEU only\Loading plot axis1.svg'

NEU_noSTEC_Clustermemberships <- (hap1noSTECdapc1$posterior) ##For mapplots



row.names(NEU_noSTEC_Clustermemberships)= pop(hap1noSTEC) ## add them to the cluster memberships
head(NEU_noSTEC_Clustermemberships, n = 20) ##check


write.csv(NEU_noSTEC_Clustermemberships, file = "/home/dan/Dropbox/PhD/Dan's PhD (Shared)/Data/Microsatellites/DAPC/Complete dataset outputs/NorthEU only/NEU_noSTEC_Clustermemberships.csv")  ## Wrote cluster memberships to csv file.

