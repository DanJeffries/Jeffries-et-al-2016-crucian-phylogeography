library(adegenet)

body <- read.fstat("C:/Users/Dan/Dropbox/PhD/Dan's PhD (Shared)/Data/Microsatellites/Purecruchecked.DAT")


coordinates <- read.csv("C:/Users/Dan/Dropbox/PhD/Dan's PhD (Shared)/Data/Microsatellites/R stats/Adegenet IBD/Mantelcoordinates.csv",header = TRUE, row.names = 1) ##table containing xy coordinates

############ ALL ############
body@other <- (coords) ## Attach coordinates to the GENIND file

Dgen <- dist.genpop(genbod,method=2) ##Generates the genetic distance semi-matrix

Dgeo <- dist(body$other) ## Generates the Euclidian distance semi matrix. Note- vignette says to use "dist(body$other$xy),I dont know why but this wont work for me!

ibd <- mantel.randtest(Dgen,Dgeo)

ibd  ## gives the output for the test  NOTE that the "observation" is the R value. Sqaure to get R^2 (obviously)

plot(ibd) ##plots histogram fgor test





####### hap1 1 with STEC ##########

coordshap1 <- coords[c(1:9,11, 14:32, 34:43),] ## Subset without STEC and Danube

body@other <- (coords) ## Attach coordinates to the GENIND file


hap1noSTEC <- body[c(1:166,176:196,214:582,614:753),] ## Subsets all Hapolotype 1 (not including STEC) populations (N. EU)

genbodhap1 <- genind2genpop(hap1noSTEC) ## Creates the genpop object for the hap1nostec subset


Dgen <- dist.genpop(genbod,method=2) ##Generates the genetic distance semi-matrix

Dgeo <- dist(body$other) ## Generates the Euclidian distance semi matrix. Note- vignette says to use "dist(body$other$xy),I dont know why but this wont work for me!

ibd <- mantel.randtest(Dgen,Dgeo)

ibd  ## gives the output for the test  NOTE that the "observation" is the R value. Sqaure to get R^2 (obviously)

plot(ibd) ##plots histogram fgor test


plot(Dgeo, Dgen)  ## Plots points as a scatter plot

abline(lm(Dgen~Dgeo), col="red",lty=2)  ## Adds regression line

dens <- kde2d(Dgeo,Dgen, n=300, lims=c(-3, 50,-.1,1))  ##Uses a 2-dimensional kernel density estimation (kde2d) to create an object tht displays desinsty of dots

myPal <- colorRampPalette(c("white","blue","gold", "orange", "red"))  ## I think sets the colour palatte

plot(Dgeo, Dgen, pch=20,cex=.5) ## Plots the scatter with smaller dots

image(dens, col=transp(myPal(300),.7), add=TRUE)

abline(lm(Dgen~Dgeo))
