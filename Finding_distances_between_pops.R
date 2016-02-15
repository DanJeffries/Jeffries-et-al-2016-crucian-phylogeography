cords <- read.csv("C:/Users/Dan/Dropbox/PhD/Dan's PhD (Shared)/Data/Microsatellites/R stats/Adegenet IBD/Mantelcoordinates.csv", header = T)  ## load my microsat coordinates file. 

RADcoordsall <- read.csv("C:/Users/Dan/Dropbox/PhD/Dan's PhD (Shared)/Data/Microsatellites/DAPC/Complete dataset outputs/Microsat_vs_RAD/RAD sample coordinates.csv", header =T)
MicnRADcoords <- RADcoordsall[c(1,3:5,15,17:21,24,25,27:29,32),]

cords

### Find a way to quantify the average distances between populations in both datasets ###

## fishers exact??
## or just rough average? 
## Euclidean distances?