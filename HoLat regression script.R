library(calibrate)

#### Regressions of Ho and latitude. 

#### All Samples ####

hoarlatlong <- read.csv("~/../Dropbox/PhD/Dans_PhD_Shared/Data/Microsatellites/R stats/Long_lat Ho regressions/hoARlatlong.csv", header = F, sep = '\t') ## read csv

names(hoarlatlong) <- c("pop", "Ho", "AR", "Longitude", "Latitude")

hoarlatlong_noNAs <- hoarlatlong[c(1,3:34,36:39,42:49), ] ## remove pops not used in Ar calculations
hoarlatlong_noNAs ## check


#reg1 <- lm(hoarlatlong$Ho ~ hoarlatlong$Latitude) ## 
#reg2 <- lm(hoarlatlong$Ho ~ hoarlatlong$Longitude) ## 

reg3 <- lm(hoarlatlong_noNAs$AR ~ hoarlatlong_noNAs$Latitude)
reg4 <- lm(hoarlatlong_noNAs$AR ~ hoarlatlong_noNAs$Longitude)

cor1 <- cor.test(hoarlatlong_noNAs$Ho , hoarlatlong_noNAs$AR) ## Ho and AR very highly correlated
plot(hoarlatlong_noNAs$Ho, hoarlatlong_noNAs$AR)

## -------------------------------------------------------------------------------

###### Northern Europe only ############

NEU_hoarlatlong <- hoarlatlong_noNAs[c(1:8,10,13:45),] ## remove lineage 2 populations (Danube)
NEU_hoarlatlong

reg5 <- lm(NEU_hoarlatlong$AR ~ NEU_hoarlatlong$Latitude)
reg6 <- lm(NEU_hoarlatlong$AR ~ NEU_hoarlatlong$Longitude)

# ---------------------------------------------------------------------------------
### RAD vs Microsats Ho correlation ##

setwd("~/Dropbox/PhD/Dans_PhD_Shared/Data/Micro_RAD_comp")
Micro_RAD_hobs<- read.delim("Micro_vs_RAD_Hobs.txt")
Micro_RAD_hobs

cor2 <- cor.test(Micro_RAD_hobs$Micro, Micro_RAD_hobs$RAD)

## ------------------------------------------------------------------------------
#### View results ###

#reg1
#reg2
summary(reg3)
summary(reg4)
summary(reg5)
summary(reg6)

cor1
cor2

## --------------------------------------------------------------------------------
####### PLOTS ####

par(mfrow = c(2,2))
par(oma = c(0,0,3,0))

## reg1
#plot(hoarlatlong$Ho ~ hoarlatlong$Latitude, xlab = "Latitude", ylab = "Observed Heterozygosity")
#abline(reg1) ## R squared = -0.02224
#textxy(62,3, labs = expression(paste("R"^"2"," = -0.007")), cex = 1)

## reg2
#plot(hoarlatlong$Ho ~ hoarlatlong$Longitude, xlab = "Longitude", ylab = "Observed Heterozygosity") ## 
#abline(reg2) ## R squared = 0.27***
#textxy(32,0.4, labs = expression(paste("R"^"2"," = 0.27***")), cex = 1)


## reg3
plot(hoarlatlong$AR ~ hoarlatlong$Latitude, xlab = "Latitude", ylab = "Allelic Richness")
abline(reg3) ###R squared = -0.0.007
textxy(62,2.7, labs = expression(paste("R"^"2"," = -0.007")), cex = 0.7)


## reg4
plot(hoarlatlong$AR ~ hoarlatlong$Longitude, xlab = "Longitude", ylab = "Allelic Richness")
abline(reg4) ## R-squared:  0.2887 ***
textxy(0,2.7, labs = expression(paste("R"^"2"," = 0.2887 ***")), cex = 0.7)

## reg5
#plot(hoarlatlong$Ho ~ hoarlatlong$AR, xlab = "Allelic Richness", ylab = "Observed Heterozygosity")
#abline(reg5) ###R-squared: 0.9012 ***
#textxy(1,0.4, labs = expression(paste("R"^"2"," = 0.9012 ***")), cex = 0.7)

## reg6
plot(NEU_hoarlatlong$AR ~ NEU_hoarlatlong$Latitude, xlab = "Latitude", ylab = "Allelic Richness")
abline(reg6) ## R-squared:  0.2887 ***
textxy(0,2.7, labs = expression(paste("R"^"2"," = -0.0223")), cex = 0.7)

## reg7
plot(NEU_hoarlatlong$AR ~ NEU_hoarlatlong$Longitude, xlab = "Longitude", ylab = "Allelic Richness")
abline(reg7) ## R-squared:  -0.0223 ***
textxy(0,2.7, labs = expression(paste("R"^"2"," = 0.3156 ***")), cex = 0.7)




mtext("Correlations between allelic richness and latitude and longitude", outer = T)



######## Without DANUBE + VOL #########

hoarlatlong









### Subsets

Balticlatlong <- hoarlatlong[c(7:9, 14:17, 26:41, 44:49), ]

NSealatlong <- hoarlatlong[c(1:6,18:25),]

reglatbalt <- lm(Balticlatlong$Ho ~ Balticlatlong$Latitude) ## does linear regression
summary(reglatbalt)

reglatbaltAR <- lm(Balticlatlong$AR ~ Balticlatlong$Latitude) ## does linear regression
summary(reglatbaltAR)

reglongbalt <- lm(Balticlatlong$Ho ~ Balticlatlong$Longitude) ## does linear regression
summary (reglongbalt)

reglongbaltAR <- lm(Balticlatlong$AR ~ Balticlatlong$Longitude) ## does linear regression
summary (reglongbaltAR)

reglatNSea <- lm(NSealatlong$Ho ~ NSealatlong$Lat) ## does linear regression
summary(reglatNSea)

reglatNSeaAR <- lm(NSealatlong$AR ~ NSealatlong$Lat) ## does linear regression
summary(reglatNSeaAR)

reglongNSea <- lm(NSealatlong$Ho ~ NSealatlong$Long) ## does linear regression
summary(reglongNSea)

reglongNSeaAR <- lm(NSealatlong$AR ~ NSealatlong$Long) ## does linear regression
summary(reglongNSeaAR)


