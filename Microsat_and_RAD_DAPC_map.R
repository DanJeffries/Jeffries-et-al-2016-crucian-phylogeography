library(maps)
library(mapdata)
library(mapplots)
par(mar = c(0,0,0,0)) ## set plot window margins
par(pin = c(4,4))
par(mfrow = c(1,1)) ## set plot window size
mycolsolid = c("black", "red", "darkgreen", "blue") # set colours
mycol <- transp(mycolsolid, alpha = .8)

setwd("~/Dropbox/PhD/Dan's PhD (Shared)/Data/Microsatellites/DAPC/Complete dataset outputs/")

map("worldHires", xlim=c(-10, 55), ylim=c(43,70), col="gray90", fill=TRUE)##plots the area of the map that I am interested in (just Europe, leave out x/ylim arguments for whole world view)


cords <- read.csv("C:/Users/Dan/Dropbox/PhD/Dan's PhD (Shared)/Data/Microsatellites/R stats/Adegenet IBD/Mantelcoordinates.csv", header = T)  ## load my microsat coordinates file. 
points(cords$lon, cords$lat, pch=19, col="red", cex=0.5)  ## Plots my microsat sample locations on the map 

pies <- read.csv("C:/Users/Dan/Dropbox/PhD/Dan's PhD (Shared)/Data/Microsatellites/DAPC/Complete dataset outputs/Whole EU 4 clusters/WholeEU_mean_clustermemberships.csv", header=T) ## have read this in again and specified that there are headers as I was having difficulty assigning headers to the object. This allows me to call populattions using the $ operator as below.

names(pies)

MicrosatPies<- pies[,c(1,3,6,7,9,15,22,23,25,30,33:35,37,41:46,50)]
names(MiscrosatPies)
MicrosatPies

## Plot map ##

map("worldHires", xlim=c(-10, 45), ylim=c(43,70), col="gray90", fill=TRUE)##plots the area of the map that I am interested in


## UK Pies ##

add.pie(pies$HOLT,x=0.4,y=54.5,labels="",radius=0.847,edges=200,clockwise=T, col = mycol)

add.pie(pies$CAKE,x= -1.9,y=53.8,labels="",radius=0.628833334,edges=200,clockwise=T, col = mycol)

add.pie(pies$BF,x=-2,y=52.3,labels="",radius=0.715833334,edges=200,clockwise=T, col = mycol)

add.pie(pies$RM,x=2.3,y=51.,labels="",radius=0.776166667,edges=200,clockwise=T, col = mycol)

add.pie(pies$MOAT,x=4.9,y=53.3 ,labels="",radius=0.719875,edges=200,clockwise=T, col = mycol)



## Baltic Pies ##


add.pie(pies	$	SK	,	x	=	13.152523	,	y	=	55.550972	,	labels	=	""	,	radius	=	0.959416667	,	edges	=	200,	clockwise	=	T, col = mycol)																																			
add.pie(pies$STYV,x=14.271862,y=57.561081,labels="",radius=0.870625,edges=200,clockwise=T, col = mycol)
aadd.pie(pies$SD,x=12.5,y=63,labels="",radius=0.903333334,edges=200,clockwise=T, col = mycol)

add.pie(pies	$	CALK	,	x	=	25.758348	,	y	=	62.262291	,	labels	=	""	,	radius	=	0.713958334	,	edges	=	200,	clockwise	=	T, col = mycol)																																			
add.pie(pies	$	OU	,	x	=	25.472832	,	y	=	65.012375	,	labels	=	""	,	radius	=	1.046208334	,	edges	=	200,	clockwise	=	T, col = mycol)																																			

## The 3 below do not have allelic richness calculated as they bought the number down too low. Have been given a standard radius, be sure to point out in the figure

add.pie(pies	$	KAP	,	x	=	18.785334	,	y	=	57.849045	,	labels	=	""	,	radius	=	0.7	,	edges	=	200,	clockwise	=	T, col = mycol)																																			


## STEC is not included as it was ommitted from DAPC analyses
add.pie(pies	$	STEC	,	x	=	17.804031	,	y	=	59.601791	,	labels	=	""	,	radius	=	0.676208334	,	edges	=	200,	clockwise	=	T, col = mycol)																																			

## Polish Pies ##


add.pie(pies	$	TU	,	x	=	20.5	,	y	=	50.5	,	labels	=	""	,	radius	=	1.477666667	,	edges	=	200,	clockwise	=	T, col = mycol)																																			

add.pie(pies$POLEN,x=25.022095,y=53,labels="",radius=1.134958334,edges=200,clockwise=T, col = mycol)

## Lower Europe Pies ##

add.pie(pies	$	PRO	,	x	=	40.46814	,	y	=	47.457809	,	labels	=	""	,	radius	=	1.279916667	,	edges	=	200,	clockwise	=	T, col = mycol)																																			

## New pies ## NEED TO DO RADIUSES!!

add.pie(pies  $ COP	,	x	=	12.55	,	y	=	55.77	,	labels	=	""	,	radius	=	1.05	,	edges	=	200,	clockwise	=	T, col = mycol)																																			
add.pie(pies  $	OBY	,	x	=	17.79	,	y	=	60.21	,	labels	=	""	,	radius	=	0.76	,	edges	=	200,	clockwise	=	T, col = mycol)																																			
add.pie(pies  $	PED	,	x	=	12.34	,	y	=	55.73	,	labels	=	""	,	radius	=	0.911	,	edges	=	200,	clockwise	=	T, col = mycol)																																			
add.pie(pies  $	TROM	,	x	=	18.95	,	y	=	69.65	,	labels	=	""	,	radius	=	0.593	,	edges	=	200,	clockwise	=	T, col = mycol)																																			
add.pie(pies  $	WEN	,	x	=	18.95	,	y	=	59.66 ,	labels	=	""	,	radius	=	1	,	edges	=	200,	clockwise	=	T, col = mycol)																																			
add.pie(pies  $	GAM	,	x	=	12.5	,	y	=	56	,	labels	=	""	,	radius	=	1.05	,	edges	=	200,	clockwise	=	T, col = mycol)																																			
