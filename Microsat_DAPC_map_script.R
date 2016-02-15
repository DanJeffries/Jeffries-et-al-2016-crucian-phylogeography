library(maps)
library(mapdata)
library(mapplots)
par(mfrow = c(1,1))
par(mar = c(1,1,1,1)) ## set plot window margins
par(pin = c(5,5)) ## set plot window size
mycolsolid = c("black", "red", "darkgreen", "blue") # set colours
mycol <- transp(mycolsolid, alpha = .8)

setwd("/home/dan//Dropbox/PhD/Dan's PhD (Shared)/Data/Microsatellites/DAPC/Complete dataset outputs/")

map("worldHires", xlim=c(-10, 55), ylim=c(43,70), col="gray90", fill=TRUE)##plots the area of the map that I am interested in (just Europe, leave out x/ylim arguments for whole world view)

cords <- read.csv("/home/dan/Dropbox/PhD/Dan's PhD (Shared)/Data/Microsatellites/R stats/Adegenet IBD/Mantelcoordinates.csv", header = T)  ## load my microsat coordinates file. 

#geocent <- read.csv("/home/dan/Dropbox/PhD/Dan's PhD (Shared)/Data/Microsatellites/R stats/mapplots/geo cluster centre.csv")  ##  Loads my geographic centre data

points(geocent$lat, geocent$lon, pch=19, col=mycol, cex=3)  ## Plots the geocent data on the map.

mtcords <- read.csv("/home/dan/Dropbox/PhD/Dan's PhD (Shared)/Data/Microsatellites/R stats/Adegenet IBD/MtDNAcoordinates.csv", header = T)  ## Loads my mtDNA sample coordinates

points(cords$lon, cords$lat, pch=19, col="red", cex=0.5)  ## Plots my microsat sample locations on the map 

points(mtcords$lat, mtcords$lon, pch=19, col="blue", cex=0.5)  ## Plots my mtDNA Genbank sample locations on the map 

## Before you proceed you need the cluster identity proportions from the dapc analyses. These are located in the @posterior category of the GENIND object outputted from the dapc analysis. so the below script accesses this:




#############################################
############### WHOLE EU ####################
#############################################


body <- read.csv("/home/dan/Dropbox/PhD/Dan's PhD (Shared)/Data/Microsatellites/DAPC/Complete dataset outputs/Whole EU 22 clusters/WholeEU_Clustermemberships.csv", header = T)  ## Load my cluster membership data file.
head(body)

head(body)
## Next create an object containing the averaged cluster memberships for each subset. 


WHOLE_EU_meanclustermembs <- sapply(split(body[2:12],body$X),colMeans) ##Very handy script, "colMeans" calculates the mean of columns in the subgroups produced as a result of performing split(body[2:13],body$population).  I then wrote this to the .csv below.

WHOLE_EU_meanclustermembs  ## check it 

write.csv(WHOLE_EU_meanclustermembs, "/home/dan/Dropbox/PhD/Dan's PhD (Shared)/Data/Microsatellites/DAPC/Complete dataset outputs/Whole EU 22 clusters/WHOLE_EU_mean_clustermemberships.csv")   ## I manually added my population names and cloumn headers to the csv 

##Just add the below script after plotting the map.

pies <- read.csv("/home/dan/Dropbox/PhD/Dan's PhD (Shared)/Data/Microsatellites/DAPC/Complete dataset outputs/Whole EU 22 clusters/WHOLE_EU_mean_clustermemberships.csv", header=T) ## have read this in again and specified that there are headers as I was having difficulty assigning headers to the object. This allows me to call populattions using the $ operator as below.

pies

library(maps)
library(mapdata)
library(mapplots)

## Plot map ##

map("worldHires", xlim=c(-10, 45), ylim=c(43,70), col="gray90", fill=TRUE)##plots the area of the map that I am interested in

## Add sample locations ## 

points(cords$lon, cords$lat, pch=19, col="red", cex=0.5)  ## Plots my sample locations on the map (just for fun). But dont know how to lable, yet!
points(mtcords$lat, mtcords$lon, pch=19, col="blue", cex=0.5)  ## Plots my mtDNA Genbank sample locations on the map (just for fun). But dont know how to lable, yet!


## Now plot the pies, note that these pies have been offset from their actual population locations so as to allow for minimum overlap ##
## Also note that all Ar's have been divided by 2 for plotting so as to reduce the size of the pies, while retaining relative sizes


myowncols = c("red", "goldenrod", "darkgreen", "darkblue", "darkred", "purple", "yellow", "turquoise", "lightgreen", "pink", "black", "gray90" )

#mycols = rainbow(12)
mycol = transp(myowncols, 0.8)

## UK Pies ##

add.pie(pies$HOLT,x=0.4,y=54.5,labels="",radius=0.847,edges=200,clockwise=T, col = mycol)

add.pie(pies$CCS,x=0.1279768,y=50.5,labels="",radius=0.664541667,edges=200,clockwise=T, col = mycol)

add.pie(pies$CAKE,x=-1.9,y=53.8,labels="",radius=0.628833334,edges=200,clockwise=T, col = mycol)

add.pie(pies$BF,x=-2,y=52.3,labels="",radius=0.715833334,edges=200,clockwise=T, col = mycol)

add.pie(pies$RM,x=2.3,y=51.,labels="",radius=0.776166667,edges=200,clockwise=T, col = mycol)

add.pie(pies$OTOM,x=3.8,y=52.1,labels="",radius=0.635333334,edges=200,clockwise=T, col = mycol)

add.pie(pies$RAIL,x=3.1,y=54.5,labels="",radius=0.77325,edges=200,clockwise=T, col = mycol)

add.pie(pies$MOAT,x=4.9,y=53.3,labels="",radius=0.719875,edges=200,clockwise=T, col = mycol)

add.pie(pies$FFF,x=-2.5,y=50.5,labels="",radius=0.711208334,edges=200,clockwise=T, col = mycol)

## Belgian Pies ##

add.pie(pies$BOK,x=7.5,y=50,labels="",radius=0.711208334,edges=200,clockwise=T, col = mycol)

add.pie(pies$MVW,x=2.8,y=49.2,labels="",radius=0.739416667,edges=200,clockwise=T, col = mycol)
add.pie(pies$MVWZ,x=5.2,y=49.6,labels="",radius=0.737166667,edges=200,clockwise=T, col = mycol)


add.pie(pies  $	FFG	,	x	=	7.558594	,	y	=	51.890053	,	labels	=	""	,	radius	=	1.186083334	,	edges	=	200,	clockwise	=	T, col = mycol)																																			

## Danubian Pies ##

add.pie(pies$GEW3,x=8.502993,y=46.5,labels="",radius=1.172458334,edges=200,clockwise=T, col = mycol)

add.pie(pies$GEW8,x=13.02993,y=46.5,labels="",radius=1.069291667,edges=200,clockwise=T, col = mycol)

add.pie(pies$CA.CR,x=17.5,y=47.8764583,labels="",radius=1.205708334,edges=200,clockwise=T, col = mycol)

## Baltic Pies ##


add.pie(pies	$	SK	,	x	=	13.152523	,	y	=	55.550972	,	labels	=	""	,	radius	=	0.959416667	,	edges	=	200,	clockwise	=	T, col = mycol)																																			
add.pie(pies$STYV,x=14.271862,y=57.561081,labels="",radius=0.870625,edges=200,clockwise=T, col = mycol)
add.pie(pies$LMCA,x=9.5,y=59.453506,labels="",radius=1.28,edges=200,clockwise=T, col = mycol)
add.pie(pies$EK,x=13,y=61,labels="",radius=1.012708334,edges=200,clockwise=T, col = mycol)
add.pie(pies$SD,x=12.5,y=63,labels="",radius=0.903333334,edges=200,clockwise=T, col = mycol)
add.pie(pies$GD,x=15.5,y=64.5,labels="",radius=1.195166667,edges=200,clockwise=T, col = mycol)
add.pie(pies	$	UMCA	,	x	=	20.405216	,	y	=	63.712364	,	labels	=	""	,	radius	=	0.945041667	,	edges	=	200,	clockwise	=	T, col = mycol)
add.pie(pies$OST,x=18.380814,y=61.8,labels="",radius=1.171291667,edges=200,clockwise=T, col = mycol)

add.pie(pies	$	AL	,	x	=	19.852339	,	y	=	60.359329	,	labels	=	""	,	radius	=	1.239541667	,	edges	=	200,	clockwise	=	T, col = mycol)																																			
add.pie(pies	$	SA	,	x	=	23.104763	,	y	=	60.371616	,	labels	=	""	,	radius	=	1.102416667	,	edges	=	200,	clockwise	=	T, col = mycol)																																			
add.pie(pies$VIIKCA,x=29,y=60.363937,labels="",radius=1.102666667,edges=200,clockwise=T, col = mycol)
add.pie(pies	$	CALK	,	x	=	25.758348	,	y	=	62.262291	,	labels	=	""	,	radius	=	0.713958334	,	edges	=	200,	clockwise	=	T, col = mycol)																																			
add.pie(pies	$	NLP	,	x	=	29.676218	,	y	=	62.680271	,	labels	=	""	,	radius	=	0.630083334	,	edges	=	200,	clockwise	=	T, col = mycol)																																			
add.pie(pies	$	OU	,	x	=	25.472832	,	y	=	65.012375	,	labels	=	""	,	radius	=	1.046208334	,	edges	=	200,	clockwise	=	T, col = mycol)																																			
add.pie(pies	$	MY20	,	x	=	30.248108	,	y	=	55.903034	,	labels	=	""	,	radius	=	1.1505	,	edges	=	200,	clockwise	=	T, col = mycol)																																			

## The 3 below do not have allelic richness calculated as they bought the number down too low. Have been given a standard radius, be sure to point out in the figure

add.pie(pies	$	KAP	,	x	=	18.785334	,	y	=	57.849045	,	labels	=	""	,	radius	=	0.7	,	edges	=	200,	clockwise	=	T, col = mycol)																																			
add.pie(pies	$	EST	,	x	=	24.3	,	y	=	57.8	,	labels	=	""	,	radius	=	0.7	,	edges	=	200,	clockwise	=	T, col = mycol)																																			
add.pie(pies$EST2,x=29.5,y=58,labels="",radius=0.7,edges=200,clockwise=T, col = mycol)


## STEC is not included as it was ommitted from DAPC analyses
add.pie(pies	$	STEC	,	x	=	17.804031	,	y	=	59.601791	,	labels	=	""	,	radius	=	0.676208334	,	edges	=	200,	clockwise	=	T, col = mycol)																																			

## Polish Pies ##

add.pie(pies$GR2,x=15.357225,y=52.9312583,labels="",radius=1.265833334,edges=200,clockwise=T, col = mycol)

add.pie(pies	$	GR1	,	x	=	19.5	,	y	=	54.897537	,	labels	=	""	,	radius	=	1.34925	,	edges	=	200,	clockwise	=	T, col = mycol)																																			

add.pie(pies	$	TU	,	x	=	20.5	,	y	=	50.5	,	labels	=	""	,	radius	=	1.477666667	,	edges	=	200,	clockwise	=	T, col = mycol)																																			

add.pie(pies$POLEN,x=25.022095,y=53,labels="",radius=1.134958334,edges=200,clockwise=T, col = mycol)

## Lower Europe Pies ##

add.pie(pies	$	UKR	,	x	=	30.52002	,	y	=	52.469398	,	labels	=	""	,	radius	=	1.448208334	,	edges	=	200,	clockwise	=	T, col = mycol)																																			
add.pie(pies	$	PRO	,	x	=	40.46814	,	y	=	47.457809	,	labels	=	""	,	radius	=	1.279916667	,	edges	=	200,	clockwise	=	T, col = mycol)																																			

## New pies ## NEED TO DO RADIUSES!!

add.pie(pies  $ COP	,	x	=	12.55	,	y	=	55.77	,	labels	=	""	,	radius	=	1.05	,	edges	=	200,	clockwise	=	T, col = mycol)																																			
add.pie(pies  $	OBY	,	x	=	17.79	,	y	=	60.21	,	labels	=	""	,	radius	=	0.76	,	edges	=	200,	clockwise	=	T, col = mycol)																																			
add.pie(pies  $	PED	,	x	=	12.34	,	y	=	55.73	,	labels	=	""	,	radius	=	0.911	,	edges	=	200,	clockwise	=	T, col = mycol)																																			
add.pie(pies  $	TROM	,	x	=	18.95	,	y	=	69.65	,	labels	=	""	,	radius	=	0.593	,	edges	=	200,	clockwise	=	T, col = mycol)																																			
add.pie(pies  $	WEN	,	x	=	18.95	,	y	=	59.66 ,	labels	=	""	,	radius	=	1	,	edges	=	200,	clockwise	=	T, col = mycol)																																			
add.pie(pies  $	GAM	,	x	=	12.5	,	y	=	56	,	labels	=	""	,	radius	=	1.05	,	edges	=	200,	clockwise	=	T, col = mycol)																																			





## Finally, need to add a Key ###

keypie1= c(0,1)

add.pie(keypie1,x=-7.5,y=45,labels="",radius=1.5, edges=200,clockwise=T, col = "black")
add.pie(keypie1,x=2,y=45,labels="",radius=.5, edges=200,clockwise=T, col = "black")
add.pie(keypie1,x=-2,y=45,labels="",radius=1., edges=200,clockwise=T, col = "black")













#############################################
############ North EU no STEC ###############
##############################################


## Note, if you are carrying on straight from DAPC analyses you can obviously skip this writing and reading process. 


body <- read.csv("/home/dan/Dropbox/PhD/Dan's PhD (Shared)/Data/Microsatellites/DAPC/Complete dataset outputs/NorthEU only/NEU_noSTEC_Clustermemberships.csv", header = T)  ## Load my cluster membership data file.
head(body)

head(body)
## Next create an object containing the averaged cluster memberships for each subset. 


NEU_noSTEC_meanclustermembs <- sapply(split(body[2:12],body$X),colMeans) ##Very handy script, "colMeans" calculates the mean of columns in the subgroups produced as a result of performing split(body[2:13],body$population).  I then wrote this to the .csv below.

NEU_noSTEC_meanclustermembs  ## check it 

write.csv(NEU_noSTEC_meanclustermembs, "/home/dan/Dropbox/PhD/Dan's PhD (Shared)/Data/Microsatellites/DAPC/Complete dataset outputs/NorthEU only/NEU_noSTEC_mean_clustermemberships.csv")   ## I manually added my population names and cloumn headers to the csv 

##Just add the below script after plotting the map.

pies <- read.csv("/home/dan/Dropbox/PhD/Dan's PhD (Shared)/Data/Microsatellites/DAPC/Complete dataset outputs/NorthEU only/NEU_noSTEC_mean_clustermemberships.csv", header=T) ## have read this in again and specified that there are headers as I was having difficulty assigning headers to the object. This allows me to call populattions using the $ operator as below.

pies


## Plot map ##

map("worldHires", xlim=c(-10, 45), ylim=c(43,70), col="gray90", fill=TRUE)##plots the area of the map that I am interested in

## Add sample locations ## 

points(cords$lon, cords$lat, pch=19, col="red", cex=0.5)  ## Plots my sample locations on the map (just for fun). But dont know how to lable, yet!
points(mtcords$lat, mtcords$lon, pch=19, col="blue", cex=0.5)  ## Plots my mtDNA Genbank sample locations on the map (just for fun). But dont know how to lable, yet!


## Now plot the pies, note that these pies have been offset from their actual population locations so as to allow for minimum overlap ##
## Also note that all Ar's have been divided by 2 for plotting so as to reduce the size of the pies, while retaining relative sizes

mycol = rainbow(12)
myowncols = c("darkgreen", "#CCFF00FF",  "blue", "purple", "darkred",  "red", "#FF6600FF","green","goldenrod", "lightgreen", "turquoise", "brown" )


mycol = transp(myowncols, 0.8)

## UK Pies ##

add.pie(pies$HOLT,x=0.4,y=54.5,labels="",radius=0.847,edges=200,clockwise=T, col = mycol)

add.pie(pies$CCS,x=0.1279768,y=50.5,labels="",radius=0.664541667,edges=200,clockwise=T, col = mycol)

add.pie(pies$CAKE,x=-1.9,y=53.8,labels="",radius=0.628833334,edges=200,clockwise=T, col = mycol)

add.pie(pies$BF,x=-2,y=52.3,labels="",radius=0.715833334,edges=200,clockwise=T, col = mycol)

add.pie(pies$RM,x=2.3,y=51.,labels="",radius=0.776166667,edges=200,clockwise=T, col = mycol)

add.pie(pies$OTOM,x=3.8,y=52.1,labels="",radius=0.635333334,edges=200,clockwise=T, col = mycol)

add.pie(pies$RAIL,x=3.1,y=54.5,labels="",radius=0.77325,edges=200,clockwise=T, col = mycol)

add.pie(pies$MOAT,x=4.9,y=53.3,labels="",radius=0.719875,edges=200,clockwise=T, col = mycol)

add.pie(pies$FFF,x=-2.5,y=50.5,labels="",radius=0.711208334,edges=200,clockwise=T, col = mycol)

## Belgian Pies ##

add.pie(pies$BOK,x=7.5,y=50,labels="",radius=0.711208334,edges=200,clockwise=T, col = mycol)

add.pie(pies$MVW,x=2.8,y=49.2,labels="",radius=0.739416667,edges=200,clockwise=T, col = mycol)
add.pie(pies$MVWZ,x=5.2,y=49.6,labels="",radius=0.737166667,edges=200,clockwise=T, col = mycol)


add.pie(pies	$	FFG	,	x	=	7.558594	,	y	=	51.890053	,	labels	=	""	,	radius	=	1.186083334	,	edges	=	200,	clockwise	=	T, col = mycol)																																			

## Danubian Pies ##

#add.pie(pies$GEW3,x=8.502993,y=46.5,labels="",radius=1.172458334,edges=200,clockwise=T, col = mycol)

#add.pie(pies$GEW8,x=13.02993,y=46.5,labels="",radius=1.069291667,edges=200,clockwise=T, col = mycol)

#add.pie(pies$CACR,x=17.5,y=47.8764583,labels="",radius=1.205708334,edges=200,clockwise=T, col = mycol)

## Baltic Pies ##


add.pie(pies	$	SK	,	x	=	13.152523	,	y	=	55.550972	,	labels	=	""	,	radius	=	0.959416667	,	edges	=	200,	clockwise	=	T, col = mycol)																																			
add.pie(pies$STYV,x=14.271862,y=57.561081,labels="",radius=0.870625,edges=200,clockwise=T, col = mycol)
add.pie(pies$LMCA,x=9.5,y=59.453506,labels="",radius=1.28,edges=200,clockwise=T, col = mycol)
add.pie(pies$EK,x=13,y=61,labels="",radius=1.012708334,edges=200,clockwise=T, col = mycol)
add.pie(pies$SD,x=12.5,y=63,labels="",radius=0.903333334,edges=200,clockwise=T, col = mycol)
add.pie(pies$GD,x=15.5,y=64.5,labels="",radius=1.195166667,edges=200,clockwise=T, col = mycol)
add.pie(pies	$	UMCA	,	x	=	20.405216	,	y	=	63.712364	,	labels	=	""	,	radius	=	0.945041667	,	edges	=	200,	clockwise	=	T, col = mycol)
add.pie(pies$OST,x=18.380814,y=61.8,labels="",radius=1.171291667,edges=200,clockwise=T, col = mycol)

add.pie(pies	$	AL	,	x	=	19.852339	,	y	=	60.359329	,	labels	=	""	,	radius	=	1.239541667	,	edges	=	200,	clockwise	=	T, col = mycol)																																			
add.pie(pies	$	SA	,	x	=	23.104763	,	y	=	60.371616	,	labels	=	""	,	radius	=	1.102416667	,	edges	=	200,	clockwise	=	T, col = mycol)																																			
add.pie(pies$VIIKCA,x=29,y=60.363937,labels="",radius=1.102666667,edges=200,clockwise=T, col = mycol)
add.pie(pies	$	CALK	,	x	=	25.758348	,	y	=	62.262291	,	labels	=	""	,	radius	=	0.713958334	,	edges	=	200,	clockwise	=	T, col = mycol)																																			
add.pie(pies	$	NLP	,	x	=	29.676218	,	y	=	62.680271	,	labels	=	""	,	radius	=	0.630083334	,	edges	=	200,	clockwise	=	T, col = mycol)																																			
add.pie(pies	$	OU	,	x	=	25.472832	,	y	=	65.012375	,	labels	=	""	,	radius	=	1.046208334	,	edges	=	200,	clockwise	=	T, col = mycol)																																			
add.pie(pies	$	MY20	,	x	=	30.248108	,	y	=	55.903034	,	labels	=	""	,	radius	=	1.1505	,	edges	=	200,	clockwise	=	T, col = mycol)																																			

## The 3 below do not have allelic richness calculated as they bought the number down too low. Have been given a standard radius, be sure to point out in the figure

add.pie(pies	$	KAP	,	x	=	18.785334	,	y	=	57.849045	,	labels	=	""	,	radius	=	0.7	,	edges	=	200,	clockwise	=	T, col = mycol)																																			
add.pie(pies	$	EST	,	x	=	24.3	,	y	=	57.8	,	labels	=	""	,	radius	=	0.7	,	edges	=	200,	clockwise	=	T, col = mycol)																																			
add.pie(pies$EST2,x=29.5,y=58,labels="",radius=0.7,edges=200,clockwise=T, col = mycol)


## STEC is not included as it was ommitted from DAPC analyses
#add.pie(pies	$	STEC	,	x	=	17.804031	,	y	=	59.601791	,	labels	=	""	,	radius	=	0.676208334	,	edges	=	200,	clockwise	=	T, col = mycol)																																			

## Polish Pies ##

add.pie(pies$GR2,x=15.357225,y=52.9312583,labels="",radius=1.265833334,edges=200,clockwise=T, col = mycol)

add.pie(pies	$	GR1	,	x	=	19.5	,	y	=	54.897537	,	labels	=	""	,	radius	=	1.34925	,	edges	=	200,	clockwise	=	T, col = mycol)																																			

add.pie(pies	$	TU	,	x	=	20.5	,	y	=	50.5	,	labels	=	""	,	radius	=	1.477666667	,	edges	=	200,	clockwise	=	T, col = mycol)																																			

add.pie(pies$POLEN,x=25.022095,y=53,labels="",radius=1.134958334,edges=200,clockwise=T, col = mycol)

## Lower Europe Pies ##

add.pie(pies	$	UKR	,	x	=	30.52002	,	y	=	52.469398	,	labels	=	""	,	radius	=	1.448208334	,	edges	=	200,	clockwise	=	T, col = mycol)																																			
#add.pie(pies	$	PRO	,	x	=	40.46814	,	y	=	47.457809	,	labels	=	""	,	radius	=	1.279916667	,	edges	=	200,	clockwise	=	T, col = mycol)																																			

## New pies ## NEED TO DO RADIUSES!!

add.pie(pies  $ COP	,	x	=	12.55	,	y	=	55.77	,	labels	=	""	,	radius	=	1.05	,	edges	=	200,	clockwise	=	T, col = mycol)																																			
add.pie(pies  $	OBY	,	x	=	17.79	,	y	=	60.21	,	labels	=	""	,	radius	=	0.76	,	edges	=	200,	clockwise	=	T, col = mycol)																																			
add.pie(pies  $	PED	,	x	=	12.34	,	y	=	55.73	,	labels	=	""	,	radius	=	0.911	,	edges	=	200,	clockwise	=	T, col = mycol)																																			
add.pie(pies  $	TROM	,	x	=	18.95	,	y	=	69.65	,	labels	=	""	,	radius	=	0.593	,	edges	=	200,	clockwise	=	T, col = mycol)																																			
add.pie(pies  $	WEN	,	x	=	18.95	,	y	=	59.66 ,	labels	=	""	,	radius	=	1	,	edges	=	200,	clockwise	=	T, col = mycol)																																			
add.pie(pies  $	GAM	,	x	=	12.5	,	y	=	56	,	labels	=	""	,	radius	=	1.05	,	edges	=	200,	clockwise	=	T, col = mycol)																																			



add.pie(pies$SD,x=12.5,y=63,labels="",radius=0.903333334,edges=200,clockwise=T, col = mycol)
add.pie(pies$GD,x=15.5,y=64.5,labels="",radius=1.195166667,edges=200,clockwise=T, col = mycol)






## Finally, need to add a Key ###

keypie1= c(0,1)

add.pie(keypie1,x=-7.5,y=45,labels="",radius=1.5, edges=200,clockwise=T, col = "black")
add.pie(keypie1,x=2,y=45,labels="",radius=.5, edges=200,clockwise=T, col = "black")
add.pie(keypie1,x=-2,y=45,labels="",radius=1., edges=200,clockwise=T, col = "black")
