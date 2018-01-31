##plot study sites maps
env<-read.csv("/Users/gabzuq/dropbox/BioGeoAmazonia/Analysis/Data/ambi14.csv")
env_spatial<-subset(env, !is.na(env$Lat))
coordinates(env_spatial)=~Long+Lat
proj4string(env_spatial) <- CRS("+proj=longlat +datum=WGS84")

mfrow <- c(1,1) 
par(mfrow=mfrow, mar=c(0.2,0.8,2,0.5), oma=c(2,1,.5,0.5), mgp=c(1.7,0.6,0))
plot(env_spatial, col="limegreen", pch=19,cex=1)
#library(rworldmap)
#newmap <- getMap(resolution = "low")
plot(newmap,add=TRUE)
maps::map.scale(x=-73,y=3,relwidth = 0.11,ratio = F,cex=1)
north.arrow(xb=-73,yb=0, len=0.5, lab="N")
map.axes(cex.axis=1.5)

library(maps)
list()
install.packages("plotKML")
library(plotKML)
plotKML(soilfern)
plot(soilfern)

mapSOILFERN<-mask(soilfern, AMZ)
plot(mapSOILFERN)
writeRaster(mapSOILFERN, "Soil-fern_crop.tif")
