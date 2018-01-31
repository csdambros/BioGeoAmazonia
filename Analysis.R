# Script created for data analyses of PPBio data
# Created by CSDambros, GZuquim, and GMoulatlet

# Load R packages ####

library(vegan) # install.packages("vegan")
library(ecodist) # install.packages("ecodist")
library(VennDiagram) # install.packages("VennDiagram") # Only required for one graph
library(MuMIn) # install.packages("MuMIn")
library(vegan) # install.packages("vegan")
library(ape) # install.packages("ape")

# Load R functions created specificaly for this paper

source("RFunctions.R") # Script need to be in the same folder as the .RData file
source("https://raw.githubusercontent.com/csdambros/R-functions/master/chaodist.R") # Used if dissimilarities corrected for undersampling are used

# Import data ####

## Species occurrence data in long format ####
occLong<-read.csv2("Data/bioticdata10.csv")

levels(occLong$plotID)

#tapply(occLong$species,occLong$group,function(x)length(unique(x)))

#sort(tapply(occLong$species,paste(occLong$group,occLong$sampl_method),function(x)length(unique(x))))


summary(occLong)

# Change order the groups appear (vertebrates, then inv., then plants)
occLong$group<-factor(occLong$group,levels = c("birds","fish","bats","ant","termite","butterflies","fern","ginger","palm"))

## Define pretty names for groups (only used for plotting)
gName<-c("Birds","Fishes (High Disp)","Fishes (Low Disp)","Bats","Ants (Bait)","Ants (Pitfall)","Ants (Winkler)","Termites","Butterflies","Ferns","Gingers","Palms")

## Environmental data ####

# Import
env<-read.csv("Data/ambi16.csv")

# Set names for rows in env data
rownames(env)<-env$plotID

# Change code in terSteege classification to names
env$class_Steege<-factor(c("Central","Eastern","Northern","Southern","NWA","SWA")[env$class_Steege],levels=c("Central","Eastern","Northern","Southern","NWA","SWA"))

# check if regions are correct in the map
plot(env$Long,env$Lat,pch=21,bg=env$class_Steege)
legend("topright",levels(env$class_Steege),pch=21,pt.bg = 1:nlevels(env$class_Steege))

# Check if the environmental data frame is correct
summary(env)


# Check if there are plots in env that are not in biotic data
env$plotID[!env$plotID%in%levels(occLong$plotID)]

# Check if there are plots in biotic data that are not in env
levels(occLong$plotID)[(!levels(occLong$plotID)%in%env$plotID)]

## Remove from biotic data plots that are not in env or when Latitude is no available
occLong<-occLong[occLong$plotID%in%env$plotID[!is.na(env$Lat)],]

## Remove from env plots without Latitude data
env<-env[!is.na(env$Lat),]

# Recreate factors in biotic and env data frames (remove labels with no occurrence)
occLong$plotID<-factor(occLong$plotID)
env$plotID<-factor(env$plotID)

# Order env data with levels from occurrence data
env<-env[match(levels(occLong$plotID),env$plotID),]
env<-droplevels(env)

# Check if plotID names in env match exactly with plotID names in biotic data (must be zero)
sum(env$plotID!=levels(occLong$plotID))

# Check if plotID levels (as a factor) in env match exactly with plotID names in biotic data (must be zero)
sum(levels(env$plotID)!=levels(occLong$plotID))

# Create column with log of soil bases
env$SumofBases_cmol.log<-log10(subset(env,select="SumofBases_cmol")+0.001)[,1]
summary(env)

### Input missing clay data from other samples (random sample from other plots)
### !Adds noise to the data, potentially reducing the power of statistical tests
env[is.na(env$clay),"clay"]<-sample(env$clay[!is.na(env$clay)],sum(is.na(env$clay)),replace = TRUE)

### Merge observed soil bases data from data estimated using the Kriger method using ferns
env$SumofBases_cmol.log.input<-ifelse(is.na(env$SumofBases_cmol.log),env$KrigeSoil_fernR,env$SumofBases_cmol.log)

### if there are still areas without soil bases, input (as in clay)
env$SumofBases_cmol.log.input[is.na(env$SumofBases_cmol.log.input)]<-sample(env$SumofBases_cmol.log.input[!is.na(env$SumofBases_cmol.log.input)],sum(is.na(env$SumofBases_cmol.log.input)),replace = TRUE)

# Data analysis

## Creating predictor variables

## Create distance matrices ####

### Geographical distance
LongLat<-subset(env,select=c("Long","Lat"))
geoDist<-dist(LongLat)

coordinates(LongLat)<-c("Long","Lat")
proj4string(LongLat)<-CRS("+proj=longlat +datum=WGS84")

UTM<-spTransform(LongLat,CRS("+proj=utm +zone=20 ellps=WGS84"))

geoDistUTM<-dist(UTM@coords)

# Check if distance in UTM = distance in LongLat
#cor(geoDist,geoDistUTM)

### Environmental distance
treeCoverDist<-dist(subset(env,select = "sa.latlong.treecover"))
tempMaxDist<-dist(subset(env,select="CHELSA_bio_5"))
precDryQDist<-dist(subset(env,select="CHELSA_bio_17"))
basesDist<-dist(subset(env,select="SumofBases_cmol"))
basesLogDist<-dist(subset(env,select="SumofBases_cmol.log.input"))
clayDist<-dist(subset(env,select="clay"))
caDist<-dist(subset(env,select="Ca_cmol"))

### Region distance (1 if same region, 0 otherwise)

regionRibasMat<-as.matrix(dist(as.integer(env$class_Ribas)))>0
range(regionRibasMat)
rownames(regionRibasMat)<-labels(clayDist)
regionRibasDist<-as.dist(regionRibasMat)

regionSteegeMat<-as.matrix(dist(as.integer(env$class_Steege)))>0
range(regionSteegeMat)
rownames(regionSteegeMat)<-labels(clayDist)
regionSteegeDist<-as.dist(regionSteegeMat)

regionWWFMat<-as.matrix(dist(as.integer(env$ECO_NAME)))>0
range(regionWWFMat)
rownames(regionWWFMat)<-labels(clayDist)
regionWWFDist<-as.dist(regionWWFMat)

# For fishes only

### Environmental distance
oxDist<-dist(subset(env,select="O2"))
phDist<-dist(subset(env,select="PH"))
igTempDist<-dist(subset(env,select="TEMP"))
condDist<-dist(subset(env,select="COND"))


## Handling biotic data by group from now on

## Separate data into groups (uses group name, sampling method, and dispersa columns)
occLongSplit<-split(occLong,paste(as.integer(occLong$group),occLong$sampl_method,occLong$dispersion))

levels(factor(paste(as.integer(occLong$group),occLong$sampl_method,occLong$dispersion)))

lapply(occLongSplit,head)


## Attribute pretty names for each group
names(occLongSplit)<-gName

# Check names
labels(occLongSplit)

# Check how the data look in each group (is everything correct?)
# For example, are the data identified as winkler trapped ants actually all collected with winkler?
# If not, then there might be a problem with the labels (not the data)
# Almost 1005 sure no problem will happen here, but is always good to check
lapply(occLongSplit,summary)

## Create short table for each group (abundance matrix)
occABSplit<-lapply(occLongSplit,function(x){tapply(x$n,list(x$plotID,x$species),sum)})

## Replace NAs with zeroes
occABSplit<-lapply(occABSplit,function(x)ifelse(is.na(x),0,x))

## Remove sites with zero occurrences (usually not sampled)
occABSplit<-lapply(occABSplit,function(x)x[rowSums(x,na.rm=TRUE)>0,])

## Create Presence/Absence matrices
occPASplit<-lapply(occABSplit,function(x)(x>0)+0)

## Calculate total abundance per plot
NSplit<-lapply(occABSplit,rowSums)

## Calculate species richness per plot
SSplit<-lapply(occPASplit,rowSums)

## Create dissimilarity matrices
ecoDistSplit<-lapply(occPASplit,vegdist,method="jac")
#ecoDistSplit<-lapply(ecoDistSplit,stepacross) # in case you want to use extended dissimilarities
#ecoDistSplit<-lapply(occPASplit,chaodist) # in case you want to use the chao method for jaccard
# The chao method corrects for non detected species, but might require some amount of high quality data (not too many rare species)

## Standardize dissimilarity matrices (sd=1, mean=0)
ecoDistSplitStd<-lapply(ecoDistSplit,decostandDist,na.rm=TRUE)

## Summarize species composition in PCoA axes
pcoaSplit<-lapply(ecoDistSplit,cmdscale)

#########

### Analyses using all the data
# Analyses using only areas that are shared with birds are available below (remove)
# Analyses using environmental data for fish are also available below

### Multiple regression models (all predictors in a single model)

### In a single model
# MRM3 conduct regression using the original predictor variables

# original response matrices
lapply(ecoDistSplit, MRM3,log(geoDist+0.01),treeCoverDist,clayDist,basesLogDist)

# standardized response matrices
lapply(ecoDistSplitStd, MRM3,log(geoDist+0.01),treeCoverDist,clayDist,basesLogDist)


# Standardizing the response and predictor variables
# MRM4 cuts and standardizes all response variables before analyses
AllCombResuStd<-lapply(ecoDistSplitStd, MRM4,log(geoDist+0.01),
                       treeCoverDist,clayDist,basesLogDist)


## Create tables to export with coefficients and p-values
AllCombCoefs<-do.call(cbind,lapply(AllCombResuStd,function(x)x$coef[,1]))[-1,]
AllCombPvals<-do.call(cbind,lapply(AllCombResuStd,function(x)x$coef[,2]))[-1,]

## Set names for rows
rownames(AllCombCoefs)<-c("Distance","TreeCover","Clay","Sum of Bases")
rownames(AllCombPvals)<-c("Distance","TreeCover","Clay","Sum of Bases")

# Export from R to csv spreadsheets
write.csv(round(AllCombCoefs,3),"Results/AllCoefsInMultipleRegression.csv")
write.csv(AllCombPvals,"Results/AllPvalsInMultipleRegression.csv")

### Variance partitioning
varpartSplit<-lapply(ecoDistSplit,varpart3,Y=log(geoDist+0.01),Z=list(treeCoverDist,clayDist,basesLogDist))

# Plot variance partitioning results (better graphs later)
#par(mfrow=c(3,4),mar=c(1,1,1,1))
#lapply(varpartSplit,plot,cex=1.5,Xnames=c("GeoDist","EnvDist"))


### Simple regression models (all variables analyzed individually)

# Create dataframe with all predictors for individual analysis
## List variables
AllVar<-list(Geo=log(geoDist+0.01),Tree=treeCoverDist,Temp=tempMaxDist,Prec=precDryQDist,Clay=clayDist,BasesLog=basesLogDist)

## Standardize predictor variables
AllVarStd<-lapply(AllVar,decostandDist,na.rm=TRUE)

## Run MRM for all variables individually in each group
### !!! It may take some time to run
AllIndResu<-lapply(AllVar,function(x){
  lapply(ecoDistSplitStd,MRM4,predictorDist=x) })

# Organize all coefficients and p-values in a single table
AllIndCoefs<-lapply(AllIndResu,function(x)do.call(cbind,lapply(x,function(y)y$coef[2,1])))
AllIndPvals<-lapply(AllIndResu,function(x)do.call(cbind,lapply(x,function(y)y$coef[2,2])))

AllIndCoefs<-do.call(rbind,AllIndCoefs)
AllIndPvals<-do.call(rbind,AllIndPvals)

# Redefine names for rows
rownames(AllIndCoefs)<-c("Distance","TreeCover","Temp","Prec","Clay","Sum of Bases")
rownames(AllIndPvals)<-c("Distance","TreeCover","Temp","Prec","Clay","Sum of Bases")

# Export tables from R to csv spreadsheets
write.csv(round(AllIndCoefs,3),"Results/AllCoefsInSimpleRegression.csv")
write.csv(AllIndPvals,"Results/AllPvalsInSimpleRegression.csv")


#### TEST RIVERS AS BARRIERS

# Remove the effect of geographical distance from data first

## Create regression model
MRMDist<-lapply(ecoDistSplit, MRM3,predictorDist=log(geoDist+0.01))

sort(unlist(lapply(MRMDist,function(x)x$r.squared[1])))

## Extract Fitted data and residuals
MRMResFit<-lapply(1:length(MRMDist),
                  function(x){
                    dist<-ecoDistSplit[[x]];
                    residuals.mrm(MRMDist[[x]],dist,predictors=adjustPredictor(dist,log(geoDist+0.01)))
                  })

# Extract residuals
MRMRes<-lapply(MRMResFit,function(x)x[[2]])

# Set original names for residual data
names(MRMRes)<-names(ecoDistSplit)

# Run regressions using residuals from distance (can region explain something that space cant?)
lapply(MRMRes,MRM3,predictorDist=regionRibasDist)
lapply(MRMRes,MRM3,predictorDist=regionSteegeDist)
#lapply(MRMRes,MRM3,predictorDist=regionWWFDist)

### Variance partitioning
varpartRibasSplit<-lapply(ecoDistSplit,varpart3,Y=log(geoDist+0.01),Z=list(regionRibasDist))
varpartWWFSplit<-lapply(ecoDistSplit,varpart3,Y=log(geoDist+0.01),Z=list(regionWWFDist))
varpartSteegeSplit<-lapply(ecoDistSplit[1:11],varpart3,Y=log(geoDist+0.01),Z=list(regionSteegeDist))

varpartDivisions<-cbind((unlist(lapply(varpartRibasSplit,function(x)x$part$fract[2,2]))),(unlist(lapply(varpartWWFSplit,function(x)x$part$fract[2,2]))),(unlist(lapply(varpartSteegeSplit,function(x)x$part$fract[2,2]))))

varpartDivisions["Palms",3]<-NA

write.csv(varpartDivisions,file = "Results/varpartDivisions.csv")



# Plot (better graphs later)
#par(mfrow=c(3,4),mar=c(1,1,1,1))
#lapply(1:length(varpartRibasSplit),function(x){
#  plot(varpartRibasSplit[[x]],cex=1.5,Xnames=c("",""))
#  title(main = gName[x])
#}
#)


#### Correlograms ####

correlog<-list()

for(i in 1:length(ecoDistSplit)){
  correlog[[i]]<-mantel.correlog(ecoDistSplit[[i]],adjustPredictor(ecoDistSplit[[i]],geoDistUTM))
}

par(mfrow=c(3,4),mar=c(2,2,.5,.5),oma=c(5,4,0,0))

for(i in 1:12){
  plot(correlog[[i]]$mantel.res[,1]/1000,correlog[[i]]$mantel.res[,3],pch=21,bg=correlog[[i]]$mantel.res[,5]<0.05,
       xlim=c(0,900),ylim=c(-0.2,.7),type="o")
  legend("topright",gName[i],bty = "n")
  abline(h=0,col=2,lty=2)
  #    points(ver$mantel.res[,1],ver$mantel.res[,3],pch=21,bg=ver$mantel.res[,5]<0.05,xlim=c(-5,0),ylim=c(0,.5),type="o")
}
mtext("Correlation",2,outer = TRUE,cex=2.5,line = 1)
mtext("Geographical Distance",1,outer = TRUE,cex=2.5,line=2)

dev.copy2pdf(file="Figures/Correlograms.pdf")

# Using PCoA axes ####

# Put Inambari as reference for contrast
env$class_Ribas<-factor(env$class_Ribas,levels=c("Inambari","Guiana","Napo","Negro","Rondonia", "Tapajos","Tapajos_South"))


options(na.action = "na.fail")

#vars<-c("clay","CHELSA_bio_5","CHELSA_bio_17","SumofBases_cmol.log.input","sa.latlong.treecover","Lat","Long")

vars<-c("clay","SumofBases_cmol.log.input","sa.latlong.treecover","Lat","Long")
pcoa1Coefs<-matrix(NA,length(vars)+length(unique(env$class_Ribas)),length(pcoaSplit))
pcoa1Pvals<-matrix(NA,length(vars)+length(unique(env$class_Ribas)),length(pcoaSplit))
pcoa1R2<-matrix(NA,3,length(pcoaSplit))

R2<-{}

AIC4<-list()
SpatialAutoc<-list()
###
colnames(pcoa1Coefs)<-gName
colnames(pcoa1Coefs)<-gName
colnames(pcoa1R2)<-gName


rownames(pcoa1Coefs)<-c("(Intercept)",paste0("class_Ribas",levels(env$class_Ribas)[-1]),vars)
rownames(pcoa1Pvals)<-c("(Intercept)",paste0("class_Ribas",levels(env$class_Ribas)[-1]),vars)
rownames(pcoa1R2)<-c("Rivers","Space","Env")


for(i in 1:length(pcoaSplit)){

envSel<-env[names(pcoaSplit[[i]][,1]),]
envSel.std<-cbind(decostand(envSel[,vars],"standardize"),class_Ribas=envSel$class_Ribas)
resp<-decostand(pcoaSplit[[i]][,1],"standardize")

#Space<-cbind(envSel$Lat,envSel$Long)

CompleteModel<-lm(resp~class_Ribas+clay+SumofBases_cmol.log.input+
                    sa.latlong.treecover+Lat+Long,data=envSel.std)

plot(resp~envSel.std$Long)

plot(envSel.std$Lat,envSel.std$Long,pch=21,bg=heat.colors(30)[cut(resp,25)])
plot(envSel.std$Lat,envSel.std$Long,pch=21,bg=envSel.std$class_Ribas)


summary(CompleteModel)

AICmodels<-dredge(CompleteModel,extra = list("R^2"),fixed = c("Lat","Long"))
BestModel<-get.models(AICmodels, 1)[[1]]
resuBest<-summary(BestModel)

Rivers<-BestModel$model$class_Ribas
if(is.null(BestModel$model$class_Ribas)){Rivers<-rnorm(length(resp),0,0.001)}

Space<-cbind(BestModel$model$Lat,BestModel$model$Long,rnorm(length(resp),0,0.001))
Env<-cbind(BestModel$model$clay,BestModel$model$SumofBases_cmol.log.input,
           BestModel$model$sa.latlong.treecover,rnorm(length(resp),0,0.001))

R2Part<-varpart(resp,~Rivers,~Space,~Env)

plot(R2Part)
title(gName[[i]])
#form<-paste0("varpart(","resp,~",paste0(colnames(BestModel$model[,-1]),collapse = ",~"),",data=BestModel$model)")
#R2Part<-eval(parse(text=form))

BestModelCoefs<-resuBest$coefficients

pcoa1Coefs[names(BestModel$coefficients),i]<-BestModelCoefs[,1]
pcoa1Pvals[names(BestModel$coefficients),i]<-BestModelCoefs[,4]
pcoa1R2[1:3,i]<-R2Part$part$indfract[1:3,3]

AIC4[[gName[i]]]<-subset(AICmodels,delta<4)
R2[gName[i]]<-AIC4[[i]][1,]$`R^2`
####

#### Test for spatial autocorrelation in residuals

GeoDist<-as.matrix(dist(cbind(envSel$Lat,envSel$Long)))
SpatialAutoc[[gName[i]]]<-Moran.I(BestModel$residuals,GeoDist)
}

# Test which groups show spatial autocorrelation in model residuals
lapply(SpatialAutoc,function(x)x$p.value)<0.05

# Termites and Gingers need further tests 
# To confirm the inflation of type I error rates due to spatial autocorrelation
# are not problematic

#### DEAL WITH SPATIAL AUTOCORRELATION IN GINGERS AND TERMITES

pcoa1Coefs<-cbind(pcoa1Coefs,termitesMEM=NA,gingersMEM=NA)
pcoa1Pvals<-cbind(pcoa1Pvals,termitesMEM=NA,gingersMEM=NA)
pcoa1R2<-cbind(pcoa1R2,termitesMEM=NA,gingersMEM=NA)

R2MEM<-{}
AIC4MEM<-list()
SpatialAutocMEM<-list()

### Gingers - Use significant MEMs as covariates (do not change results)
i=11

envSel<-env[env$plotID%in%names(pcoaSplit[[i]][,1]),]
envSel.std<-cbind(decostand(envSel[,vars],"standardize"),class_Ribas=envSel$class_Ribas)
resp<-decostand(pcoaSplit[[i]][match(envSel$plotID,rownames(pcoaSplit[[i]])),1],"standardize")

# Create geographic distance matrix
GeoDist<-as.matrix(dist(cbind(envSel$Lat,envSel$Long)))

# Run Eigen analysis to generate eigenvectors
E<-eigen(GeoDist)

# Calculate spatial autocorrelation in vectors
MoranPval<-{}
for(m in 1:ncol(E$vectors)){
  MoranPval[m]<-Moran.I(E$vectors[,m],GeoDist)$p.value
  print(m)
}


# What are the associated eigenvalues for those vectors with spatial autocorrelation?
E$values[MoranPval<0.05] # Most negative = fine scale spatial autocorrelation (eg. within grids)

# Run model with significant vectors as covariates
CompleteModel<-lm(resp~class_Ribas+clay+SumofBases_cmol.log.input+
                    sa.latlong.treecover+E$vectors[,MoranPval<0.05],data=envSel.std)

AICmodels<-dredge(CompleteModel,extra = list("R^2"))
BestModel<-get.models(AICmodels, 1)[[1]]
resuBest<-summary(BestModel)

# varpart
Rivers<-BestModel$model$class_Ribas
if(is.null(BestModel$model$class_Ribas)){Rivers<-rnorm(length(resp),0,0.001)}
Space<-cbind(E$vectors[,MoranPval<0.05],rnorm(length(resp),0,0.001))
Env<-cbind(BestModel$model$clay,BestModel$model$SumofBases_cmol.log.input,
           BestModel$model$sa.latlong.treecover,rnorm(length(resp),0,0.001))

R2Part<-varpart(resp,~Rivers,~Space,~Env)
#

BestModelCoefs<-resuBest$coefficients
pcoa1Coefs[,"gingersMEM"]

pcoa1Coefs[na.omit(match(names(BestModel$coefficients),rownames(pcoa1Coefs))),"gingersMEM"]<-BestModelCoefs[names(BestModel$coefficients)%in%rownames(pcoa1Coefs),1]
pcoa1Pvals[na.omit(match(names(BestModel$coefficients),rownames(pcoa1Coefs))),"gingersMEM"]<-BestModelCoefs[names(BestModel$coefficients)%in%rownames(pcoa1Coefs),4]
pcoa1R2[1:3,"gingersMEM"]<-R2Part$part$indfract[1:3,3]

AIC4MEM[["gigersMEM"]]<-subset(AICmodels,delta<4)
R2MEM["gigersMEM"]<-AIC4MEM[["gigersMEM"]][1,]$`R^2`
####

#### Test for spatial autocorrelation in residuals
GeoDist<-as.matrix(dist(cbind(envSel$Lat,envSel$Long)))
SpatialAutocMEM[["gingersMEM"]]<-Moran.I(BestModel$residuals,GeoDist)

####

####
### Termites
i=8

envSel<-env[env$plotID%in%names(pcoaSplit[[i]][,1]),]
envSel.std<-cbind(decostand(envSel[,vars],"standardize"),class_Ribas=envSel$class_Ribas)
resp<-decostand(pcoaSplit[[i]][match(envSel$plotID,rownames(pcoaSplit[[i]])),1],"standardize")

# Create geographic distance matrix
GeoDist<-as.matrix(dist(cbind(envSel$Lat,envSel$Long)))

# Run Eigen analysis to generate eigenvectors
E<-eigen(GeoDist)

# Calculate spatial autocorrelation in vectors
MoranPval<-{}
for(m in 1:ncol(E$vectors)){
  MoranPval[m]<-Moran.I(E$vectors[,m],GeoDist)$p.value
  print(m)
}


# What are the associated eigenvalues for those vectors with spatial autocorrelation?
E$values[MoranPval<0.05] # Most negative = fine scale spatial autocorrelation (eg. within grids)

# Run model with significant vectors as covariates
CompleteModel<-lm(resp~class_Ribas+clay+SumofBases_cmol.log.input+
                    sa.latlong.treecover+E$vectors[,MoranPval<0.05],data=envSel.std)


AICmodels<-dredge(CompleteModel,extra = list("R^2"))
BestModel<-get.models(AICmodels, 1)[[1]]
resuBest<-summary(BestModel)

# varpart
Rivers<-BestModel$model$class_Ribas
if(is.null(BestModel$model$class_Ribas)){Rivers<-rnorm(length(resp),0,0.001)}
Space<-cbind(E$vectors[,MoranPval<0.05],rnorm(length(resp),0,0.001))
Env<-cbind(BestModel$model$clay,BestModel$model$SumofBases_cmol.log.input,
           BestModel$model$sa.latlong.treecover,rnorm(length(resp),0,0.001))

R2Part<-varpart(resp,~Rivers,~Space,~Env)
#


BestModelCoefs<-resuBest$coefficients
pcoa1Coefs[,"termitesMEM"]

pcoa1Coefs[na.omit(match(names(BestModel$coefficients),rownames(pcoa1Coefs))),"termitesMEM"]<-BestModelCoefs[names(BestModel$coefficients)%in%rownames(pcoa1Coefs),1]
pcoa1Pvals[na.omit(match(names(BestModel$coefficients),rownames(pcoa1Coefs))),"termitesMEM"]<-BestModelCoefs[names(BestModel$coefficients)%in%rownames(pcoa1Coefs),4]
pcoa1R2[1:3,"termitesMEM"]<-R2Part$part$indfract[1:3,3]

AIC4MEM[["termitesMEM"]]<-subset(AICmodels,delta<4)
R2MEM["termitesMEM"]<-AIC4MEM[["termitesMEM"]][1,]$`R^2`
####

#### Test for spatial autocorrelation in residuals
GeoDist<-as.matrix(dist(cbind(envSel$Lat,envSel$Long)))
SpatialAutocMEM[["termitesMEM"]]<-Moran.I(BestModel$residuals,GeoDist)


######### End of MEMs ####

### Write PCOA results to csv files

write.csv(round(pcoa1Coefs,3),"Results/PCoACoefs.csv")
write.csv(pcoa1Pvals,"Results/PCoAPvals.csv")

## AIC tables

AIC4g<-lapply(1:length(AIC4),function(x)cbind(Group=names(AIC4)[x],AIC4[[x]]))
AICtable<-do.call(rbind,AIC4g)

write.csv(AICtable,"Results/AICTable.csv")

AIC4gMEM<-lapply(1:length(AIC4MEM),function(x)cbind(Group=names(AIC4MEM)[x],AIC4MEM[[x]]))
AICtableMEM<-do.call(rbind,AIC4gMEM)

write.csv(AICtableMEM,"Results/AICTableMEM.csv")

Importance<-lapply(AIC4,importance)
ImportanceTable<-do.call(rbind,Importance)

write.csv(ImportanceTable,"Results/ImportanceTable.csv")

ImportanceMEM<-lapply(AIC4MEM,importance)
ImportanceTableMEM<-do.call(rbind,ImportanceMEM)

write.csv(ImportanceTableMEM,"Results/ImportanceTableMEM.csv")

write.csv(pcoa1R2,"Results/VarpartPCOA1IndFrac.csv")

#########



#### Only for fishes ####
library(MuMIn)
library(vegan)
library(ape)

options(na.action = "na.fail")

vars<-c("O2","PH","TEMP","COND","Lat","Long")
pcoa1CoefsFish<-matrix(NA,length(vars)+length(unique(env$class_Ribas)),length(pcoaSplit))
pcoa1PvalsFish<-matrix(NA,length(vars)+length(unique(env$class_Ribas)),length(pcoaSplit))
pcoa1R2Fish<-matrix(NA,3,length(pcoaSplit))

R2Fish<-{}

AIC4Fish<-list()
SpatialAutocFish<-list()
###
colnames(pcoa1CoefsFish)<-gName
colnames(pcoa1CoefsFish)<-gName
colnames(pcoa1R2Fish)<-gName


rownames(pcoa1CoefsFish)<-c("(Intercept)",paste0("class_Ribas",levels(env$class_Ribas)[-1]),vars)
rownames(pcoa1PvalsFish)<-c("(Intercept)",paste0("class_Ribas",levels(env$class_Ribas)[-1]),vars)
rownames(pcoa1R2Fish)<-c("Rivers","Space","Env")


for(i in c(2,3)){

  envSel<-env[names(pcoaSplit[[i]][,1]),]
  envSel.std<-cbind(decostand(envSel[,vars],"standardize"),class_Ribas=envSel$class_Ribas)
  resp<-decostand(pcoaSplit[[i]][,1],"standardize")
  
  #Space<-cbind(envSel$Lat,envSel$Long)
  
  CompleteModel<-lm(resp~class_Ribas+O2+
                      PH+TEMP+COND+Lat+Long,data=envSel.std)
  
  
  AICmodels<-dredge(CompleteModel,extra = list("R^2"),fixed = c("Lat","Long"))
  BestModel<-get.models(AICmodels, 1)[[1]]
  resuBest<-summary(BestModel)
  
  Rivers<-BestModel$model$class_Ribas
  if(is.null(BestModel$model$class_Ribas)){Rivers<-rnorm(length(resp),0,0.001)}
  
  Space<-cbind(BestModel$model$Lat,BestModel$model$Long,rnorm(length(resp),0,0.001))
  Env<-cbind(BestModel$model$O2,BestModel$model$PH,
             BestModel$TEMP,BestModel$COND,rnorm(length(resp),0,0.001))
  
  R2Part<-varpart(resp,~Rivers,~Space,~Env)
  #form<-paste0("varpart(","resp,~",paste0(colnames(BestModel$model[,-1]),collapse = ",~"),",data=BestModel$model)")
  #R2Part<-eval(parse(text=form))
  
  BestModelCoefs<-resuBest$coefficients
  
  pcoa1CoefsFish[names(BestModel$coefficients),i]<-BestModelCoefs[,1]
  pcoa1PvalsFish[names(BestModel$coefficients),i]<-BestModelCoefs[,4]
  pcoa1R2Fish[1:3,i]<-R2Part$part$indfract[1:3,3]
  
  AIC4Fish[[gName[i]]]<-subset(AICmodels,delta<4)
  R2Fish[gName[i]]<-AIC4[[i]][1,]$`R^2`
  ####
  
  #### Test for spatial autocorrelation in residuals
  
  GeoDist<-as.matrix(dist(cbind(envSel$Lat,envSel$Long)))
  SpatialAutocFish[[gName[i]]]<-Moran.I(BestModel$residuals,GeoDist)
  
}


## AIC tables

AIC4gFish<-lapply(1:length(AIC4Fish),function(x)cbind(Group=names(AIC4Fish)[x],AIC4Fish[[x]]))
AICtableFish<-do.call(rbind,AIC4gFish)

write.csv(AICtableFish,"Results/AICTableFish.csv")

ImportanceFish<-lapply(AIC4Fish,importance)
ImportanceTableFish<-do.call(rbind,ImportanceFish)

write.csv(ImportanceTableFish,"Results/ImportanceTableFish.csv")

write.csv(pcoa1R2Fish,"Results/VarpartPCOA1IndFracFish.csv")


### Write PCOA results to csv files

write.csv(round(pcoa1CoefsFish,3),"Results/PCoACoefsFish.csv")
write.csv(pcoa1PvalsFish,"Results/PCoAPvalsFish.csv")

##


par(mfrow=c(3,4),mar=c(3,4,.5,0.5))

names(pcoaSplit)

summary(lm(pcoaSplit[[i]][,1]~envSel$clay+
             envSel$SumofBases_cmol.log.input+
             envSel$class_Ribas+
             envSel$sa.latlong.treecover+
             envSel$Lat+envSel$Long))


for(i in 1:12){
  envSel<-env[names(pcoaSplit[[i]][,1]),]
  plot(pcoaSplit[[i]][,1]~envSel$Long,pch=21,bg=envSel$class_Ribas)
  
  #plot(pcoaSplit[[i]][,1]~envSel$class_Ribas)
}

for(i in 1:12){
  plot(pcoaSplit[[i]][,1]~env[names(pcoaSplit[[i]][,1]),]$Lat)
}


########################################
########################################
#### TUKEY TEST COMPARING REGIONS ######

TukeyResu<-list()
combs<-combn(levels(envSel.std$class_Ribas),2)
names<-paste(combs[2,],combs[1,],sep="-")
names2<-paste(combs[1,],combs[2,],sep="-")

TukeyResuMean<-matrix(NA,length(names),length(pcoaSplit))
TukeyResuLow<-TukeyResuMean
TukeyResuHigh<-TukeyResuMean

rownames(TukeyResuMean)<-names

for(i in 1:length(pcoaSplit)){

    envSel<-env[names(pcoaSplit[[i]][,1]),]
  envSel.std<-cbind(decostand(envSel[,vars],"standardize"),class_Ribas=envSel$class_Ribas)
  resp<-decostand(pcoaSplit[[i]][,1],"standardize")
  
  #Space<-cbind(envSel$Lat,envSel$Long)
  
  TukeyResu[[i]]<-TukeyHSD(aov(resp~class_Ribas,data=envSel.std),ordered=TRUE)$class_Ribas

  placeInMat<-rowSums(cbind(match(rownames(TukeyResu[[i]]),names),match(rownames(TukeyResu[[i]]),names2)),na.rm=TRUE)
  
  TukeyResuMean[placeInMat,i]<-TukeyResu[[i]][,1]
  TukeyResuLow[placeInMat,i]<-TukeyResu[[i]][,2]
  TukeyResuHigh[placeInMat,i]<-TukeyResu[[i]][,3]

}

TukeyResu2<-lapply(1:length(pcoaSplit),function(x){
  cbind(gName[x],TukeyResu[[x]])})

TukeyResu3<-do.call(rbind,TukeyResu2)

write.csv(TukeyResu3,file="Results/TukeyTests.csv")



par(mar=c(3,1,1,1),mfrow=c(3,4),oma=c(4,15,0,0))

for(i in 1:ncol(TukeyResuMean)){
plot(NA,xlim=range(c(TukeyResuMean,TukeyResuLow,TukeyResuHigh),na.rm=TRUE),ylim=c(1,length(names)),axes=FALSE)

points(TukeyResuMean[,i],1:length(names))
arrows(TukeyResuLow[,i],1:length(names),TukeyResuHigh[,i],1:length(names),angle = 90,code = 3,length = 0.04)
abline(v=0,lty=2)
axis(1)
legend("topright",gName[i],border = FALSE,bty = "n")
box()

if(i%in%c(1,5,9))axis(2,at = 1:length(names),labels = names,las=2)

}



# GRAPHS ####
## Similarity vs Geographical distance ####

##### Define colors and line types for each group
groupColors<-c(1,"grey30","grey50","gold3",1,"red2","red3","dark blue","purple","dark green","dark blue","red","gold")
groupLty<-c(1,1,1,1,2,2,2,2,4,4,4,4,4,4)
groupLwd<-c(2,2,2,2,2,2,2,2,2,2,2,2,2,2)


# Plot frame
par(mar=c(5,5,3,2),mfrow=c(1,1))
plot(NA,xlim=c(0.01,16),ylim=c(0.01,.5),xlab="Distance (degree)",ylab="Similarity (Jaccard)",cex.lab=2.4)


plot(NA,xlim=c(0.01*111,16*111),ylim=c(0.01,.5),xlab="Distance (Km)",ylab="Similarity (Jaccard)",cex.lab=2.4)


### Plot lines
lapply(1:length(ecoDistSplit),function(x){
  pred<-adjustPredictor(ecoDistSplit[[x]],geoDist*111)
  glm1<-glm(1-ecoDistSplit[[x]]~pred,family = binomial(link="log"),start=c(-0.08,-0.009))
  points(seq(min(pred,na.rm=TRUE),max(pred,na.rm=TRUE),0.05),
         predict(glm1,newdata = list(pred=seq(min(pred,na.rm=TRUE),max(pred,na.rm=TRUE),0.05)),type="response"),
         type="l",col=groupColors[x],lty = groupLty[x],lwd=groupLwd[x]); return(glm1)
})


## Add legend
legend("topright",names(ecoDistSplit),col=groupColors,lty = groupLty,cex = 1,lwd=groupLwd,bty="n",title = "Groups")


##### Similarity vs Environmental distance
# Plot frame
plot(NA,xlim=c(0,4),ylim=c(0,.4),xlab="Dif in Bases (mm)",ylab="Similarity (Jaccard)")

# add lines
lapply(1:length(ecoDistSplit),function(x){
  pred<-adjustPredictor(ecoDistSplit[[x]],basesLogDist)
  glm1<-glm(1-ecoDistSplit[[x]]~pred,family = binomial(link="log"),start=c(-0.08,-0.009))
  points(seq(min(pred,na.rm=TRUE),max(pred,na.rm=TRUE),0.05),predict(glm1,newdata = list(pred=seq(min(pred,na.rm=TRUE),max(pred,na.rm=TRUE),0.05)),type="response"),type="l",col=x,lwd=2,lty=x)
})

## Add legend
legend("topright",names(ecoDistSplit),col=1:length(ecoDistSplit),lty = 1:length(ecoDistSplit),cex=0.7,bty="n")


## Group by group
par(mfrow=c(3,4),mar=c(2,2,.5,.5),oma=c(5,4,0,0))

for(x in 1:12){
  x=3
    plot(NA,xlim=c(0.01*111,16*111),ylim=c(0.01,.5),xlab="Distance (Km)",ylab="Similarity (Jaccard)",cex.lab=2.4)
  pred<-adjustPredictor(ecoDistSplit[[x]],geoDist*111)
  glm1<-glm(1-ecoDistSplit[[x]]~pred,family = binomial(link="log"),start=c(-0.08,-0.009))
  points(seq(min(pred,na.rm=TRUE),max(pred,na.rm=TRUE),0.05),
         predict(glm1,newdata = list(pred=seq(min(pred,na.rm=TRUE),max(pred,na.rm=TRUE),0.05)),type="response"),
         type="l",lwd=groupLwd[x])

  legend("topright",gName[x],bty = "n")
}


mtext("Similarity (Jaccard)",2,outer = TRUE,cex=2.5,line = 1)
mtext("Geographical Distance (Km)",1,outer = TRUE,cex=2.5,line=2)

dev.copy2pdf(file="Figures/Correlograms.pdf")












### Variance partitioning

# Environment vs distance

# Increase plot margins to fit graphs and labels
par(mfrow=c(1,1),mar=c(1,1,3,1))

# Plot and export diagrams for all groups
lapply(1:length(varpartSplit),function(x){
  plot.new()
  resu<-round(varpartSplit[[x]]$part$indfract[,3],2)[c(1,2,3)]
  venn.plot <- draw.pairwise.venn(
    area1 = sum(resu[c(1,2)]),
    area2 = sum(resu[c(2,3)]),
    cross.area = sum(resu[c(2)]),
    category = c("", ""),
    fill = c("black", "green"),
    lty = "blank",
    cex = 4,
    cat.cex = 2,
    cat.pos = c(285, 105),
    cat.dist = 0.09,
    cat.just = list(c(-1, -1), c(1, 1)),
    ext.pos = 30,
    ext.dist = -0.05,
    ext.length = 0.85,
    ext.line.lwd = 2,
    ext.line.lty = 1,sigdigs = 2
  );
  title(main=labels(occLongSplit)[x],cex.main=4)
  dev.copy2pdf(file=paste0("Figures/",labels(occLongSplit)[x],".pdf"))
})


## Biogeographical region vs distance


# Increase plot margins to fit graphs and labels
par(mfrow=c(1,1),mar=c(1,1,3,1))

# Plot and export diagrams for all groups
lapply(1:length(varpartRibasSplit),function(x){
  plot.new()
  resu<-round(varpartRibasSplit[[x]]$part$indfract[,3],2)[c(1,2,3)]
  venn.plot <- draw.pairwise.venn(
    area1 = sum(resu[c(1,2)]),
    area2 = sum(resu[c(2,3)]),
    cross.area = sum(resu[c(2)]),
    category = c("", ""),
    fill = c("black", "blue"),
    lty = "blank",
    cex = 4,
    cat.cex = 2,
    cat.pos = c(285, 105),
    cat.dist = 0.09,
    cat.just = list(c(-1, -1), c(1, 1)),
    ext.pos = 30,
    ext.dist = -0.05,
    ext.length = 0.85,
    ext.line.lwd = 2,
    ext.line.lty = 1,sigdigs = 2
  );
  title(main=labels(occLongSplit)[x],cex.main=4)
  dev.copy2pdf(file=paste0("Figures/",labels(occLongSplit)[x],"RibasxGeo",".pdf"))
})


### Ordination plots comparing biogeographical regions

## Ribas

# Define colors
colors<-rainbow(8)
colors[4]<-"black"  



# Define plot frame
par(mfrow=c(3,4),mar=c(1,1,2,1))

# Plot graphs
lapply(1:length(pcoaSplit),function(y){
  x<-pcoaSplit[[y]]
  groups<-env[rownames(x),]$class_Ribas
  ordiplot(x,type="n",axes=FALSE)
  ordispider(x,groups,draw = "polygon",col=colors,pch=21,bg=1)
  points(x,pch=21,bg=colors[as.integer(groups)],col=1)
  title(labels(pcoaSplit)[y],cex.main=1.5)
  })


dev.copy2pdf(file=paste0("Figures/","OrdiRibas",".pdf"))

par(mfrow=c(1,1),mar=c(1,1,2,1))

# Add legend
plot.new()
legend("top",levels(env$class_Ribas),fill = colors,cex=1.5,bty = "n")
title("Biogeographic region",cex.main=1.5)

dev.copy2pdf(file=paste0("Figures/","OrdiRibasLegend",".pdf"))

## terSteege

# Define colors
colors<-rainbow(8)
colors[4]<-"black"  
colors[c(3,2)]<-colors[c(2,3)]

# Define plot frame
par(mfrow=c(3,4),mar=c(1,1,2,1))

# Add plots
lapply(1:length(pcoaSplit),function(y){
  x<-pcoaSplit[[y]]
  groups<-env[rownames(x),]$class_Steege
  ordiplot(x,type="n",axes=FALSE)
  ordispider(x,groups,draw = "polygon",col=colors,pch=21,bg=1)
  points(x,pch=21,bg=colors[as.integer(groups)],col=1)
  title(labels(pcoaSplit)[y],cex.main=1.5)
})

dev.copy2pdf(file=paste0("Figures/","OrdiSteege",".pdf"))


par(mfrow=c(1,1),mar=c(1,1,2,1))

# Add legend
plot.new()
legend("top",levels(env$class_Steege),fill = colors,cex=1.5,bty = "n")
title("Biogeographic region",cex.main=1.5)

dev.copy2pdf(file=paste0("Figures/","OrdiSteegeLegend",".pdf"))


## WWF

# Define colors
colors<-rainbow(17)
colors[4]<-"black"  

# Define plot frame
par(mfrow=c(3,4),mar=c(1,1,2,1))

# Add plots
lapply(1:length(pcoaSplit),function(y){
  x<-pcoaSplit[[y]]
  groups<-env[rownames(x),]$ECO_NAME
  ordiplot(x,type="n",axes=FALSE)
  ordispider(x,groups,col=colors,pch=21,bg=1)
  points(x,pch=21,bg=colors[as.integer(groups)],col=1)
  title(labels(pcoaSplit)[y],cex.main=1.5)
})

dev.copy2pdf(file=paste0("Figures/","OrdiWWF",".pdf"))


par(mfrow=c(1,1),mar=c(1,1,2,1))

# Add legend
plot.new()
legend("top",levels(env$ECO_NAME),fill = colors,cex=.9,bty = "n")
title("Biogeographic region (WWF)",cex.main=1.2)

dev.copy2pdf(file=paste0("Figures/","OrdiWWFLegend",".pdf"))

#############################
# BIRD  BIRD  BIRD  BIRD  BIRD

# Run analyses only for those sites where birds were sampled.

# Make subset of the data for all groups based on data from birds
ecoDistSplitBird<-lapply(ecoDistSplit, function(x)adjustPredictor(ecoDistSplit$Birds,x))

# Remove groups with no data
ecoDistSplitBird$`Fishes (High Disp)`<-NULL # Set null for fish (no shared sites with birds)
ecoDistSplitBird$`Fishes (Low Disp)`<-NULL # Set null for fish (no shared sites with birds)

# For each group, remove plots with no occurrences (is this necessary?)
ecoDistSplitBird<-lapply(ecoDistSplitBird,function(x){
as.dist(as.matrix(x)[colSums(!is.na(as.matrix(x)))>1,
                                 colSums(!is.na(as.matrix(x)))>1])
})

# Restandardize data using only plots with subset of plots
ecoDistSplitBirdStd<-lapply(ecoDistSplitBird, decostandDist,na.rm=TRUE)

# Summarize data into PCoA axes
pcoaSplitBird<-lapply(ecoDistSplitBird,cmdscale)


## Multiple regression model (all predictors in a single model)
AllCombResuBirdStd<-lapply(ecoDistSplitBirdStd, MRM4,
                           log(geoDist+0.01),
                           treeCoverDist,
                           clayDist,
                           basesLogDist)

# Create table with all coefficients and p-values
AllCombCoefsBird<-do.call(cbind,lapply(AllCombResuBirdStd,function(x)x$coef[,1]))[-1,]
AllCombPvalsBird<-do.call(cbind,lapply(AllCombResuBirdStd,function(x)x$coef[,2]))[-1,]

# Put variable names in rows
rownames(AllCombCoefsBird)<-c("Distance","TreeCover","Clay","Sum of Bases")
rownames(AllCombPvalsBird)<-c("Distance","TreeCover","Clay","Sum of Bases")

# Export R tables to .csv file
write.csv(round(AllCombCoefs,3),"Results/AllCoefsInMultipleRegressionBird.csv")
write.csv(AllCombPvals,"Results/AllPvalsInMultipleRegressionBird.csv")

## Simple regression model
AllVar<-list(Geo=log(geoDist+0.01),Tree=treeCoverDist,Temp=tempMaxDist,Prec=precDryQDist,Clay=clayDist,BasesLog=basesLogDist)

## Standardize predictor variables
AllVarStd<-lapply(AllVar,decostandDist,na.rm=TRUE)

## Run MRM for all variables individually in each group
### !!! It may take some time to run
AllIndResuBird<-lapply(AllVar,function(x){
  lapply(ecoDistSplitBirdStd,MRM4,predictorDist=x) })

# Coefs
AllIndCoefsBird<-lapply(AllIndResuBird,function(x)do.call(cbind,lapply(x,function(y)y$coef[2,1])))

# Pvals
AllIndPvalsBird<-lapply(AllIndResuBird,function(x)do.call(cbind,lapply(x,function(y)y$coef[2,2])))

AllIndCoefsBird<-do.call(rbind,AllIndCoefsBird)
AllIndPvalsBird<-do.call(rbind,AllIndPvalsBird)

rownames(AllIndCoefsBird)<-c("Distance","TreeCover","Temp","Prec","Clay","Sum of Bases")
rownames(AllIndPvalsBird)<-c("Distance","TreeCover","Temp","Prec","Clay","Sum of Bases")

write.csv(round(AllIndCoefsBird,3),"Results/AllCoefsInSimpleRegressionBird.csv")
write.csv(AllIndPvalsBird,"Results/AllPvalsInSimpleRegressionBird.csv")

##########
# GRAPHS

colors<-rainbow(8)
colors[4]<-"black"  

par(mfrow=c(3,4),mar=c(1,1,2,1))

#lapply(pcoaSplit,plot,col=rainbow(70)[as.integer(cut(env$ECO_NUM[match],40))])
lapply(1:length(pcoaSplitBird),function(y){
  x<-pcoaSplitBird[[y]]
  groups<-env[rownames(x),]$class_Ribas
  ordiplot(x,type="n",axes=FALSE)
  #ordispider(x,groups,draw = "polygon",col=colors,pch=21,bg=1)
  ordispider(x,groups,draw = "polygon",col=colors,pch=21,bg=1)
  points(x,pch=21,bg=colors[as.integer(groups)],col=1)
  #title(gName[y],cex.main=1.5)
  title(labels(pcoaSplitBird)[y],cex.main=1.5)
})

plot.new()
legend("top",levels(env$class_Ribas)[levels(env$class_Ribas)%in%env[rownames(pcoaSplitBird$Birds),]$class_Ribas],fill = colors,cex=1.5,bty = "n")
title("Biogeographic region",cex.main=1.5)


###

plot(NA,xlim=c(0,5),ylim=c(0,.9),xlab="Distance (degree)",ylab="Similarity (Jaccard)")

lapply(1:length(ecoDistSplitBird),function(x){
  pred<-adjustPredictor(ecoDistSplitBird[[x]],geoDist)
  glm1<-glm(1-ecoDistSplitBird[[x]]~pred,family = binomial(link="log"))
  points(seq(min(pred,na.rm=TRUE),max(pred,na.rm=TRUE),0.05),predict(glm1,newdata = list(pred=seq(min(pred,na.rm=TRUE),max(pred,na.rm=TRUE),0.05)),type="response"),type="l",col=x,lwd=2,lty=x)
})

legend("topright",names(ecoDistSplitBird),col=1:length(ecoDistSplit),lty = 1:length(ecoDistSplitBird))

####species area curve
occAcccurve<-lapply(occPASplit, specaccum) 
par(mfrow=c(3,4),mar=c(2,2,.5,.5),oma=c(5,4,0,0))

lapply(names(occAcccurve), function(x) {plot(occAcccurve[[x]],xlab="",ylab="",cex.lab=1.5,ci.type="poly",ci.col="grey70");legend("bottomright",x,bty = "n")})

mtext("Number of species",2,outer = TRUE,cex=2.5,line = 1)
mtext("Number of plots",1,outer = TRUE,cex=2.5,line=2)



##### FISHES

AllCombFishResuStd<-lapply(ecoDistSplitStd, MRM4,
                       decostandDist(log(geoDist+0.01),na.rm=TRUE),
                       decostandDist(oxDist,na.rm=TRUE),
                       decostandDist(phDist,na.rm=TRUE),
                       decostandDist(igTempDist,na.rm=TRUE),
                       decostandDist(condDist,na.rm=TRUE))

AllCombFishCoefs<-do.call(cbind,lapply(AllCombFishResuStd,function(x)x$coef[,1]))[-1,]
AllCombFishPvals<-do.call(cbind,lapply(AllCombFishResuStd,function(x)x$coef[,2]))[-1,]

rownames(AllCombFishCoefs)<-c("Distance","Oxygen","pH","Water Temperature","Conductivity")
rownames(AllCombFishPvals)<-c("Distance","Oxygen","pH","Water Temperature","Conductivity")

write.csv(round(AllCombFishCoefs,3),"Results/AllCoefsInMultipleRegressionFish.csv")
write.csv(AllCombFishPvals,"Results/AllPvalsInMultipleRegressionFish.csv")

### Simple regression model (still on matrices)

AllFishVar<-list(Geo=log(geoDist+0.01),Temp=tempMaxDist,Prec=precDryQDist,Oxygen=oxDist,pH=phDist,igTemp=igTempDist,Conductivity=condDist)

## Data inputation!!!!
# CAUTION: ONLY WRITEN SO THAT THE CODE WORKS FOR ALL GROUPS (HELPS WHEN BINDING RESULTS). DO NOT INTERPRET RESULTS THAT ARE NOT FOR FISH!!
AllFishVar<-lapply(AllFishVar,function(x){x[is.na(x)]<-sample(x[!is.na(x)],sum(is.na(x)),replace = TRUE);x})

## 
## Standardize predictor variables
AllFishVarStd<-lapply(AllFishVar,decostandDist,na.rm=TRUE)

## Run MRM for all variables individually in each group
### !!! It may take some time to run
AllIndFishResu<-lapply(AllFishVar,function(x){
  lapply(ecoDistSplitStd,MRM4,predictorDist=x) })

# Coefs
AllIndFishCoefs<-lapply(AllIndFishResu,function(x)do.call(cbind,lapply(x,function(y)y$coef[2,1])))

# Pvals
AllIndFishPvals<-lapply(AllIndFishResu,function(x)do.call(cbind,lapply(x,function(y)y$coef[2,2])))

AllIndFishCoefs<-do.call(rbind,AllIndFishCoefs)
AllIndFishPvals<-do.call(rbind,AllIndFishPvals)

rownames(AllIndFishCoefs)<-c("Distance","Temp","Prec","Oxygen","pH","igTemp","Conductivity")
rownames(AllIndFishPvals)<-c("Distance","Temp","Prec","Oxygen","pH","igTemp","Conductivity")

write.csv(round(AllIndFishCoefs,3),"Results/AllCoefsInSimpleRegressionFish.csv")
write.csv(AllIndFishPvals,"Results/AllPvalsInSimpleRegressionFish.csv")




##
varpartFishSplit<-lapply(ecoDistSplit,varpart3,Y=log(geoDist+0.01),Z=list(oxDist,phDist,condDist,igTempDist))

fishNum<-grep("Fish", labels(varpartFishSplit))

par(mfrow=c(1,1),mar=c(1,1,3,1))

lapply(fishNum,function(x){
  plot.new()
  resu<-round(varpartFishSplit[[x]]$part$indfract[,3],2)[c(1,2,3)]
  venn.plot <- draw.pairwise.venn(
    area1 = sum(resu[c(1,2)]),
    area2 = sum(resu[c(2,3)]),
    cross.area = sum(resu[c(2)]),
    category = c("", ""),
    fill = c("black", "green"),
    lty = "blank",
    cex = 4,
    cat.cex = 2,
    cat.pos = c(285, 105),
    cat.dist = 0.09,
    cat.just = list(c(-1, -1), c(1, 1)),
    ext.pos = 30,
    ext.dist = -0.05,
    ext.length = 0.85,
    ext.line.lwd = 2,
    ext.line.lty = 1,sigdigs = 2
  );
  title(main=labels(occLongSplit)[x],cex.main=4)
  dev.copy2pdf(file=paste0("Figures/",labels(occLongSplit)[x],"FishVars",".pdf"))
})


#### MAPS
# 
# for(i in 1:12){
# 
# plot(env[names(pcoaSplit[[i]][,1]),"Long"],env[names(pcoaSplit[[i]][,1]),"Lat"],pch=21,bg=heat.colors(25)[cut(pcoaSplit[[i]][,1],20)])
# 
#   }
# 
# 
# dim(pcoaSplit[[1]])
# 
# 
# plot()
# 



### 

###  Summary statistics

unlist(lapply(1:12,function(x)length(unique(env[names(pcoaSplit[[x]][,1]),]$ECO_NAME))))

