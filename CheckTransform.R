##### Script to check and transform data from the PPBio sites

# Import plot data (PlotID, Env, LatLong, etc)
PlotData<-read.csv("Data/PlotData.csv")

#BioFiles<-list.files(recursive = TRUE)
BioFiles<-c("Data/Isoptera.csv","Data/Zingiberales.csv")

for(i in 1:length(BioFiles)){

BioData<-read.csv(files[i])

# Check if any of the biological data is in an unidentified plot in the plot data (PlotID)
# Should return empty
cat(paste0(unique(BioData$PlotID[!BioData$PlotID%in%PlotData$PlotID]),collapse="\n"))

print(i)
}
##

