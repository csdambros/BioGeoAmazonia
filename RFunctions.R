## Auxiliary functions

adjustPredictor<-function(ecoDist,predictorDist){
  
  predMatrixFull<-as.matrix(predictorDist)
  match1<-match(labels(ecoDist),rownames(predMatrixFull))
  predMatrix<-predMatrixFull[match1,match1]
  
  predDist<-as.dist(predMatrix)
  
  return(predDist)
}

# Function to standardize variables in distance format

decostandDist<-function(dist,method="standardize",...){
  
  require(vegan)
  
  lab<-labels(dist)
  
  Mat<-as.matrix(dist)
  nMat<-nrow(Mat)
  resuMat<-matrix(decostand(as.vector(Mat),method=method,...),nMat,nMat)
  
  rownames(resuMat)<-lab
  colnames(resuMat)<-lab
  
  resu<-as.dist(resuMat)
  return(resu)
    
}



MRM2<-function(ecoDist,predictorDist){
  
  predDist<-adjustPredictor(ecoDist,predictorDist)
  
  mrmResu<-MRM (ecoDist~predDist)
  
  return(mrmResu)
}


MRM3<-function(ecoDist,...){
 
  x <- list(...)

  predDist<-lapply(x,adjustPredictor,ecoDist=ecoDist)
  
  formula2MRM<-formula(paste("ecoDist~",paste0(paste0("predDist[[",1:length(predDist),"]]"),collapse="+")))
  
  if(sum(unlist(lapply(predDist,function(x)sum(x,na.rm=TRUE)==0)))){}else{
  
  mrmResu<-MRM (formula2MRM)

  return(mrmResu)
}
}

# Function to run regression with standardized variables
MRM4<-function(ecoDist,...){
  
  x <- list(...)
  
  predDist<-lapply(x,adjustPredictor,ecoDist=ecoDist)
  
  predDist<-lapply(predDist,decostandDist,na.rm=TRUE)
  
  formula2MRM<-formula(paste("ecoDist~",paste0(paste0("predDist[[",1:length(predDist),"]]"),collapse="+")))
  
  if(sum(unlist(lapply(predDist,function(x)sum(x,na.rm=TRUE)==0)))){}else{
    
    mrmResu<-MRM (formula2MRM)
    
    return(mrmResu)
  }
}

# Extract fitted values and residuals from MRM
# MRM is the resulf of the MRM fit (what the MRM function spills)
# Response is a distance matrix
# Predictors is a data frame with all predictor variables as columns or a distance matrix

residuals.mrm<-function(MRM,response,predictors){
  
  predictors<-cbind(1,predictors)
  MRMFit<-predictors%*%MRM$coef[,1] # Calculate fitted (predicted) values
  MRMRes<-response-MRMFit # Calculate residuals
  
  resu<-list(predicted.values=MRMFit,residuals=MRMRes)
  return(resu)
}

#MRM5
# Same as 4, but includes original data as an object for posterior use





varpart3<-function(X,Y,Z){
#X<-BirdDist
#Y<-BirdGeoDist
#Z<-list(BirdTempDist)
  

  if(class(Y)=="dist"){
    
    Y<-list(Y)
    
  }
  
  if(class(Z)=="dist"){
    
    Z<-list(Z)
    
  }
  
  Y<-lapply(Y,adjustPredictor,ecoDist = X)
  Y<-do.call(cbind,Y)
  
  Z<-lapply(Z,adjustPredictor,ecoDist = X)
  Z<-do.call(cbind,Z)

rm.na<-!(rowSums(is.na(cbind(X,Y,Z)))>0)

varpart(cbind(X)[rm.na,],Y[rm.na,],Z[rm.na,])

}
  




#plot(varpart3(BirdDist,geoDist,list(treeCoverDist)))

# L = (1+sqrt(1+8*n))/45 # If a vector is provided, calculate the size of original matrix

# Rarefy data for those groups/sites with more sampling than others
# Example, data was sampled in 3 years in one site and 2 years in others
# Produces the same results as randomly sampling 2 years and counting species occurrences thousands of times
# Can be used as input for jaccard, sorensen, etc.

# N is the maximum number of samples taken in a site (ex. 3 years for sites with 3 years and 2 years for the rest)
# Ni Number of samples a species was found (ex. max o 3 for sites sampled 3 times/years)
# n number of samples in the standardized data (ie. all sites will be compared as if they were sampled in 2 years, then n=2)
# # If N and n are the same, then results must be the same as the original data
# 
# POccur<-function(N,Ni,n){1-(factorial(N-n)/factorial(N))*(factorial(N-Ni)/factorial(N-Ni-n))}# Function to calculate prob of occurrence in n sections given Ni occurrences in N sections (without replacement)
# 
# 
# 
# 
# POccur(10,,10)
# 
# N=4
# Ni=3
# n=5
# 
# 1-(factorial(N-n)/factorial(N))*(factorial(N-Ni)/factorial(N-Ni-n))
# 
# specaccum()
# rarefy()
# 
# 
# 5*(1-(factorial(Ni)/factorial(Ni-(N-n+1)))/(factorial(N)/factorial(N-n+1)))
# 
# POccur2<-function(A,a,ni){
#   
#   (factorial(A-ni)/(factorial(a)*factorial(A-ni-a)))/
#     (factorial(A)/1)
#   
#   
# }
# 
# resu<-{}
# 
# for(i in 1:10000){
# resu[i]<-sum((sample(1:5)<=4)[1:4])
# }
# mean(resu)
# 
# for(i in 1:10000){
#   resu[i]<-sum((sample(1:5,replace = TRUE)<=3)[1:4])
# }
# mean(resu)
# 
# POccur(N = 5,Ni = 5,n = 1)
# 
# 
# pbinom(4,5,3/5)
# 
# P
# 
# 
# POccur<-function(N,Ni,n){sum((sample(1:N)<=Ni)[1:n])}# Function to calculate prob of occurrence in n sections given Ni occurrences in N sections (without replacement)
# 
# 



