

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
  
  distCopy<-dist
#  lab<-labels(dist)
  
#  Mat<-as.matrix(dist)
#  nMat<-nrow(Mat)
#  resuMat<-matrix(decostand(as.vector(Mat),method=method,...),nMat,nMat)
  distStd<-decostand(as.vector(dist),"standardize",...)
  distCopy[1:length(distCopy)]<-distStd
#  rownames(resuMat)<-lab
#  colnames(resuMat)<-lab
  
#  resu<-as.dist(resuMat)
  return(distCopy)
    
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
  
  names(x)<-paste0("p",1:length(x),".",names(x))

  predDist<-lapply(x,adjustPredictor,ecoDist=ecoDist)
  
  predDist<-lapply(predDist,decostandDist,na.rm=TRUE)
  
  
  formula2MRM<-formula(paste("ecoDist~",paste0(names(predDist),collapse="+")))

  predFrame<-as.data.frame(do.call(cbind,predDist))
  
  if(sum(unlist(lapply(predDist,function(x)sum(x,na.rm=TRUE)==0)))){}else{
    
    mrmResu<-MRM (formula2MRM,data =predFrame)
    
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



varpart4<-function(X,Y,Z,H){
#X<-BirdDist
#Y<-BirdGeoDist
#Z<-list(BirdTempDist)
  require(ecodist)
  require(vegan)
  
  

  if(class(Y)=="dist"){
    
    Y<-list(Y)
    
  }
  
  if(class(Z)=="dist"){
    
    Z<-list(Z)
    
  }
  
  if(!is.null(H)){
  
  if(class(H)=="dist"){
    
    H<-list(H)
    
  }
    
    H<-lapply(H,adjustPredictor,ecoDist = X)
    H<-do.call(cbind,H)
    
  }
  
  
  Y<-lapply(Y,adjustPredictor,ecoDist = X)
  Y<-do.call(cbind,Y)
  
  Z<-lapply(Z,adjustPredictor,ecoDist = X)
  Z<-do.call(cbind,Z)

rm.na<-!(rowSums(is.na(cbind(X,Y,Z,H)))>0)

varpart(cbind(X)[rm.na,],Y[rm.na,],Z[rm.na,],H[rm.na,])

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

## plot varpart

shArea<-function(d,r,R){
  d2<-d^2
  r2<-r^2
  R2<-R^2
  alpha = acos((d2 + r2 - R2) / (2*d*r))
  beta = acos((d2 + R2 - r2) / (2*d*R))
  SharedArea<-r2 * alpha + R2 * beta -0.5 * (r2 * sin(2*alpha) + R2 * sin(2*beta))
  return(SharedArea)
}

optshArea<-function(d,sharea,r,R){
  abs(sharea-shArea(d,r,R))
}

optShd<-function(sharea,r,R){
  rmax<-max(r,R)
  rmin<-min(r,R)
  interval<-c(rmax-rmin+(rmax-rmin)*0.00001,r+R-(r+R)*0.00001) 
  
  resu<-optimize(interval = interval,optshArea,sharea=sharea,r=r,R=R)
  
  return(resu$minimum)
}



plot.varpart2<-function(varpart,col=c(1,2),lty=2,lty1=0,xlim=c(4.5,5.5),ylim=c(4,6),main="",border=col,labels=TRUE,...){
  
  library(plotrix)
  
  
  p1<-round(varpart$part$indfract[c(1,2),3],2)
  sp1<-sum(p1[p1>0])
  
  p2<-round(varpart$part$indfract[c(2,3),3],2)
  sp2<-sum(p2[p2>0])
  
  p12<-round(varpart$part$indfract[c(2),3],2)
  sp12<-sum(p12[p12>0])
  
  r1<-sqrt(sp1/pi)
  r2<-sqrt(sp2/pi)

  d1<-optShd(sp12,r1,r2)
  ###
  
  plot.new();plot.window(xlim=xlim,ylim=ylim)
  title(main=main)
  draw.circle(5,5,r1,col=col[1],lty=lty1,border=border[1],...)
  draw.circle(5+d1,5,r2,col=col[2],lty=lty1,border=border[2],...)
  
  draw.circle(5,5,r1,lty=lty)
  draw.circle(5+d1,5,r2,lty=lty)

  text((10+d1-r2+r1)/2,5,sp12)
  text(5-r1-0.1,5+0.1,sp1-sp12)
  text(5+d1+r2+0.1,5+0.1,sp2-sp12)
  
  
   #distances<-dist(cbind(c(5+d1,5,5-x),c(5,5,5-y)))
  
  #return(distances)
  
}

plot.varpart3<-function(varpart,col=c(1,2),lty=2,lty1=0,xlim=c(4.5,5.5),ylim=c(4,6),main="",border=col,labels=TRUE,tdist=0.1,values=TRUE,...){
  
# varpart<-varpartAllSplit[[9]]
  
  vtable<-round(varpart$part$indfract[,3],2)
  vtable[vtable<=0]<-NA
  
  p1<-vtable[c(1,4,6,7)]
  sp1<-sum(p1[p1>0],na.rm = TRUE)
  
  p2<-vtable[c(2,4,5,7)]
  sp2<-sum(p2[p2>0],na.rm=TRUE)
  
  p3<-vtable[c(3,5,6,7)]
  sp3<-sum(p3[p3>0],na.rm=TRUE)
  
  p12<-vtable[c(4,7)]
  sp12<-sum(p12[p12>0],na.rm=TRUE)
  
  p13<-vtable[c(6,7)]
  sp13<-sum(p13[p13>0],na.rm=TRUE)
  
  p23<-vtable[c(5,7)]
  sp23<-sum(p23[p23>0],na.rm=TRUE)
  
  p123<-vtable[7]
  sp123<-sum(p123[p123>0],na.rm = TRUE)
  
  r1<-sqrt(sp1/pi)
  r2<-sqrt(sp2/pi)
  r3<-sqrt(sp3/pi)
  
  d1<-optShd(sp12,r1,r2)
  shArea(d1,r1,r2)
  
  d2<-optShd(sp13,r1,r3)
  shArea(d2,r1,r3)
  
  d3<-optShd(sp23,r2,r3)
  shArea(d3,r2,r3)
  
  ###
  
  x<-(d2^2-d1^2-d3^2)/(2*d1)
  if((x^2)>(d3^2)){x<-d3}
  y<-sqrt(d3^2-x^2)
  
  library(plotrix)
  plot.new();plot.window(xlim=xlim,ylim=ylim)
#  plot.new();plot.window(xlim=c(4,6),ylim=c(4,6))
  title(main=main)
  #draw.circle(5,5,r2,col=col[2],lty=lty1,border=border[2],...)
  #draw.circle(5+d1,5,r1,col=col[1],lty=lty1,border=border[1],...)
  #draw.circle(5-x,5-y,r3,col=col[3],lty=lty1,border=border[3],...)
  
  polygon(sin(0:360*pi/180)*r1+5+d1,cos(0:360*pi/180)*r1+5,col=col[1],lty=lty1,border=border[1],...)
  polygon(sin(0:360*pi/180)*r2+5,cos(0:360*pi/180)*r2+5,col=col[2],lty=lty1,border=border[2],...)
  polygon(sin(0:360*pi/180)*r3+5-x,cos(0:360*pi/180)*r3+5-y,col=col[3],lty=lty1,border=border[3],...)

  # draw.circle(5,5,r2,lty=lty)
  # draw.circle(5+d1,5,r1,lty=lty)
  # draw.circle(5-x,5-y,r3,lty=lty)

  polygon(sin(0:360*pi/180)*r1+5+d1,cos(0:360*pi/180)*r1+5,lty=lty,border=border[1])
  polygon(sin(0:360*pi/180)*r2+5,cos(0:360*pi/180)*r2+5,lty=lty,border=border[2])
  polygon(sin(0:360*pi/180)*r3+5-x,cos(0:360*pi/180)*r3+5-y,lty=lty,border=border[3])
  
  a<-atan(0)
  a2<-atan2(y,x)
  a3<-atan2(-y,(-d1-x))
  
  #points(5-x,5-y)
  
  if(values){
  text(cos(a)*(-r2-tdist)+5,sin(a)*(-r2-tdist)+5,vtable[1])
  text(cos(a)*(r1+tdist)+5+d1,sin(a)*(r1+tdist)+5,vtable[2])
  text(cos(a2)*-(r3+tdist)+5-x,sin(a2)*-(r3+tdist)+5-y,vtable[3])
  
  text((cos(a)*(r2-r1)+d1)/2+5,5,vtable[4])
  text((cos(a2)*(-r2+r3)-x)/2+5,(sin(a2)*(-r2+r3)-y)/2+5,vtable[5])
  text((cos(a3)*(r1-r3)+d1-x)/2+5,(sin(a3)*(r1-r3)-y)/2+5,vtable[6])

}
  #points(5+d1+cos(a3)*seq(0,1,0.01),5+sin(a3)*seq(0,1,0.01),cex=.1)
  

  distances<-dist(cbind(c(5+d1,5,5-x),c(5,5,5-y)))
  
  return(distances)
  
}


