###########################################################
#### the code below is to associate the uni_plot_name column to the tables of different tables
#### unifying the names of ferns

dir()
fernsOri <- read.csv("ferndataPPBioMay2017.csv",stringsAsFactors = FALSE)
fernsGMM <- read.csv("fernpointspPPBio_afterGMM.csv",sep =",",stringsAsFactors = FALSE)
  fernsGMM <- na.omit(fernsGMM)

fernsOri$uni_plot_name <- NA

for (i in seq(nrow(fernsGMM))){
  fernsOri[fernsOri$plotcode==fernsGMM$plot[i],"uni_plot_name"] <- fernsGMM$uni_plot_name[i]
}

write.csv(fernsOri,"ferndataPPBioMay2017_afterGMM.csv")


#### unifying the name of fishes

dir()

fishOri <- read.csv("StreamFish_FishData_ToCris_afterGMM.csv",stringsAsFactors = FALSE)
fishGMM <- read.csv("StreamFish_PlotData_ToCris_afterGMM.csv",stringsAsFactors = FALSE )

################################################################################

fishOri$uni_plot_name <- NA

for (i in seq(nrow(fishGMM))){
  fishOri[fishOri$PlotID==fishGMM$IG[i],"uni_plot_name"] <- fishGMM$uni_plot_name[i]
}

write.csv(fishOri,"StreamFish_FishData_ToCris_afterGMM2.csv")

#### unifying the cupim data

cupimOri <- read.csv("PlotData_isoptera.csv",stringsAsFactors = FALSE)
cupimGMM <- read.csv("Isoptera_afterGMM.csv",stringsAsFactors = FALSE)

cupimOri$uni_plot_name <- NA

for (i in seq(nrow(cupimGMM))){
  cupimOri[cupimOri$PlotID==cupimGMM$PlotID[i],"uni_plot_name"] <- cupimGMM$uni_plot_name[i]
}

cupimOri <- cupimOri[which(!is.na(cupimOri$uni_plot_name)),]


write.csv(cupimOri,"PlotData_isoptera_afterGMM.csv")


#####################################################################################
##### fill the lat/long information
##### the code below is to associate coordinates of different groups to the PLOT_ID column
if(F){
dir()
abiotic <- read.csv("abioticdata.csv",stringsAsFactors = FALSE)
   abiotic<- abiotic[-1,-c(3:4)]
   

#ants coordinates   
ants    <- read.csv("unificada 06_07_2016_afterGMM.csv",stringsAsFactors = FALSE)
  abiotic$AntsLat  <- NA
  abiotic$AntsLong <- NA

  for (i in seq(nrow(ants))){
    abiotic[abiotic$plotID==ants$uni_plot_name[i],"AntsLat"] <- ants$Latitude_Geo[i]
  }
  for (i in seq(nrow(ants))){
    abiotic[abiotic$plotID==ants$uni_plot_name[i],"AntsLong"] <- ants$Longitude_Geo[i]
  }

#zing coordinates
zing <- read.csv("zing_all_dbase_afterGMM.csv",stringsAsFactors = FALSE)
  abiotic$ZingLat  <- NA
  abiotic$ZingLong <- NA
  
  for (i in seq(nrow(zing))){
    abiotic[abiotic$plotID==zing$uni_plot_name[i],"ZingLat"]  <- zing$lat[i]
  }
  for (i in seq(nrow(zing))){
    abiotic[abiotic$plotID==zing$uni_plot_name[i],"ZingLong"] <- zing$long[i]
  }

#fern coordinates
fern <- read.csv("fernpointspPPBio_afterGMM.csv",stringsAsFactors = FALSE)
  abiotic$fernLat <- NA
  abiotic$fernLong <- NA
  
  for (i in seq(nrow(fern))){
    abiotic[abiotic$plotID==fern$uni_plot_name[i],"fernLat"]  <- fern$LAT[i]
  }
  for (i in seq(nrow(fern))){
    abiotic[abiotic$plotID==fern$uni_plot_name[i],"fernLong"] <- fern$LONG[i]
  }
  
#fish coordinates
fish <- read.csv("StreamFish_PlotData_ToCris_afterGMM.csv",stringsAsFactors = FALSE )  
  abiotic$fishLat <- NA
  abiotic$fishLong <- NA
  
  for (i in seq(nrow(fish))){
    abiotic[abiotic$plotID==fish$uni_plot_name[i],"fishLat"]  <- fish$LAT[i]
  }
  for (i in seq(nrow(fish))){
    abiotic[abiotic$plotID==fish$uni_plot_name[i],"fishLong"] <- fish$LONG[i]
  }
  

#cupim coordinates
cupim <- read.csv("PlotData_isoptera_afterGMM2.csv",stringsAsFactors = FALSE)
  abiotic$cupimLat <- NA
  abiotic$cupimLong <- NA
  
  for (i in seq(nrow(cupim))){
    abiotic[abiotic$plotID==cupim$uni_plot_name[i],"cupimLat"]  <- cupim$LAT[i]
  }
  for (i in seq(nrow(cupim))){
    abiotic[abiotic$plotID==cupim$uni_plot_name[i],"cupimLong"] <- cupim$LONG[i]
  }
  
  write.csv(abiotic,"coord_all.csv",row.names = FALSE)
} 




###########################################################################  
##### the code below is to associate the coordinates to the soil data of the
##### table "Standardize_soil_database_AMAZONIA_ppbio_05"
  
  dir()
  abiotic_full <- read.csv("ambi.csv",stringsAsFactors = FALSE)

    soils <- read.csv("standardized_soil_database_AMAZONIA_ppbio_0a5.csv",stringsAsFactors = FALSE)
    
    teste <- abiotic_full[which(abiotic_full$plotID %in% soils$plotID),]

    #aggregate Chandless data
    cha <- subset(soils,soils$site=="chandless")
      colnames(cha)
      head(cha)
      str(cha)
    cha_agg <- aggregate(cha[,c("pH_KCl","pH_H2O","P_available.mg_kg.",
                                "Ca_cmol","Mg_cmol","Al_cmol","K_cmol","Fe_mg",
                                "Zn_mg","clay","silt","sand","SumofBases_cmol")],by=list(cha[,4]),"mean",na.action=na.omit)
    cha_agg <- as.data.frame(cha_agg) 
    colnames(cha_agg)[1] <- "plotID"
    
    mergeChandless <- merge(abiotic_full,cha_agg,by="plotID",incomparables = NA,all.x = TRUE)
    

    #average duplicates
      
    dupnames <- names(which(table(soilexceptCha$plotID)==2))
    dupsoil <- soilexceptCha[soilexceptCha$plotID %in% dupnames,]
    dupsoil$id <- seq(length(dupsoil[,1]))
    
    
    dupsoil_agg <- aggregate(dupsoil[,c("pH_KCl","pH_H2O","P_available.mg_kg.",
                                "Ca_cmol","Mg_cmol","Al_cmol","K_cmol","Fe_mg",
                                "Zn_mg","clay","silt","sand","SumofBases_cmol")],by=list(dupsoil$plotID),"mean",na.action=na.omit)
    colnames(dupsoil_agg)[1] <- "plotID"
    
    for (i in unique(dupsoil_agg$plotID)){
      mergeChandless[mergeChandless$plotID==i,c("pH_KCl","pH_H2O","P_available.mg_kg.",
                                                "Ca_cmol","Mg_cmol","Al_cmol","K_cmol","Fe_mg",
                                                "Zn_mg","clay","silt","sand","SumofBases_cmol")] <- dupsoil_agg[dupsoil_agg$plotID==i,c("pH_KCl","pH_H2O","P_available.mg_kg.",
                                                                                                                                            "Ca_cmol","Mg_cmol","Al_cmol","K_cmol","Fe_mg",
                                                                                                                                            "Zn_mg","clay","silt","sand","SumofBases_cmol")]
    }
    
    #gambiarra para inserir as colunas que viraram NA no aggregate
    duptoadd <- dupsoil[dupsoil$id[c(9,11,13,15,18,20,22,24,26,28)],c("plotID","clay","silt","sand")]
    
    for (i in unique(duptoadd$plotID)){
      mergeChandless[mergeChandless$plotID==i,c("clay","silt","sand")] <- duptoadd[duptoadd$plotID==i,c("clay","silt","sand")]
    }
    
    
    #now join the table for the rest of the data (except chandless and duplicate data
    
    soilexceptCha <- subset(soils,soils$site!="chandless")
    soilexceptCha <- soilexceptCha[,c("plotID","pH_KCl","pH_H2O","P_available.mg_kg.",
                                      "Ca_cmol","Mg_cmol","Al_cmol","K_cmol","Fe_mg",
                                      "Zn_mg","clay","silt","sand","SumofBases_cmol")]
    
    colnames(mergeChandless)
    colnames(soilexceptCha)
    soilexceptCha$plotID
    str(soilexceptCha)
    str(mergeChandless)
    
    soilexceptChadup <- subset(soilexceptCha,!soilexceptCha$plotID %in% dupsoil$plotID)
    
    soilexceptChadup$plotID
    
    for (i in unique(soilexceptChadup$plotID)){
      mergeChandless[mergeChandless$plotID==i,c("pH_KCl","pH_H2O","P_available.mg_kg.",
                                                "Ca_cmol","Mg_cmol","Al_cmol","K_cmol","Fe_mg",
                                                "Zn_mg","clay","silt","sand","SumofBases_cmol")] <- soilexceptChadup[soilexceptChadup$plotID==i,c("pH_KCl","pH_H2O","P_available.mg_kg.",
                                                                                                                                            "Ca_cmol","Mg_cmol","Al_cmol","K_cmol","Fe_mg",
                                                                                                                                            "Zn_mg","clay","silt","sand","SumofBases_cmol")]
    }
    
    #write.csv(mergeChandless,"ambi2..csv")
    
###################################################################################
# join the fish info
    dir()
    
    fishambi <- read.csv("StreamFish_PlotData_ToCris_afterGMM.csv",stringsAsFactors = FALSE)
      colnames(fishambi)
      colnames(fishambi)[3] <- "plotID"
      fishambi <- fishambi[,-c(10:12)]
      fishambi <- fishambi[!duplicated(fishambi$plotID),]

    
    mergeChandlessFish <- merge(mergeChandless, fishambi,by="plotID",incomparables = NA,all.x = TRUE )
    
    #write.csv(mergeChandlessFish,"ambi3.csv",row.names = FALSE)
    
#################################################################################    
    

    dir()
    cupimambi <- read.csv("PlotData_isoptera_afterGMM2.csv",stringsAsFactors = FALSE)
      colnames(cupimambi)
      colnames(cupimambi)[15] <- "clay"
      colnames(cupimambi)[22] <- "P_available.mg_kg."
      
      colnames(mergeChandlessFish)
      
      unique(cupimambi$GridID)
      subsetcupim <- subset(cupimambi,cupimambi$GridID=="campusUFAM")
      
      for (i in unique(subsetcupim$plotID)){
        mergeChandlessFish[mergeChandlessFish$plotID==i,c("P_available.mg_kg.",
                                                          "Ca_cmol","Mg_cmol","K_cmol","SumofBases_cmol",
                                                          "clay")] <- subsetcupim[subsetcupim$plotID==i,c("P_available.mg_kg.",
                                                                                                          "Ca_cmol","Mg_cmol","K_cmol","SumofBases_cmol","clay")]
      }
   
      unique(cupimambi$GridID)
      subsetcupim <- subset(cupimambi,cupimambi$GridID=="JP")
      
      for (i in unique(subsetcupim$plotID)){
        mergeChandlessFish[mergeChandlessFish$plotID==i,c("P_available.mg_kg.",
                                                          "Ca_cmol","Mg_cmol","K_cmol","SumofBases_cmol",
                                                          "clay")] <- subsetcupim[subsetcupim$plotID==i,c("P_available.mg_kg.",
                                                                                                          "Ca_cmol","Mg_cmol","K_cmol","SumofBases_cmol","clay")]
      } 
   
      unique(cupimambi$GridID)
      subsetcupim <- subset(cupimambi,cupimambi$GridID=="NJ")
      
      for (i in unique(subsetcupim$plotID)){
        mergeChandlessFish[mergeChandlessFish$plotID==i,c("P_available.mg_kg.",
                                                          "Ca_cmol","Mg_cmol","K_cmol","SumofBases_cmol",
                                                          "clay")] <- subsetcupim[subsetcupim$plotID==i,c("P_available.mg_kg.",
                                                                                                          "Ca_cmol","Mg_cmol","K_cmol","SumofBases_cmol","clay")]
      } 
      
      unique(cupimambi$GridID)
      subsetcupim <- subset(cupimambi,cupimambi$GridID=="IP")
      
      for (i in unique(subsetcupim$plotID)){
        mergeChandlessFish[mergeChandlessFish$plotID==i,c("P_available.mg_kg.",
                                                          "Ca_cmol","Mg_cmol","K_cmol","SumofBases_cmol",
                                                          "clay")] <- subsetcupim[subsetcupim$plotID==i,c("P_available.mg_kg.",
                                                                                                          "Ca_cmol","Mg_cmol","K_cmol","SumofBases_cmol","clay")]
      } 
      
      
      unique(cupimambi$GridID)
      subsetcupim <- subset(cupimambi,cupimambi$GridID=="IB")
      
      for (i in unique(subsetcupim$plotID)){
        mergeChandlessFish[mergeChandlessFish$plotID==i,c("P_available.mg_kg.",
                                                          "Ca_cmol","Mg_cmol","K_cmol","SumofBases_cmol",
                                                          "clay")] <- subsetcupim[subsetcupim$plotID==i,c("P_available.mg_kg.",
                                                                                                          "Ca_cmol","Mg_cmol","K_cmol","SumofBases_cmol","clay")]
      } 
      
      unique(cupimambi$GridID)
      subsetcupim <- subset(cupimambi,cupimambi$GridID=="TE")
      
      for (i in unique(subsetcupim$plotID)){
        mergeChandlessFish[mergeChandlessFish$plotID==i,c("P_available.mg_kg.",
                                                          "Ca_cmol","Mg_cmol","K_cmol","SumofBases_cmol",
                                                          "clay")] <- subsetcupim[subsetcupim$plotID==i,c("P_available.mg_kg.",
                                                                                                          "Ca_cmol","Mg_cmol","K_cmol","SumofBases_cmol","clay")]
      } 
      
    
      unique(cupimambi$GridID)
      subsetcupim <- subset(cupimambi,cupimambi$GridID=="FEUFAM")
      
      for (i in unique(subsetcupim$plotID)){
        mergeChandlessFish[mergeChandlessFish$plotID==i,c("P_available.mg_kg.",
                                                          "Ca_cmol","Mg_cmol","K_cmol","SumofBases_cmol",
                                                          "clay")] <- subsetcupim[subsetcupim$plotID==i,c("P_available.mg_kg.",
                                                                                                          "Ca_cmol","Mg_cmol","K_cmol","SumofBases_cmol","clay")]
      }  
      
      
      unique(cupimambi$GridID)
      subsetcupim <- subset(cupimambi,cupimambi$GridID=="ParnaJau")
      
      for (i in unique(subsetcupim$plotID)){
        mergeChandlessFish[mergeChandlessFish$plotID==i,c("P_available.mg_kg.",
                                                          "Ca_cmol","Mg_cmol","K_cmol","SumofBases_cmol",
                                                          "clay")] <- subsetcupim[subsetcupim$plotID==i,c("P_available.mg_kg.",
                                                                                                          "Ca_cmol","Mg_cmol","K_cmol","SumofBases_cmol","clay")]
      }  
      
      
      unique(cupimambi$GridID)
      subsetcupim <- subset(cupimambi,cupimambi$GridID=="RFAD")
      
      for (i in unique(subsetcupim$plotID)){
        mergeChandlessFish[mergeChandlessFish$plotID==i,c("P_available.mg_kg.",
                                                          "Ca_cmol","Mg_cmol","K_cmol","SumofBases_cmol",
                                                          "clay")] <- subsetcupim[subsetcupim$plotID==i,c("P_available.mg_kg.",
                                                                                                          "Ca_cmol","Mg_cmol","K_cmol","SumofBases_cmol","clay")]
      }  
      
      write.csv(mergeChandlessFish,"ambi3.3.csv",row.names = FALSE)
      
######################################################################################################################
#### gabi produce a ambi4 file and I will merge it with the margechandlessFish to creat the ambi5
      
      dir()
      ambi4gabi <- read.csv("ambi4_gabi.csv",stringsAsFactors = FALSE)
        colnames(ambi4gabi)
        ambi4gabi <- ambi4gabi[,c("plotID","sa.latlong.treecover", "sa.latlong.evergreen","chelsa6","bioregion")]
      
      mergeChandlessFishGabi <- merge(mergeChandlessFish, ambi4gabi,by="plotID",incomparables = NA,all.x = TRUE )
      
      
      write.csv(mergeChandlessFish,"ambi5.csv",row.names = FALSE)
      
      ambi5 <- mergeChandlessFish


######################################################################################################################
##### adding missing information about RFAD soil texture

      dir()
      ambi6 <- read.csv("ambi6.csv",stringsAsFactors = FALSE)
      granuRFAD <- read.csv("granu_ducke_afterGMM.csv",stringsAsFactors = FALSE)
        colnames(granuRFAD)
        colnames(granuRFAD)[3:5] <- c("clay","silt","sand"); colnames(granuRFAD)
      
      
      
      for(i in unique(granuRFAD$plotID))  {
      ambi6[ambi6$plotID==i,c("clay","silt","sand")] <- granuRFAD[granuRFAD$plotID==i,c("clay","silt","sand")]
      }
      
        ambi6[grep("RFAD_L",ambi6$plotID),c("plotID","clay","silt","sand")]
        
        write.csv(ambi6,"ambi7.csv",row.names = FALSE)
        ambi7 <- ambi6
        colnames(ambi7)["pH_H2O"]
        
        ##### adding missing information about RFAD coordinates
      
       ambi7[ambi7$plotID==c("RFAD_L1_3000"),c("Lat","Long")] <- c(-2.92367,-59.94254)
       ambi7[ambi7$plotID==c("RFAD_L1_7500"),c("Lat","Long")] <- c(-2.91735,-59.90273)
       ambi7[ambi7$plotID==c("RFAD_L2_5500"),c("Lat","Long")] <- c(-2.929,-59.91908)
       ambi7[ambi7$plotID==c("RFAD_L3_7500"),c("Lat","Long")] <- c(-2.93511,-59.89992)
       ambi7[ambi7$plotID==c("RFAD_L5_7500"),c("Lat","Long")] <- c(-2.95278,-59.89711)
       ambi7[ambi7$plotID==c("RFAD_L7_5500"),c("Lat","Long")] <- c(-2.97326,-59.91206)
       ambi7[ambi7$plotID==c("RFAD_L8_7500"),c("Lat","Long")] <- c(-2.97929,-59.8929)
       ambi7[ambi7$plotID==c("RFAD_L9_5500"),c("Lat","Long")] <- c(-2.99102,-59.90933)
       ambi7[ambi7$plotID==c("RFAD_L9_7500"),c("Lat","Long")] <- c(-2.98813,-59.89149)
       
       ##### adding missing information about RFAD soil nutrients
       
       dir()
       duckequimica <- read.csv("ducke_quimica_revisado2a_aftergmm.csv",stringsAsFactors = FALSE)
       colnames(duckequimica)
       
       colnames(ambi7)
       teste <- ambi7[grep("RFAD_L",ambi7$plotID),c("plotID","pH_H2O","P_available.mg_kg.","Ca_cmol",
                                           "Mg_cmol","Al_cmol","K_cmol","Fe_mg","Zn_mg")]
       length(teste$plotID)
       
       for(i in unique(duckequimica$plotID)){
         ambi7[ambi7$plotID==i,c("pH_H2O","P_available.mg_kg.","Ca_cmol",
                                 "Mg_cmol","Al_cmol","K_cmol","Fe_mg","Zn_mg")] <- duckequimica[duckequimica$plotID==i,c("pH_H2O","P_available.mg_kg.","Ca_cmol",
                                                                                                                                                                                                                                                  "Mg_cmol","Al_cmol","K_cmol","Fe_mg","Zn_mg")]
       }
       
       write.csv(ambi7,"ambi7.csv",row.names = FALSE)
  
##########################################################################################################################
#### adding madeireira mil data
       
       dir()
       ambi7.2 <- read.csv("ambi7.2.csv",stringsAsFactors = FALSE)
       colnames(ambi7.2)
        
       ambi7.2 <- ambi7.2[!duplicated(ambi7.2$plotID),]
       length(ambi7.2[,1])
       
       write.csv(ambi7.2,"ambi7.3.csv",row.names = FALSE)
       
##########################################################################################################################
####### adding Cauame envi data
       
       dir()
       cauameF <- read.csv("menger.169.1-cauame_fisica_original2_0a5.csv",stringsAsFactors = FALSE)
        colnames(cauameF)
       cauameQ <- read.csv("menger.168.1-cauame_quimica_revisado2_0a5.csv", stringsAsFactors = FALSE)
        colnames(cauameQ)
       
       ambi8 <- read.csv("ambi8.csv", stringsAsFactors = FALSE)
       
       for(i in unique(cauameF$plotID)){
         ambi8[ambi8$plotID==i,c("clay","silt" ,"sand")] <- cauameF[cauameF$plotID==i,c("clay" ,"silt" ,"sand")]
       }
       
       for(i in unique(cauameQ$plotID)){
         ambi8[ambi8$plotID==i,c("pH_H2O","P_available.mg_kg.","Ca_cmol",
                                 "Mg_cmol","Al_cmol","K_cmol","Fe_mg","Zn_mg")] <- cauameQ[cauameQ$plotID==i,c("pH_H2O","P_available.mg_kg.","Ca_cmol",
                                                                                                               "Mg_cmol","Al_cmol","K_cmol","Fe_mg","Zn_mg")]
       }
       
      colnames(ambi8)
      write.csv(ambi8, "ambi8.csv",row.names = FALSE)
       
      #write.csv(ambi8,"ambi8.csv",row.names = FALSE)
      
###########################################################################################################
#### reproject data
library(sp)
library(rgdal)      
      
      ambi10 <- read.csv("ambi10.csv",stringsAsFactors = F)
      colnames(ambi10)
      ambi10_temp <- na.omit(ambi10[,c("Lat","Long")])
      str(ambi10_temp)
      ambi10_temp <- ambi10_temp[,c("Lat","Long")]
      
      
      coordinates(ambi10_temp) <- ~Lat+Long
      utm20S <- CRS("+proj=utm +zone=20 +south +ellps=WGS84 +datum=WGS84 +units=m +no_defs" )
      proj4string(ambi10_temp) <- utm20S

###########################################################################################################
# adicionar bats
      
      dir()
      
      batssp <- read.csv("Bats_SpeciesData.csv",stringsAsFactors = F)
        colnames(batssp)
        batssp$plotID <- NA
      batsambi <- read.csv("Bats_PlotID_afterGMM.csv",stringsAsFactors = F)
        colnames(batsambi)
      
      for(i in unique(batsambi$PlotID_old)){
        batssp[batssp$PlotID==i,"plotID"] <- batsambi[batsambi$PlotID_old==i,"plotID"]
      }
        unique(batssp$plotID)
      
      write.csv(batssp,"Bats_SpeciesData_afterGMM2.csv",row.names = F)
      
# find bat duplicates in the ambi11
      
      ambi11 <- read.csv("ambi11.csv",stringsAsFactors = F)
      length(ambi11$plotID)
      
      teste <- ambi11[duplicated(ambi11$plotID),]
      length(teste$plotID)
      
      write.csv(ambi11, "ambi11.csv",row.names = F)
      
      
      ambi11_temp <- read.csv("ambi11_temp.csv",stringsAsFactors = F)
      
      teste <- ambi11_temp[ambi11_temp$grid=="UHE_Santo_Antonio",]
      
      ambi11.2 <- ambi11_temp[duplicated(ambi11_temp$plotID),"plotID"]
      
      teste <- ambi11.2[ambi11.2$grid=="UHE_Santo_Antonio",]
      