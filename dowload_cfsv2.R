#Script to download SST from CFSv2 from IRI data library in CPT format

#Created by: Lizeth Llanos and Diego Agudelo
#Date: February 2018

options(timeout=180)

if(require(stringr)==FALSE){install.packages("stringr",dependencies = TRUE)}
library("stringr")

#Function to download SST for one area
download_CFSV2_CPT_1=function(firs_year,last_year,i_month,ic,dir_save,area1){
  
  lead <- i_month-ic
  if(lead<0)lead <- lead + 12
  route <- paste0("http://iridl.ldeo.columbia.edu/SOURCES/.NOAA/.NCEP/.EMC/.CFSv2/.ENSEMBLE/.OCNF/.surface/.TMP/SOURCES/.NOAA/.NCEP/.EMC/.CFSv2/.REALTIME_ENSEMBLE/.OCNF/.surface/.TMP/appendstream/350/maskge/S/%280000%201%20",month.abb[ic],"%20",firs_year,"-",last_year,"%29VALUES/L/",lead,".5/",lead+2,".5/RANGE%5BL%5D//keepgrids/average/M/1/24/RANGE%5BM%5Daverage/Y/%28",area1[4],"%29%28",area1[3],"%29RANGEEDGES/X/%28",area1[1],"%29%28",area1[2],"%29RANGEEDGES/-999/setmissing_value/%5BX/Y%5D%5BS/L/add%5Dcptv10.tsv.gz")
  
  trimestrel <- (ic+lead):(ic+lead+2)
  if(sum(trimestrel>12)>0)trimestrel[which(trimestrel>12)]=trimestrel[which(trimestrel>12)]-12
  path_save <- paste0(dir_save,"/",month.abb[ic],"_",paste(month.abb[trimestrel],collapse = "-"),".tsv.gz")
  download.file(route,path_save)
  return(paste("Successful download",path_save))
}



#Ejemplo para descargar Feb_Abr-May-Jun (trimestre AMJ con condición inicial 
#en Febrero para el área xmin =0, xmax =359, ymin= -30, ymax =30)

area1 <- c(0,359,-30,30)

i_month <- 4
ic <- 2
firs_year <- 1981
last_year <- 2018
dir_save <- "C:/Users/lllanos/Desktop/ejercicios cpt/CFSv2"

download_CFSV2_CPT_1(firs_year,last_year,i_month,ic,dir_save,area1)


# Function to download SST from CFSv2 for two areas -----------------------


download_CFSV2_CPT=function(firs_year,last_year,i_month,ic,dir_save,area1,area2){
  
  lead <- i_month-ic
  if(lead<0)lead <- lead + 12
  route <- paste0("http://iridl.ldeo.columbia.edu/SOURCES/.NOAA/.NCEP/.EMC/.CFSv2/.ENSEMBLE/.OCNF/.surface/.TMP/SOURCES/.NOAA/.NCEP/.EMC/.CFSv2/.REALTIME_ENSEMBLE/.OCNF/.surface/.TMP/appendstream/350/maskge/S/%280000%201%20",month.abb[ic],"%20",firs_year,"-",last_year,"%29/VALUES/L/",lead,".5/",lead+2,".5/RANGE/%5BL%5D//keepgrids/average/M/1/24/RANGE/%5BM%5Daverage/X/",area1[1],"/",area1[2],"/flagrange/Y/",area1[3],"/",area1[4],"/flagrange/add/1/flaggt/X/",area2[1],"/",area2[2],"/flagrange/Y/",area2[3],"/",area2[4],"/flagrange/add/1/flaggt/add/mul/0/setmissing_value/-999/replaceNaN/%5BX/Y%5D%5BS/L/add/%5Dcptv10.tsv.gz")
  trimestrel <- (ic+lead):(ic+lead+2)
  if(sum(trimestrel>12)>0)trimestrel[which(trimestrel>12)]=trimestrel[which(trimestrel>12)]-12
  path_save <- paste0(dir_save,"/",month.abb[ic],"_",paste(month.abb[trimestrel],collapse = "-"),".tsv.gz")
  download.file(route,path_save)
  return("Successful download")
  
}

######## Run #############

area1 <- c(180, 270,-13,13)
area2 <- c(45 , 105,-13,13)
firs_year <- 1981
last_year <- 2016
dir_save <- "C:/Users/dagudelo/Desktop"
i_month <- 3
ic <- 2

download_CFSV2_CPT(firs_year,last_year,i_month,ic,dir_save,area1,area2)

