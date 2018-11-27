#Script to download the ERSST from IRI data library for a specific season 
#and specific area in CPT format

#Created by: Lizeth Llanos and Diego Agudelo
#Date: February 2018

options(timeout=180)

if(require(stringr)==FALSE){install.packages("stringr",dependencies = TRUE)}
library("stringr")
if(require(R.utils)==FALSE){install.packages("R.utils",dependencies = TRUE)}
library(R.utils)

# Generación de estructura de carpetas ------------------------------------------

main_dir <- "D:/OneDrive - CGIAR/Desktop/test" # Modifique esta línea de acuerdo a su directorio de trabajo
nom_c <- "Corrida_1" # Modifique esta línea con el nombre que desee

dir.create(paste0(main_dir,"/",  nom_c,"/input/sst_ersst"), recursive = T)
dir.create(paste0(main_dir,"/", nom_c,"/input/stations"), recursive = T)
dir_save <- paste0(main_dir,"/", nom_c,"/sst_ersst")


download_ERSST_CPT=function(firs_year,last_year,i_month,l_season,dir_save,m_for,l_for,area1){
  
  trimestrel <- i_month:(i_month+l_season-1)
  fores <- m_for:(m_for+l_for-1)
  if(sum(fores>12)>0)fores[which(fores>12)]=fores[which(fores>12)]-12
  if(sum(trimestrel>12)>0)trimestrel[which(trimestrel>12)]=trimestrel[which(trimestrel>12)]-12
  route <- paste0("http://iridl.ldeo.columbia.edu/SOURCES/.NOAA/.NCDC/.ERSST/.version4/.sst/T/%28", month.abb[trimestrel[1]] ,"%20", firs_year ,"%29%28",  month.abb[trimestrel[l_season]] ,"%20", last_year ,"%29RANGEEDGES/T/", l_season ,"/boxAverage/T/12/STEP/Y/%28",area1[4],"%29%28",area1[3],"%29RANGEEDGES/X/%28",area1[1],"%29%28",area1[2],"%29RANGEEDGES/-999/setmissing_value/Y/high/low/RANGE/%5BX/Y%5D%5BT%5Dcptv10.tsv.gz")
  path_save <- paste0(dir_save,"/",paste(month.abb[fores],collapse = "-"),"_",paste(month.abb[trimestrel],collapse = "-"),".tsv.gz")
  download.file(route,path_save)
  gunzip(path_save)
  
  return("Successful download")
  
}

#Ejemplo para descargar SON para pronosticar DEF de la ERSST v4.0 

area1 <- c(0,359,-30,30) #xmin, xmax, ymin, ymax
i_month <- 9 #First month to download
l_season <- 3 #Length season (1, 2, 3. meses)
m_for <- 12 #Month to forecast
l_for <- 3 #Length forescast 
firs_year <- 1981 #Initial year
last_year <- 2018 #Last year

download_ERSST_CPT(firs_year,last_year,i_month ,l_season ,dir_save,m_for,l_for,area1)
