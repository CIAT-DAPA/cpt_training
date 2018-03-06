#Script to download the ERSST from IRI data library for a specific season 
#and specific area in CPT format

#Created by: Lizeth Llanos and Diego Agudelo
#Date: February 2018


if(require(stringr)==FALSE){install.packages("stringr",dependencies = TRUE)}
library("stringr")

download_ERSST_CPT=function(firs_year,last_year,i_month,l_season,dir_save,m_for,area1){
  
  trimestrel <- i_month:(i_month+l_season-1)
 
   if(sum(trimestrel>12)>0)trimestrel[which(trimestrel>12)]=trimestrel[which(trimestrel>12)]-12
  route <- paste0("http://iridl.ldeo.columbia.edu/SOURCES/.NOAA/.NCDC/.ERSST/.version4/.sst/T/%28", month.abb[trimestrel[1]] ,"%20", firs_year ,"%29%28",  month.abb[trimestrel[l_season]] ,"%20", last_year ,"%29RANGEEDGES/T/", l_season ,"/boxAverage/T/12/STEP/Y/%28",area1[4],"%29%28",area1[3],"%29RANGEEDGES/X/%28",area1[1],"%29%28",area1[2],"%29RANGEEDGES/-999/setmissing_value/Y/high/low/RANGE/%5BX/Y%5D%5BT%5Dcptv10.tsv.gz")
  m_for_final=str_pad(m_for, 2, pad = "0")
  path_save <- paste0(dir_save,"/",m_for_final,"_",paste(month.abb[trimestrel],collapse = "-"),".tsv.gz")
  download.file(route,path_save)
 
   return("Successful download")
  
}


#Example to download January from ERSST v4.0 

area1 <- c(0,359,-30,30) #xmin, xmax, ymin, ymax
i_month <- 1 #First month to download
l_season <- 1 #Length season (1, 2, 3. meses)
m_for <- 4 #Month to forecast
firs_year <- 1981 #Initial year
last_year <- 2018 #Last year
dir_save <- "C:/Users/dagudelo/Desktop/ERSST" #Path 

download_ERSST_CPT(firs_year,last_year,i_month ,l_season ,dir_save,m_for ,area1)
