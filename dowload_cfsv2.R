#Script to download SST from CFSv2 from IRI data library in CPT format

#Created by: Lizeth Llanos and Diego Agudelo
#Date: February 2018

options(timeout=180)

if(require(stringr)==FALSE){install.packages("stringr",dependencies = TRUE)}
library("stringr")
if(require(R.utils)==FALSE){install.packages("R.utils",dependencies = TRUE)}
library(R.utils)


# Generación de estructura de carpetas ------------------------------------------

main_dir <- "D:/OneDrive - CGIAR/Desktop/test" # Modifique esta línea de acuerdo a su directorio de trabajo
nom_c <- "Corrida_2" # Modifique esta línea con el nombre que desee

dir.create(paste0(main_dir,"/", nom_c,"/input/sst_cfsv2"), recursive = T)
dir.create(paste0(main_dir,"/", nom_c,"/input/stations"), recursive = T)
dir_save <- paste0(main_dir,"/", nom_c,"/input/sst_cfsv2")


# Función para descargar una área de la TSM del modelo CFSv2 --------------
download_CFSV2_CPT_1=function(firs_year,last_year,i_month,ic,dir_save,area1, lg){
  lg_s <- lg -1
  lead <- i_month-ic
  if(lead<0)lead <- lead + 12
  route <- paste0("http://iridl.ldeo.columbia.edu/SOURCES/.NOAA/.NCEP/.EMC/.CFSv2/.ENSEMBLE/.OCNF/.surface/.TMP/SOURCES/.NOAA/.NCEP/.EMC/.CFSv2/.REALTIME_ENSEMBLE/.OCNF/.surface/.TMP/appendstream/350/maskge/S/%280000%201%20",month.abb[ic],"%20",firs_year,"-",last_year,"%29VALUES/L/",lead,".5/",lead+lg_s,".5/RANGE%5BL%5D//keepgrids/average/M/1/24/RANGE%5BM%5Daverage/Y/%28",area1[4],"%29%28",area1[3],"%29RANGEEDGES/X/%28",area1[1],"%29%28",area1[2],"%29RANGEEDGES/-999/setmissing_value/%5BX/Y%5D%5BS/L/add%5Dcptv10.tsv.gz")
  
  trimestrel <- (ic+lead):(ic+lead+lg_s)
  if(sum(trimestrel>12)>0)trimestrel[which(trimestrel>12)]=trimestrel[which(trimestrel>12)]-12
  path_save <- paste0(dir_save,"/",month.abb[ic],"_",paste(month.abb[trimestrel],collapse = "-"),".tsv.gz")
  download.file(route,path_save)
  gunzip(path_save)
  
  return(paste("Successful download",path_save))
}


#Ejemplo para descargar Feb_Abr-May-Jun (trimestre AMJ con condición inicial 
#en Febrero para el área xmin =0, xmax =359, ymin= -30, ymax =30)

area1 <- c(0,359,-30,30) # Modifique esta línea con el Área a descargar xmin, xmax, ymin, ymax)

i_month <- 3 # Modifique esta línea con el Mes de inicio del trimestre de interés
lg <- 3 # Modifique esta línea con el Tamaño de la temporada a descargar. Por ej. AMJ lg=3/ AM lg=2
ic <- 2 # Modifique esta línea con la Condición inicial de la corrida del pronóstico o también conocido como lead time
firs_year <- 1981 # Modifique esta línea con el Año de inicio de la descarga
last_year <- 2019 # Modifique esta línea con el Año de final de la descarga

# Aquí se ejecuta la función con los parámetros dados anteriormente
download_CFSV2_CPT_1(firs_year,last_year,i_month,ic,dir_save,area1,lg=lg)


# Función para descargar dos áreas de la TSM del modelo CFSv2 --------------
download_CFSV2_CPT_2=function(firs_year,last_year,i_month,ic,dir_save,area1,area2){
  
  lg_s <-lg-1
  lead <- i_month-ic
  if(lead<0)lead <- lead + 12
  route <- paste0("http://iridl.ldeo.columbia.edu/SOURCES/.NOAA/.NCEP/.EMC/.CFSv2/.ENSEMBLE/.OCNF/.surface/.TMP/SOURCES/.NOAA/.NCEP/.EMC/.CFSv2/.REALTIME_ENSEMBLE/.OCNF/.surface/.TMP/appendstream/350/maskge/S/%280000%201%20",month.abb[ic],"%20",firs_year,"-",last_year,"%29/VALUES/L/",lead,".5/",lead+lg_s,".5/RANGE/%5BL%5D//keepgrids/average/M/1/24/RANGE/%5BM%5Daverage/X/",area1[1],"/",area1[2],"/flagrange/Y/",area1[3],"/",area1[4],"/flagrange/add/1/flaggt/X/",area2[1],"/",area2[2],"/flagrange/Y/",area2[3],"/",area2[4],"/flagrange/add/1/flaggt/add/mul/0/setmissing_value/-999/replaceNaN/%5BX/Y%5D%5BS/L/add/%5Dcptv10.tsv.gz")
  trimestrel <- (ic+lead):(ic+lead+lg_s)
  if(sum(trimestrel>12)>0)trimestrel[which(trimestrel>12)]=trimestrel[which(trimestrel>12)]-12
  path_save <- paste0(dir_save,"/",month.abb[ic],"_",paste(month.abb[trimestrel],collapse = "-"),".tsv.gz")
  download.file(route,path_save)
  gunzip(path_save)
  
  return("Successful download")
  
}

#Ejemplo para descargar Feb_Abr-May-Jun (trimestre AMJ con condición inicial 
#en Febrero para el área xmin =0, xmax =359, ymin= -30, ymax =30)
area1 <- c(180, 270,-13,13) # Modificar
area2 <- c(45 , 105,-13,13) # Modificar

lg <- 2 # Modificar
firs_year <- 1981 # Modificar
last_year <- 2016 # Modificar
i_month <- 3 # Modificar
ic <- 2 # Modificar

# Aquí se ejecuta la función con los parámetros dados anteriormente para 2 áreas
download_CFSV2_CPT_2(firs_year,last_year,i_month,ic,dir_save,area1,area2)

