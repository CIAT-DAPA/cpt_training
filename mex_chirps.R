# Creado por: Lizeth Llanos y Diego Agudelo
# Ejercicios para manejar archivos raster y shapefile

# Cargar librerias. Recuerde instalarlas previamente con install.packages ---
# Link para descargar datos de chirps : https://iridl.ldeo.columbia.edu/SOURCES/.UCSB/.CHIRPS/.v2p0/.monthly/.global/.precipitation/Y/%2812%29/%2835%29/RANGEEDGES/X/%28-120%29/%28-80%29/RANGEEDGES/datafiles.html
suppressMessages(if(!require(raster)){install.packages('raster'); library(raster)} else {library(raster)})
suppressMessages(if(!require(ncdf4)){install.packages('ncdf4'); library(ncdf4)} else {library(ncdf4)})
suppressMessages(if(!require(maptools)){install.packages('maptools'); library(maptools)} else {library(maptools)})
suppressMessages(if(!require(rgdal)){install.packages('rgdal'); library(rgdal)} else {library(rgdal)})

# establecer directorio de trabajo

setwd("E:/")

# Cargar stack de chirps --------------------------------------------------


chirps_all = stack("E:/boyaca.nc")
chirps_all

plot(chirps_all[[1]])
#plot(chirps_all)

avg_all = mean(chirps_all)
plot(avg_all)

# Descargar shapefile hnd------------------------------------------------------

shp_col <- getData(name = "GADM", country = "Colombia", level = 1)
#shp_mex <- shapefile("C:/Users/lllanos/Google Drive/PNUD_HND/Taller_1/Jornada_1/ejercicios_R/regiones_shapefile/Regiones_Desarrollo_prj_v2.shp")

plot(chirps_all[[1]])
plot(shp_col, add=T)

# Ajustar mapa a un shapefile ---------------------------------------------

chirps_col = crop(x = chirps_all, y = extent(shp_col))
plot(chirps_col[[1]])

chirps_col = mask(chirps_col,shp_col)

plot(chirps_mex[[1]])
plot(shp_col, add =T)

# Filtrar shapefile -------------------------------------------------------

View(shp_col@data)
shp_boy = shp_col[shp_col$NAME_1=="Boyacá", ]
plot(chirps_all[[1]])
plot(shp_boy, add=T)


chirps_boy = crop(x = chirps_all, y =extent(shp_boy) )
chirps_boy = mask(chirps_boy,shp_boy)
plot(chirps_boy[[1]])
plot(shp_boy, add=T)

# Exportar nuevo raster ---------------------------------------------------

writeRaster(chirps_boy,"boyacá_chirps.tiff")

# Exportar datos del raster -----------------------------------------------

data_boy = as.data.frame(t(rasterToPoints(chirps_boy)))
data_boy$dates = c("x","y",substring(as.Date(seq(as.Date("1981/01/01"),as.Date("2020/06/28"),"month"),format="%Y-%m"),1,7))
data_boy = data_boy[,c(ncol(data_boy),1:(ncol(data_boy)-1))]
write.csv(data_boy,"boyacá_chirps_data.csv",quote = F,row.names = F)

