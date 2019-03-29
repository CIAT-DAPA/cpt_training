# Creado por: Lizeth Llanos
# Ejercicios para manejar archivos raster y shapefile

# Cargar librerias. Recuerde instalarlas previamente con install.packages ---
# Link para descargar datos de chirps para México: https://iridl.ldeo.columbia.edu/SOURCES/.UCSB/.CHIRPS/.v2p0/.monthly/.global/.precipitation/Y/%2812%29/%2835%29/RANGEEDGES/X/%28-120%29/%28-80%29/RANGEEDGES/datafiles.html

library(raster)
library(ncdf4)
library(maptools)
library(rgdal)


# Cargar stack de chirps --------------------------------------------------


chirps_all = stack("D:/OneDrive - CGIAR/Tobackup/CIAT/Projects/Capacitaciones/México/Ejercicio CHIRPS/mexico.nc")
chirps_all

plot(chirps_all[[1]])
#plot(chirps_all)

avg_all = mean(chirps_all)
plot(avg_all)

# Descargar shapefile hnd------------------------------------------------------


shp_mex <- getData(name = "GADM", country = "Mexico", level = 1)
#shp_mex <- shapefile("C:/Users/lllanos/Google Drive/PNUD_HND/Taller_1/Jornada_1/ejercicios_R/regiones_shapefile/Regiones_Desarrollo_prj_v2.shp")

plot(chirps_all[[1]])
plot(shp_mex, add=T)

# Ajustar mapa a un shapefile ---------------------------------------------


chirps_mex_m = crop(x = chirps_all, y = extent(shp_mex))
plot(chirps_mex_m[[1]])

chirps_mex = mask(chirps_mex_m,shp_mex)

plot(chirps_mex[[1]])
plot(shp_mex, add =T)


# Filtrar shapefile -------------------------------------------------------

View(shp_mex@data)
shp_mex_1 = shp_mex[shp_mex$ID_1==1 | shp_mex$ID_1==2, ]
plot(chirps_all[[1]])
plot(shp_mex_1, add=T)


chirps_mex_m1 = crop(x = chirps_all, y =extent(shp_mex_1) )
chirps_mex_1 = mask(chirps_mex_m1,shp_mex_1)
plot(chirps_mex_1[[1]])

# Exportar nuevo raster ---------------------------------------------------

writeRaster(chirps_mex_1,"mexico_chirps.tiff")

# Exportar datos del raster -----------------------------------------------

data_mex = as.data.frame(t(rasterToPoints(chirps_mex_1)))
data_mex$dates = c("x","y",substring(as.Date(seq(as.Date("1981/01/01"),as.Date("2019/02/28"),"month"),format="%Y-%m"),1,7))
data_mex = data_mex[,c(ncol(data_mex),1:(ncol(data_mex)-1))]
write.csv(data_mex,"mexico_chirps_data.csv",quote = F,row.names = F)

