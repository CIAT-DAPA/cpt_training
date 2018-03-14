# Script preliminar CPT graphs 
# Created by: Alejandra Esquivel (a.esquivel@cgiar.org)
# Date: March 2018


## =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= ##
## Packages
## =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= ##
suppressMessages(if(!require(raster)){install.packages('raster'); library(raster)} else {library(raster)})
suppressMessages(if(!require(ggplot2)){install.packages('ggplot2'); library(ggplot2)} else {library(ggplot2)})
suppressMessages(if(!require(rasterVis)){install.packages('rasterVis'); library(rasterVis)} else {library(rasterVis)})
suppressMessages(if(!require(sf)){install.packages('sf'); library(sf)} else {library(sf)})
suppressMessages(if(!require(grid)){install.packages('grid'); library(grid)} else {library(grid)})
suppressMessages(if(!require(dplyr)){install.packages('dplyr'); library(dplyr)} else {library(dplyr)})
suppressMessages(if(!require(tidyr)){install.packages('tidyr'); library(tidyr)} else {library(tidyr)})
suppressMessages(if(!require(tidyr)){install.packages('rgeos'); library(rgeos)} else {library(rgeos)})



path <- 'D:/graphs_data/' # Main folder


## =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= ##
## CCA Maps
## =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= ##
predictor <- "SST_CFSv2"


# files dor stations, SST, and CPT
dir.Stations <- paste0(path, 'Inputs/dep') 
dir.O.SST <-  paste0(path, 'Inputs/SST')  
dir.O.CPT <- paste0(path, 'Inputs/Cross_validated')  


# =-=-=-= final year 
final_year<-2013 # año final del periodo de entrenamiento. #(modificar)




### Transformar un archivo .tsv(CPT) en raster de CSFv2

transform_raster=function(x){
  mapa_base=raster()
  val=c(as.matrix(t(x),ncol=1,byrow = T))
  val=as.numeric(val)
  val[val==-999.000]=NA
  values(mapa_base)=val
  return(mapa_base)
}


rasterize=function(dates) { 
  
  if(require(raster)==FALSE){install.packages("raster")}
  library("raster")
  pos_years=!is.na(dates[1,])
  year_month=dates[1,][pos_years]
  if(substr(year_month[2],6,7)=="12"){year=as.numeric(substr(year_month[-1],1,4))+1
  }else{year=as.numeric(substr(year_month[-1],1,4))}
  total_row_delete=c(-1,-3,-(which(dates[,1]=="90.0")-2),-which(dates[,1]=="-90.0"),-which(dates[,1]==""))
  dates=dates[total_row_delete,-1]
  list_dates=split(dates,sort(rep(year,180)))
  all_raster=lapply(list_dates,transform_raster)
  layers=stack(all_raster)
  layers_crop=crop(layers,extent(-180, 180, -30, 30))
  
  return(layers_crop)
}


### Esta función convierte los datos de las estacines en datos trimestrales
data_trim=function(Estaciones_C, a){ #Los argumentos son el conjunto de las estaciones 
  ## y el mes de incio del periodo (a)
  stations=Estaciones_C 
  stations=stations[-1:-2,] # Quite las dos primeras filas (coordenadas)
  year=sort(rep(1981:2013,12)) # cree un vector para los años
  month=rep(1:12,length(1981:2013)) # cree el vector de meses
  data_station=cbind.data.frame(year,month,stations, row.names=NULL) # cree un data frame 
  pos=seq(a,dim(data_station)[1],12) #  posiciones ne las que se encuentra el mes a
  pos_select=sort(c(pos,pos+1,pos+2)) # muestre las posiciones del trimestre
  # Agregue los datos del trimestre y luego sumelos. 
  data_out=aggregate(data_station[pos_select,-1:-2],by=list(sort(rep(1:(length(pos)),3))),sum)
  data_out_final=na.omit(data_out[,-1]) # Elimine los NA
  years_y=na.omit(year[pos+1]) # Elimine los NA
  data_out_final=data.frame(years_y,data_out_final) # Cree un data frame con los datos finales
  return(data_out_final)
} # devuelva los datos finales





# =-=-=-=-=- read files 

xserie <- read.csv(paste0(dir.O.CPT, "/X_CCA_Map_Series.txt"),skip =2, header=T, sep="")
yserie <- read.csv(paste0(dir.O.CPT,"/Y_CCA_Map_Series.txt"),skip =2, header=T, sep="")

SST<-read.table(paste(dir.O.SST,"/DEF_Nov.tsv",sep=""),sep="\t",dec=".",skip =2,fill=TRUE,na.strings =-999)
## Conversión a raster
SST<-rasterize(SST)

var_ocanoAt <-SST[[1:31]]
b<-extent(-180, 180, -30, 30)
var_ocanoAt=crop(var_ocanoAt, b)

# estaciones
Estaciones_C <- read.delim(paste0(dir.Stations,"/precip_valle.txt"),skip =3, header=T, sep="")
a <- 12 # mes de inicio del trimestre


ruta <- 'outputs/' # modificar
lead <- 'DEF_Nov'



### Lectura del shp
colombia <- sf::st_read(dsn = paste0(path, 'Inputs/shp/colombia_depts.shp')) %>%
  as('Spatial') %>% crop(., extent( -77.5, -75.6, 3, 5 ))

plot(colombia)


# var_oceanoAt= variable oceano atmosferica
# y serie = el modo en y
# Estaciones_C= archivo de estaciones en el cual se realizo CPT
# xserie = el modo en x
# a = mes de inicio del trimestre de las estaciones
cca_maps<-function(var_ocanoAt, yserie, Estaciones_C, xserie, a){
  
  ocean=which(!is.na(var_ocanoAt[[1]][])) # tome las posiciones en las que la variable sea diferente de NA
  correl=array(NA,length(ocean)) # relice un arreglo del tamaño de oceano 
  var_table=var_ocanoAt[] # Realice una tabla de la variable
  
  for(i in 1:length(ocean)){ # En todos los pixeles diferentes de NA
    var_pixel=var_table[ocean[i],] # Extraiga el pixel i 
    correl[i]=cor(xserie$X1,var_pixel) # realice la correlación entre el pixel i el modo 1 de x
  } # 
  
  correl_map=var_ocanoAt[[1]] # Cree un raster vacio 
  correl_map[]=NA 
  correl_map[ocean]=correl # Almacene en el raster los NA 
  
  
  extent(correl_map) <- extent(0,360,-30,30)
  
  myPalette <-  colorRampPalette(c("navyblue","#2166AC", "dodgerblue3","lightblue", "lightcyan",  "white",  "yellow","orange", "orangered","#B2182B", "red4"))
  ewbrks <- c(seq(0,180,45), seq(225, 360, 45))
  nsbrks <- seq(-30,30,15)
  ewlbls <- unlist(lapply(ewbrks, function(x) ifelse(x <= 180, paste(abs(x), "°E"), ifelse(x > 180, paste( abs(360-x), "°W"),x))))
  nslbls <- unlist(lapply(nsbrks, function(x) ifelse(x < 0, paste(abs(x), "°S"), ifelse(x > 0, paste(abs(x), "°N"),x))))
  
  
  # Realice el mapa de Correlaciones entre la variable y el modo 1 de x
  Map_x<- rasterVis::gplot(correl_map) + geom_tile(aes(fill = value)) + coord_equal() + 
    labs(title="X Spatial Loadings (Mode 1)",x="",y=" ", fill = " ")  + theme(legend.key.height=unit(0.5,"cm"),legend.key.width=unit(2,"cm"),
                                                                              legend.text=element_text(size=12),
                                                                              panel.background=element_rect(fill="white",colour="black"),
                                                                              axis.text=element_text(colour="black",size=12),
                                                                              axis.title=element_text(colour="black",size=12,face="bold"),
                                                                              legend.position = "bottom", 
                                                                              legend.title = element_text(size = 12.5))  +
    scale_fill_gradientn(colours =myPalette(100), limits=c(-1,1))  + 
    scale_x_continuous(breaks = ewbrks, labels = ewlbls, expand = c(0, 0)) +
    scale_y_continuous(breaks = nsbrks, labels = nslbls, expand = c(0, 0))
  #scale_fill_gradient2(low="#2166AC",mid = "white", high="#B2182B",name = " ",  limits=c(-1,1)) 
  
  
  ###### Graficos de y
  # Convierta los datos de las estaciones en trimestrales 
  data<-data_trim(Estaciones_C, a)
  
  # La organización de la información se hace de acuerdo al mes de estudio.
  if(a == 12 | lead== "MAM_Nov"| lead== "MAM_Sep"| lead== "JJA_Dec"){ 
    data<-data[data$years_y!=1981 & data$years_y!=1982,]
  } else  data<-data[data$years_y!=1981,]
  
  correl_y=0 # inicialice las correlaciones con x
  for(i in 2:length(data)){ # realice las correlaciones para todas las estaciones
    correl_y[i-1]<-cor(data[,i],yserie$X1) # correlaciones entre la estación i y el modo 1 de x
  }
  
  Estacion=names(Estaciones_C) # extraiga los nombres de las estaciones
  coor<-data.frame(t(Estaciones_C[1:2,]), row.names = NULL) # extraiga las coordenadas
  # Cree un data frame con la información de las estaciones y las correlaciones
  datos2<-data.frame(Estacion,Long=coor$cpt.X, Lat=coor$cpt.Y,  Correly=correl_y, row.names = NULL)
  datos2$Correly=round(datos2$Correly ,3) # redondee el valor de las correlaciones a tres cifras
  
  
  # geom_polygon(data = shp_colombia, aes(x=long, y = lat, group = group), color = "black", fill = "white")
  
  
  # Realice el gráfico de las correlaciones entre las estaciones y el modo 1 de y 
  p <- ggplot(colombia, aes(x=long,y=lat)) # gráfique el país
  p <- p + geom_polygon(aes(fill=hole,group=group),fill="snow") + 
    scale_fill_manual(values=c("grey 80","grey 80")) + 
    geom_path(aes(long,lat,group=group,fill=hole),color="black",size=0.3)
  
  
  ewbrks <-  seq (round(min(datos2$Long), 1) ,  round(max(datos2$Long), 1), 0.8)
  nsbrks <-  seq (round(min(datos2$Lat), 1) ,  round(max(datos2$Lat), 1), 0.8)
  ewlbls <- unlist(lapply(ewbrks, function(x)  paste(abs(x), "°W")))
  nslbls <- unlist(lapply(nsbrks, function(x) ifelse(x < 0, paste(abs(x), "°S"), ifelse(x > 0, paste(abs(x), "°N"),x))))
  
  
  # Aqui se ingresan los datos de las estaciones
  p <- p + geom_point(data=datos2, aes(x=Long, y=Lat, map_id=Estacion,col=Correly),size=3)
  p <- p + scale_color_gradientn(colours =myPalette(100), limits=c(-1,1)) + coord_equal()
  p<-  p + theme(legend.key.height=unit(1,"cm"),legend.key.width=unit(0.5,"cm"),
                 legend.text=element_text(size=12),
                 panel.background=element_rect(fill="white",colour="black"),
                 axis.text=element_text(colour="black",size=12),
                 axis.title=element_text(colour="black",size=12,face="bold"),
                 #legend.position = "bottom", 
                 legend.title = element_text(size = 12)) + labs(title="Y Spatial Loadings (Mode 1)",  x= "", y= "", colour="") +
    scale_x_continuous(breaks = ewbrks, labels = ewlbls, expand = c(0, 0)) +
    scale_y_continuous(breaks = nsbrks, labels = nslbls, expand = c(0, 0))
  
  
  ## Gráficos Componentes 
  
  # Se crea una trama de datos con la fecha y las componentes 
  datos<-data.frame(date=data$years_y, X=xserie$X1, Y=yserie$X1, row.names = NULL)
  datos$X=round(datos$X ,4) # redondee los modos 
  datos$Y=round(datos$Y ,4) # redondee los modos 
  datos[datos$X==-999.0000,2:3]=0 # quite los valore NA
  datos[,2:3]=datos[,2:3]*100 # multipliquelos * 100
  
  # gráfico de los modos 
  modos<-  ggplot(datos, aes(date)) +   geom_line(aes(y = X ),  colour="#B2182B" ) + 
    geom_line(aes(y = Y),  colour="chartreuse4")  + 
    geom_hline(yintercept = 0, colour="gray") + theme_bw() + 
    theme( title =element_text(size=12, face='bold'),axis.text.y = element_text(size=12),  legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1, size = 12)) +
    guides(colour = guide_legend(title = " ")) + labs(subtitle=paste("Canonical Correlation = ", round(cor(datos$X,datos$Y),3),sep = ""),x="",y="Scores (X red; Y green) (*100)",
                                                      title = "Temporal Scores (Mode 1)") 
  modos <- modos  +   scale_x_continuous(breaks = seq(1982,final_year,3))
  
  layt<-grid.layout(nrow=1,ncol=3,widths=c(4/9,2.5/9, 2.5/9),default.units=c('null','null'))
  #View the layout of plots
  #grid.show.layout(layt)
  
  
  # 5 * 11
  png(filename = paste0(path, ruta, "cca_maps.png"), width = 1700, height = 400,res=100)
  grid.newpage()
  pushViewport(viewport(layout=layt))
  print(Map_x,vp=viewport(layout.pos.row=1,layout.pos.col=1))
  print(modos,vp=viewport(layout.pos.row=1,layout.pos.col=2))
  print(p,vp=viewport(layout.pos.row=1,layout.pos.col=3))
  dev.off()
  
  
}

# run function 
cca_maps(var_ocanoAt, yserie, Estaciones_C, xserie, a)





#### Indicadores

all_ind <-read.csv(paste0(path, 'metrics.csv')) 


ind <- all_ind %>% 
  select(- file, -goodness) %>% 
  gather(ind, value,   pearson:kendall)




myPalette <-  colorRampPalette(c("navyblue","#2166AC", "dodgerblue3","lightblue", "lightcyan",  "white",  "yellow","orange", "orangered","#B2182B", "red4"))

### Shp Perú
shp <- sf::st_read(dsn = paste0(path, 'Inputs/shp/PERU_DEPART/DEPARTAMENTO.shp')) %>%
  as('Spatial') 


# si tentemos varios indicadores con la misma escala de medida 
ggplot(shp, aes(x=long,y=lat)) + 
  geom_polygon(aes(fill=hole,group=group),fill="snow") + 
  scale_fill_manual(values=c("grey 80","grey 80")) + 
  geom_path(aes(long,lat,group=group,fill=hole),color="black",size=0.3)  + 
  geom_point(data=ind, aes(x= Longitud, y= Latitud, group= id ,col=value),
             size=3) + facet_wrap(~ind) +
  theme_bw() + scale_color_gradientn(colours =myPalette(100), limits=c(-1,1)) + coord_equal() 




# =-=-=-= si la escala de medidas del indicador es diferente 

ggplot(shp, aes(x=long,y=lat)) + 
  geom_polygon(aes(fill=hole,group=group),fill="snow") + 
  scale_fill_manual(values=c("grey 80","grey 80")) + 
  geom_path(aes(long,lat,group=group,fill=hole),color="black",size=0.3)  + 
  geom_point(data=all_ind, aes(x= Longitud, y= Latitud, group= id ,col=kendall),
             size=3)  +
  theme_bw() + scale_color_gradientn(colours =myPalette(100), limits=c(0,100)) + coord_equal() 





