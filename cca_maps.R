path <- 'D:/OneDrive - CGIAR/Desktop/test/cpt_r/corrida/output/complete/' # Main folder


suppressMessages(if(!require(rworldmap)){install.packages('rworldmap'); library(rworldmap)} else {library(rworldmap)})
suppressMessages(if(!require(raster)){install.packages('raster'); library(raster)} else {library(raster)})
suppressMessages(if(!require(ggplot2)){install.packages('ggplot2'); library(ggplot2)} else {library(ggplot2)})
suppressMessages(if(!require(rasterVis)){install.packages('rasterVis'); library(rasterVis)} else {library(rasterVis)})
suppressMessages(if(!require(sf)){install.packages('sf'); library(sf)} else {library(sf)})
suppressMessages(if(!require(grid)){install.packages('grid'); library(grid)} else {library(grid)})
suppressMessages(if(!require(dplyr)){install.packages('dplyr'); library(dplyr)} else {library(dplyr)})
suppressMessages(if(!require(tidyr)){install.packages('tidyr'); library(tidyr)} else {library(tidyr)})
suppressMessages(if(!require(rgeos)){install.packages('rgeos'); library(rgeos)} else {library(rgeos)})


cca_map <- function(path = path) {
  
  xserie <- read.csv(paste0(dir.O.CPT, "/Feb_Apr-May-Jun_0_cca_serie_x.txt"),skip =2, header=T, sep="")
  yserie <- read.csv(paste0(dir.O.CPT,"/Feb_Apr-May-Jun_0_cca_serie_y.txt"),skip =2, header=T, sep="")
  
  xloadcca <-  read.csv(paste0(dir.O.CPT, "/Feb_Apr-May-Jun_0_cca_load_x.txt"),skip =2, header=T, sep="",check.names = F)
  yloadcca <-  read.csv(paste0(dir.O.CPT, "/Feb_Apr-May-Jun_0_cca_load_y.txt"),skip =2, header=T, sep="")
  
  
  
  ## =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= ##
  ## CCA Maps
  ## =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= ##
  x = xloadcca
  
  ext = c(min(as.numeric(rownames(x))),max(as.numeric(rownames(x))),min(as.numeric(colnames(x))),max(as.numeric(colnames(x))))
  mapa_base=raster(nrow=dim(x)[1],ncol=dim(x)[2])
  extent(mapa_base) <- extent(ext[3],ext[4],ext[1],ext[2])
  val=c(as.matrix(t(x),ncol=1,byrow = T))
  val=as.numeric(val)
  val[val==-999.000]=NA
  val = val/max(abs(val),na.rm=T)
  values(mapa_base)=val
  
  
  #myPalette <-  colorRampPalette(c("navyblue","#2166AC", "dodgerblue3","lightblue", "lightcyan",  "white",  "yellow","orange", "orangered","#B2182B", "red4"))
  myPalette <-  colorRampPalette(c("dodgerblue4", "dodgerblue1","deepskyblue","darkslategray1", "lightcyan",  "white",  "lemonchiffon1","khaki","sandybrown", "darkorange2","firebrick2"))
  
  Map_x<- rasterVis::gplot(mapa_base) + geom_tile(aes(fill = value)) + coord_equal() + 
    labs(title="X Spatial Loadings (Mode 1)",x="",y=" ", fill = " ")  + theme(legend.key.height=unit(0.5,"cm"),legend.key.width=unit(2,"cm"),
                                                                              legend.text=element_text(size=10),
                                                                              panel.background=element_rect(fill="white",colour="black"),
                                                                              axis.text=element_text(colour="black",size=12),
                                                                              axis.title=element_text(colour="black",size=12,face="bold"),
                                                                              legend.position = "bottom", 
                                                                              legend.title = element_text(size = 12.5))  +
    scale_fill_gradientn(colours =myPalette(100), limits=c(-1,1))   
  
  yloadcca$val = yloadcca$X1/max(abs(yloadcca$X1),na.rm=T)
  y = yloadcca
  ext_y = c(min(as.numeric(y$Latitude))-0.5,max(as.numeric(y$Latitude))+0.5,min(as.numeric(y$Longitude))-0.5,max(as.numeric(y$Longitude))+0.5)
  
  sPDF <- getMap()  
  new_map = fortify(sPDF)
  sel = new_map[new_map$lat>ext_y[1] & new_map$lat<ext_y[2] & new_map$long>ext_y[3] & new_map$long<ext_y[4],]
  
   
  # Realice el gráfico de las correlaciones entre las estaciones y el modo 1 de y 
  # IT <- getData(name = "GADM", country = "Honduras", level = 1)
  p <- ggplot(sel, aes(x=long,y=lat)) # gráfique el pa�???s
  p <- p + geom_polygon(aes(fill=hole,group=group),fill="snow") + 
    scale_fill_manual(values=c("grey 80","grey 80")) + 
    geom_path(aes(long,lat,group=group),color="black",size=0.3)

  p <- p + geom_point(data=y, aes(x=Longitude, y=Latitude, col=val),size=1)
  p <- p + scale_color_gradientn(colours =myPalette(100), limits=c(-1,1)) + coord_equal()
  p <-  p + theme(legend.key.height=unit(1,"cm"),legend.key.width=unit(0.5,"cm"),
                 legend.text=element_text(size=12),
                 panel.background=element_rect(fill="white",colour="black"),
                 axis.text=element_text(colour="black",size=12),
                 axis.title=element_text(colour="black",size=12,face="bold"),
                 #legend.position = "bottom", 
                 legend.title = element_text(size = 12)) + labs(title="Y Spatial Loadings (Mode 1)",  x= "", y= "", colour="") 
  
  
  ## Gráficos Componentes 
  
  # Se crea una trama de datos con la fecha y las componentes 
  datos <- data.frame(X=xserie$X1, Y=yserie$X1, row.names = NULL)
  datos$X = round(datos$X ,4) # redondee los modos 
  datos$Y = round(datos$Y ,4) # redondee los modos 
  datos$date = as.numeric(substring(rownames(xserie),1,4))
  
  datos[datos$X==-999.0000,1:2]=0
  datos[datos$Y==-999.0000,1:2]=0 # quite los valore NA
  # quite los valore NA
  datos[,1:2]=datos[,1:2]*100 # multipliquelos * 100
  
  # gráfico de los modos 
  modos <-  ggplot(datos, aes(date,group = 1)) +   geom_line(aes(y = X ),  colour="firebrick3" ) + 
    geom_line(aes(y = Y),  colour="forestgreen")  + 
    geom_hline(yintercept = 0, colour="gray") + theme_bw() + 
    theme( title =element_text(size=12, face='bold'),axis.text.y = element_text(size=12),  legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1, size = 12)) +
    guides(colour = guide_legend(title = " ")) + labs(subtitle=paste("Canonical Correlation = ", round(cor(datos$X,datos$Y),3),sep = ""),x="",y="Scores (X red; Y green) (*100)",
                                                      title = "Temporal Scores (Mode 1)") + scale_x_continuous(breaks = seq(min(datos$date),max(datos$date),3))
  
  layt<-grid.layout(nrow=1,ncol=3,widths=c(4/9,2.5/9, 2.5/9),default.units=c('null','null'))
 
  tiff(filename = paste0(path, "cca_maps.tif"), width = 1700, height = 400,res=100,compression = 'lzw')
  grid.newpage()
  pushViewport(viewport(layout=layt))
  print(Map_x,vp=viewport(layout.pos.row=1,layout.pos.col=1))
  print(modos,vp=viewport(layout.pos.row=1,layout.pos.col=2))
  print(p,vp=viewport(layout.pos.row=1,layout.pos.col=3))
  dev.off()
  cat("Mapas CCA realizados...")
  
}

