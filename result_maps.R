

suppressMessages(if(!require(rworldmap)){install.packages('rworldmap'); library(rworldmap)} else {library(rworldmap)})
suppressMessages(if(!require(raster)){install.packages('raster'); library(raster)} else {library(raster)})
suppressMessages(if(!require(ggplot2)){install.packages('ggplot2'); library(ggplot2)} else {library(ggplot2)})
suppressMessages(if(!require(rasterVis)){install.packages('rasterVis'); library(rasterVis)} else {library(rasterVis)})
suppressMessages(if(!require(sf)){install.packages('sf'); library(sf)} else {library(sf)})
suppressMessages(if(!require(grid)){install.packages('grid'); library(grid)} else {library(grid)})
suppressMessages(if(!require(dplyr)){install.packages('dplyr'); library(dplyr)} else {library(dplyr)})
suppressMessages(if(!require(tidyr)){install.packages('tidyr'); library(tidyr)} else {library(tidyr)})
suppressMessages(if(!require(rgeos)){install.packages('rgeos'); library(rgeos)} else {library(rgeos)})

path_metric <-  "D:/OneDrive - CGIAR/Tobackup/CIAT/Projects/TNC-Honduras/zone/output/all_domain"
path_output <-  "D:/OneDrive - CGIAR/Tobackup/CIAT/Projects/TNC-Honduras/zone/output/all_domain"
path_raw <- "D:/OneDrive - CGIAR/Tobackup/CIAT/Projects/TNC-Honduras/zone/output/raw_output/Aug_Dec-Jan-Feb_0"

cca_map <- function(path_raw , path_output) {
  
 
  xserie <- read.csv(paste0(path_raw, "_cca_scores_x.txt"),skip =2, header=T, sep="")
  yserie <- read.csv(paste0(path_raw,"_cca_scores_y.txt"),skip =2, header=T, sep="")
  
  xloadcca <-  read.csv(paste0(path_raw, "_cca_load_x.txt"),skip =2, header=T, sep="",check.names = F)
  yloadcca <-  read.csv(paste0(path_raw, "_cca_load_y.txt"),skip =2, header=T, sep="")
  
  xeigen<- read.csv(paste0(path_raw, "_pca_eigen_x.txt"),skip =2, header=T, sep="")
  yeigen <- read.csv(paste0(path_raw,"_pca_eigen_y.txt"),skip =2, header=T, sep="")
  
  
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
                                                                              legend.title = element_text(size = 12.5))  +
                                                                        scale_fill_gradientn(colours =myPalette(100), limits=c(-1,1))   
  
  yloadcca$val = yloadcca$X1/max(abs(yloadcca$X1),na.rm=T)
  y = yloadcca
  ext_y = c(min(as.numeric(y$Latitude))-0.5,max(as.numeric(y$Latitude))+0.5,min(as.numeric(y$Longitude))-0.5,max(as.numeric(y$Longitude))+0.5)
  
  sPDF <<- getMap()  
  new_map = fortify(sPDF)
  sel <<-new_map[new_map$lat>ext_y[1] & new_map$lat<ext_y[2] & new_map$long>ext_y[3] & new_map$long<ext_y[4],]
  
   
  p <- ggplot(sel, aes(x=long,y=lat)) # gráfique el pa�???s
  p <- p + geom_polygon(aes(fill=hole,group=group),fill="snow") + 
    scale_fill_manual(values=c("grey 80","grey 80"))  
   

  p <- p + geom_point(data=y, aes(x=Longitude, y=Latitude, col=val),size=1)
  p <- p + scale_color_gradientn(colours =myPalette(100), limits=c(-1,1)) + coord_equal()+ geom_path(aes(long,lat,group=group),color="black",size=0.3)
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
  
  
  datos_e <- data.frame(modes=xeigen$Mode, eigenx=xeigen$variance, eigeny=yeigen$variance, row.names = NULL)
  
  # gráfico de los modos 
  modosx <-  ggplot(datos_e, aes(modes,group = 1)) +   geom_line(aes(y = eigenx ),  colour="firebrick3" ) + geom_point(aes(y = eigenx ),  colour="firebrick3" ) +
     theme_bw() + theme( title =element_text(size=12, face='bold'),axis.text.y = element_text(size=12),  legend.position = "none", axis.text.x = element_text(angle = 0, hjust = 1, size = 12)) +
    guides(colour = guide_legend(title = " ")) + labs(x="Mode",y="% variance",title = "X Scree Plot") 
  
  
  modosy <-  ggplot(datos_e, aes(modes,group = 1)) +   geom_line(aes(y = eigeny ),  colour="firebrick3" ) + geom_point(aes(y = eigeny ),  colour="firebrick3" ) +
    theme_bw() + theme( title =element_text(size=12, face='bold'),axis.text.y = element_text(size=12),  legend.position = "none", axis.text.x = element_text(angle = 0, hjust = 1, size = 12)) +
    guides(colour = guide_legend(title = " ")) + labs(x="Mode",y="% variance",title = "Y Scree Plot") 
  
  layt<-grid.layout(nrow=1,ncol=3,widths=c(4/9,2.5/9, 2.5/9),default.units=c('null','null'))
 
  tiff(filename = paste0(path_output, "/cca_maps.tif"), width = 1700, height = 400,res=100,compression = 'lzw')
  grid.newpage()
  pushViewport(viewport(layout=layt))
  print(Map_x,vp=viewport(layout.pos.row=1,layout.pos.col=1))
  print(modos,vp=viewport(layout.pos.row=1,layout.pos.col=2))
  print(p,vp=viewport(layout.pos.row=1,layout.pos.col=3))
  dev.off()
  cat("Mapas CCA realizados...\n")
  
  layt<-grid.layout(nrow=1,ncol=2)
  
  tiff(filename = paste0(path_output, "/eigen_plot.tif"), width = 1500, height = 800,res=150,compression = 'lzw')
  grid.newpage()
  pushViewport(viewport(layout=layt))
  print(modosx,vp=viewport(layout.pos.row=1,layout.pos.col=1))
  print(modosy,vp=viewport(layout.pos.row=1,layout.pos.col=2))
  dev.off()
  cat("Scree plots realizados...\n")

}

cca_map(path_raw,path_output)

#### Indicadores
metric_map <- function(path_metric, path_output){
  
   
  all_ind_complete <-read.csv(paste0(path_metric, '/metrics.csv')) 
  
  trimesters <- unique(all_ind_complete$file)
  
  for(i in 1:length(trimesters)){
    
    
    
    all_ind <- all_ind_complete[all_ind_complete$file==trimesters[i],]
    
    myPalette <-  colorRampPalette(c("dodgerblue4", "dodgerblue1","deepskyblue","darkslategray1", "lightcyan",  "lemonchiffon1","khaki","sandybrown", "darkorange2","firebrick2"))
    
    ind <- ggplot(sel, aes(x=long,y=lat)) + 
      geom_polygon(aes(group=group),fill="snow") + 
      geom_point(data=all_ind, aes(x= longitud, y= latitud, group= 1 ,col=kendall),size=3) + 
      geom_path(aes(long,lat,group=group,fill=hole),color="black",size=0.3)  + 
      theme_bw() + scale_color_gradientn(colours =myPalette(100), limits=c(0,100)) + coord_equal()+
      labs(title="2AFC Score", x=" ", y=" ", col=" ")+theme( legend.position = "bottom",legend.key.width  = unit(1.5, "cm"))
    

    ind_1 <- ggplot(sel, aes(x=long,y=lat)) + 
      geom_polygon(aes(fill=hole,group=group),fill="snow") +  
      geom_point(data=all_ind, aes(x= longitud, y= latitud, group= 1 ,col=pearson),size=3) + 
      geom_path(aes(long,lat,group=group),color="black",size=0.3)  + 
      theme_bw() + scale_color_gradientn(colours =myPalette(100), limits=c(-1,1)) + coord_equal() + 
      labs(title="Pearson's Correlation", x=" ", y=" ", col=" ")+ theme(legend.position = "bottom",legend.key.width  = unit(1.5, "cm"))
    
    
    ind_2 <-ggplot(sel, aes(x=long,y=lat)) + 
      geom_polygon(aes(fill=hole,group=group),fill="snow") +  
      geom_point(data=all_ind, aes(x= longitud, y= latitud, group= 1 ,col=roc_b),size=3) + 
      geom_path(aes(long,lat,group=group),color="black",size=0.3)  + 
      theme_bw() + scale_color_gradientn(colours =myPalette(100), limits=c(0,1)) + coord_equal() + 
      labs(title="ROC Area (Below-Normal)", x=" ", y=" ", col=" ")+ theme(legend.position = "bottom",legend.key.width  = unit(1.5, "cm"))
    
    
    ind_3 <-ggplot(sel, aes(x=long,y=lat)) + 
      geom_polygon(aes(fill=hole,group=group),fill="snow") +  
      geom_point(data=all_ind, aes(x= longitud, y= latitud, group= 1 ,col=roc_a),size=3) + 
      geom_path(aes(long,lat,group=group),color="black",size=0.3)  + 
      theme_bw() + scale_color_gradientn(colours =myPalette(100), limits=c(0,1)) + coord_equal() + 
      labs(title="ROC Area (Above-Normal)", x=" ", y=" ", col=" ")+ theme(legend.position = "bottom",legend.key.width  = unit(1.5, "cm"))
    
    ind_4 <-ggplot(sel, aes(x=long,y=lat)) + 
      geom_polygon(aes(fill=hole,group=group),fill="snow") +  
      geom_point(data=all_ind, aes(x= longitud, y= latitud, group= 1 ,col=hit_s),size=3) + 
      geom_path(aes(long,lat,group=group),color="black",size=0.3)  + 
      theme_bw() + scale_color_gradientn(colours =myPalette(100), limits=c(0,100)) + coord_equal() + 
      labs(title="Hit Score", x=" ", y=" ", col=" ")+ theme(legend.position = "bottom",legend.key.width  = unit(1.5, "cm"))
    
    
    ind_5 <-ggplot(sel, aes(x=long,y=lat)) + 
      geom_polygon(aes(fill=hole,group=group),fill="snow") +  
      geom_point(data=all_ind, aes(x= longitud, y= latitud, group= 1 ,col=hit_ss),size=3) + 
      geom_path(aes(long,lat,group=group),color="black",size=0.3)  + 
      theme_bw() + scale_color_gradientn(colours =myPalette(100), limits=c(-100,100)) + coord_equal() + 
      labs(title="Hit Skill Score", x=" ", y=" ", col=" ")+ theme(legend.position = "bottom",legend.key.width  = unit(1.5, "cm"))
    

    layt<-grid.layout(nrow=2,ncol=1)
    
    tiff(filename = paste0(path_output,"/" ,trimesters[i], "_metrics_maps.tif"), width = 800, height = 800,res=130,compression = 'lzw')
    grid.newpage()
    pushViewport(viewport(layout=layt))
    print(ind,vp=viewport(layout.pos.row=1,layout.pos.col=1))
    print(ind_1,vp=viewport(layout.pos.row=2,layout.pos.col=1))
        
    dev.off()
    cat("Mapas Metricas realizados...\n")
    
      
    layt<-grid.layout(nrow=2,ncol=1)
    
    tiff(filename = paste0(path_output,"/" ,trimesters[i], "_roc_maps.tif"), width = 800, height = 800,res=130,compression = 'lzw')
    grid.newpage()
    pushViewport(viewport(layout=layt))
    print(ind_2,vp=viewport(layout.pos.row=1,layout.pos.col=1))
    print(ind_3,vp=viewport(layout.pos.row=2,layout.pos.col=1))
    
    dev.off()
    cat("Mapas ROC realizados...\n")
    
    layt<-grid.layout(nrow=2,ncol=1)
    
    tiff(filename = paste0(path_output,"/" ,trimesters[i], "_hit_maps.tif"), width = 800, height = 800,res=130,compression = 'lzw')
    grid.newpage()
    pushViewport(viewport(layout=layt))
    print(ind_4,vp=viewport(layout.pos.row=1,layout.pos.col=1))
    print(ind_5,vp=viewport(layout.pos.row=2,layout.pos.col=1))
    
    dev.off()
    cat("Mapas Hit realizados...\n")
    
    max_C<-apply(all_ind[,10:12], 1, max)
    cat<-ifelse(all_ind[,10]==max_C, "Below", ifelse(all_ind[,11]==max_C, "Normal", ifelse(all_ind[,12]==max_C, "Above",0)))
    maximos<-cbind.data.frame(all_ind$id, all_ind$longitud, all_ind$latitud, max_C, cat)
    
    # Aqui se ingresan los datos de las estaciones
    maximos$cat = factor(maximos$cat, levels = c("Above", "Normal", "Below"))
    p <- ggplot(sel, aes(x=long,y=lat)) + 
      geom_polygon(aes(fill=hole,group=group),fill="snow")
    
    maxi <- p + geom_point(data=maximos, aes(x=all_ind$longitud, y=all_ind$latitud,size=max_C, colour=cat))+ scale_size_continuous(name = " ",
                                                                                                                                   breaks = seq(25,100,25),
                                                                                                                                   limits = c(0, 100),
                                                                                                                                   labels = c("0-25", "25-50", "50-75","75-100"),
                                                                                                                                   range = c(0, 5) )+
      geom_path(aes(long,lat,group=group),color="black",size=0.3) + scale_colour_manual(values = c("steelblue3","lightgreen","tomato2")) +
      coord_equal() + theme( legend.key.height=unit(1,"cm"),legend.key.width=unit(0.5,"cm"),
                             legend.text=element_text(size=8),
                             panel.background=element_rect(fill="white",colour="black"),
                             axis.text=element_text(colour="black",size=10),
                             axis.title=element_text(colour="black",size=10,face="bold"),
                             #legend.position = "bottom", 
                             legend.title=element_blank())  + labs(title="Probabilistic Forecast")+labs( x=" ", y=" ", size=" ")
    
    seq(0,100,25)  
    seq(25,75,25)
    c("0-25", "25-50", "50-75","75-100")
    tiff(paste0(path_output,"/" ,trimesters[i], "_probabilistic_maps.tif"),width = 3000, height = 3000, units = "px", res = 400,compression = 'lzw')
    print(maxi)
    dev.off()  
    cat("Mapas Probabilistico realizados...\n")
    
    
    myPalette2 = colorRampPalette(c("steelblue2" , "deepskyblue", "lightcyan","khaki","yellow", "orange1","darkorange3", "red", "firebrick4"))
    
    above <-ggplot(sel, aes(x=long,y=lat)) + 
      geom_polygon(aes(fill=hole,group=group),fill="snow") +  
      geom_point(data=all_ind, aes(x= longitud, y= latitud, group= 1 ,col=above),size=3) + 
      geom_path(aes(long,lat,group=group),color="black",size=0.3)  + 
      theme_bw() + scale_color_gradientn(colours =myPalette2(10), limits=c(0,100)) + coord_equal() + 
      labs(title="Above", x=" ", y=" ", col=" ")+ theme(legend.position = "bottom",legend.key.width  = unit(1.5, "cm"))
    
    normal <-ggplot(sel, aes(x=long,y=lat)) + 
      geom_polygon(aes(fill=hole,group=group),fill="snow") +  
      geom_point(data=all_ind, aes(x= longitud, y= latitud, group= 1 ,col=normal),size=3) + 
      geom_path(aes(long,lat,group=group),color="black",size=0.3)  + 
      theme_bw() + scale_color_gradientn(colours =myPalette2(10), limits=c(0,100)) + coord_equal() + 
      labs(title="Normal", x=" ", y=" ", col=" ")+ theme(legend.position = "bottom",legend.key.width  = unit(1.5, "cm"))
    
    below <-ggplot(sel, aes(x=long,y=lat)) + 
      geom_polygon(aes(fill=hole,group=group),fill="snow") +  
      geom_point(data=all_ind, aes(x= longitud, y= latitud, group= 1 ,col=below),size=3) + 
      geom_path(aes(long,lat,group=group),color="black",size=0.3)  + 
      theme_bw() + scale_color_gradientn(colours =myPalette2(10), limits=c(0,100)) + coord_equal() + 
      labs(title="Below", x=" ", y=" ", col=" ")+ theme(legend.position = "bottom",legend.key.width  = unit(1.5, "cm"))
    
    layt<-grid.layout(nrow=1,ncol=3)
    
    tiff(filename = paste0(path_output,"/" ,trimesters[i], "_probabilities_maps.tif"), width = 1400, height = 800,res=130,compression = 'lzw')
    grid.newpage()
    pushViewport(viewport(layout=layt))
    print(below,vp=viewport(layout.pos.row=1,layout.pos.col=1))
    print(normal,vp=viewport(layout.pos.row=1,layout.pos.col=2))
    print(above,vp=viewport(layout.pos.row=1,layout.pos.col=3))
    
    dev.off()
    cat("Mapas Probabilidades realizados...\n")
  }
 

}

metric_map(path_metric, path_output)
