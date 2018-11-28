# Script to automate runs from CPT(CFSV2) 
#  making predictor area selection 

# Created by: Diego Fernando Agudelo (d.agudelo@cgiar.org)
# Date: July 2018


############# parametros variables ###############

main_dir <- "D:/OneDrive - CGIAR/Desktop/Codigos_TNC/cpt_r"
modes_x <- 10
modes_y <- 10
modes_cca <- 5
trans <- 0       ###### 1 si quiere hacer transformacion y 0 si no quiere hacer transformacion
type_trans <- 2  ###### 1 transformacion normal y 2 transformacion gamma
  
########### Packages ###############

suppressMessages(if(!require(rworldmap)){install.packages('rworldmap'); library(rworldmap)} else {library(rworldmap)})
suppressMessages(if(!require(raster)){install.packages('raster'); library(raster)} else {library(raster)})
suppressMessages(if(!require(ggplot2)){install.packages('ggplot2'); library(ggplot2)} else {library(ggplot2)})
suppressMessages(if(!require(rasterVis)){install.packages('rasterVis'); library(rasterVis)} else {library(rasterVis)})
suppressMessages(if(!require(sf)){install.packages('sf'); library(sf)} else {library(sf)})
suppressMessages(if(!require(grid)){install.packages('grid'); library(grid)} else {library(grid)})
suppressMessages(if(!require(dplyr)){install.packages('dplyr'); library(dplyr)} else {library(dplyr)})
suppressMessages(if(!require(tidyr)){install.packages('tidyr'); library(tidyr)} else {library(tidyr)})
suppressMessages(if(!require(rgeos)){install.packages('rgeos'); library(rgeos)} else {library(rgeos)})
suppressMessages(if(require(stringr)==FALSE){install.packages("stringr",dependencies = TRUE)}) ;library("stringr")
suppressMessages(if(require(corpcor)==FALSE){install.packages("corpcor")}); library("corpcor")
suppressMessages(if(require(pcaPP)==FALSE){install.packages("pcaPP")}); library("pcaPP")
suppressMessages(if(require(RColorBrewer)==FALSE){install.packages("RColorBrewer")});library("RColorBrewer")

########### Functions ##############

transform_raster=function(x,y){
  mapa_base=raster(ext=y, res=c(1,1))
  val=c(as.matrix(t(x),ncol=1,byrow = T))
  val=as.numeric(val)
  val[val==-999|val== 0]=NA
  values(mapa_base)=val
  return(mapa_base)
}

data_raster=function(dates){
  
  year_month=dates[1,][!is.na(dates[1,])]
  year=ifelse(substr(year_month[-1],6,7)=="12",substr(year_month[-1],9,12),substr(year_month[-1],1,4))
  data_cpt1=na.omit(dates)
  pos=which(data_cpt1[,1]=="")
  pos=sort(rep(year,pos[2]-pos[1]))
  list_dates=split(data_cpt1,pos)
  lon=as.numeric(as.character(list_dates[[1]][1,-1]))
  lat=as.numeric(as.character(list_dates[[1]][-1,1]))
  cos_lat=diag(sqrt(cos((pi/180)*lat)))
  tables=lapply(list_dates,"[",-1,-1)
  tables_numeric=lapply(tables,function(x)sapply(x,function(y)as.numeric(as.character(y))))
  ex=extent(min(lon)-0.5,max(lon)+0.5,min(lat)-0.5,max(lat)+0.5)
  all_raster=lapply(tables_numeric,transform_raster,ex)
  layers=stack(all_raster)
  return(layers)
  
}

run_cpt=function(x,y,run,output,modes_x,modes_y,modes_cca,trans,type_trans){
  
  file_y=read.table(y,sep="\t",dec=".",skip =3,fill=TRUE,na.strings =-999,stringsAsFactors=FALSE)
  p=dim(file_y)[2]-1
  mode_y=modes_y
  if(p<10)mode_y=p
  mode_cca=modes_cca
  if(p<5)mode_cca=p
  
  t=ifelse(trans==1,541," ")
  
  GI=paste0(output,"_GI.txt"); pear=paste0(output,"_pearson.txt"); afc=paste0(output,"_2afc.txt")
  prob=paste0(output,"_prob.txt");roc_a=paste0(output,"_roc_a.txt");roc_b=paste0(output ,"_roc_b.txt")
  pca_eigen_x=paste0(output,"_pca_eigen_x.txt"); pca_load_x=paste0(output,"_pca_load_x.txt"); pca_scores_x=paste0(output,"_pca_scores_x.txt")
  pca_eigen_y=paste0(output,"_pca_eigen_y.txt"); pca_load_y=paste0(output,"_pca_load_y.txt"); pca_scores_y=paste0(output,"_pca_scores_y.txt")
  cca_load_x=paste0(output,"_cca_load_x.txt"); cca_cc=paste0(output,"_cca_cc.txt"); cca_scores_x=paste0(output,"_cca_scores_x.txt")
  cca_load_y=paste0(output,"_cca_load_y.txt"); cca_scores_y=paste0(output,"_cca_scores_y.txt")
  
  hit_s=paste0(output,"_hit_s.txt")
  hit_ss=paste0(output,"_hit_ss.txt")
  
  cmd <- "@echo off
  (
  echo 611
  echo 545
  echo 1
  echo %path_x% 
  echo /
  echo /
  echo /
  echo /
  echo 1
  echo %modex%
  echo 2
  echo %path_y%
  echo /
  echo /
  echo /
  echo /
  echo 1
  echo %modey%
  echo 1
  echo %modecca%
  echo 9
  echo 1
  echo 532
  echo /
  echo /
  echo N
  echo 2
  echo 554
  echo %typetrans%
  echo %trans%
  echo 112
  echo %path_GI%
  echo 311
  echo 451
  echo 455
  echo 413
  echo 1
  echo %path_pear%
  echo 413
  echo 3
  echo %path_2afc%
  echo 413
  echo 4 
  echo %path_hit_s%
  echo 413  
  echo 5
  echo %path_hit_ss% 
  echo 413
  echo 10
  echo %path_roc_b%
  echo 413
  echo 11
  echo %path_roc_a%
  echo 111
  echo 301
  echo %path_pca_eigen_x%
  echo 302
  echo %path_pca_load_x%
  echo 303
  echo %path_pca_scores_x%
  echo 311
  echo %path_pca_eigen_y%
  echo 312
  echo %path_pca_load_y%
  echo 313
  echo %path_pca_scores_y%
  echo 401
  echo %path_cca_cc%
  echo 411
  echo %path_cca_load_x%
  echo 412
  echo %path_cca_scores_x%
  echo 421
  echo %path_cca_load_y%
  echo 422
  echo %path_cca_scores_y%
  echo 501
  echo %path_prob%
  echo 0
  echo 0
  ) | CPT_batch.exe"
  
  cmd<-gsub("%path_x%",x,cmd)
  cmd<-gsub("%path_y%",y,cmd)
  cmd<-gsub("%path_GI%",GI,cmd)
  cmd<-gsub("%path_pear%",pear,cmd)
  cmd<-gsub("%path_2afc%",afc,cmd)
  cmd<-gsub("%path_roc_b%",roc_b,cmd)
  cmd<-gsub("%path_roc_a%",roc_a,cmd)
  cmd<-gsub("%path_prob%",prob,cmd)
  cmd<-gsub("%modey%",mode_y,cmd)
  cmd<-gsub("%modex%",modes_x,cmd)
  cmd<-gsub("%modecca%",mode_cca,cmd)
  cmd<-gsub("%typetrans%",type_trans,cmd)
  cmd<-gsub("%trans%",t,cmd)
  cmd<-gsub("%path_cca_load_x%",cca_load_x,cmd)
  cmd<-gsub("%path_cca_cc%",cca_cc,cmd)
  
  cmd<-gsub("%path_pca_eigen_x%",pca_eigen_x,cmd)
  cmd<-gsub("%path_pca_load_x%",pca_load_x,cmd)
  cmd<-gsub("%path_pca_scores_x%",pca_scores_x,cmd)
  cmd<-gsub("%path_pca_eigen_y%",pca_eigen_y,cmd)
  cmd<-gsub("%path_pca_load_y%",pca_load_y,cmd)
  cmd<-gsub("%path_pca_scores_y%",pca_scores_y,cmd)
  cmd<-gsub("%path_cca_scores_x%",cca_scores_x,cmd)
  cmd<-gsub("%path_cca_scores_y%",cca_scores_y,cmd)
  cmd<-gsub("%path_cca_load_y%",cca_load_y,cmd)
  
  cmd<-gsub("%path_hit_s%",hit_s,cmd)
  cmd<-gsub("%path_hit_ss%",hit_ss,cmd)
  
  
  write(cmd,run)
  #shell.exec(run)
  system2(run)
  
}

correl <- function(x,y){
  
  y[1,1]=""
  loadings <- na.omit(y)
  loadings[loadings==-999]=NA
  pos=which(loadings[,1]=="")
  if(length(pos)==1){list_dates=list(loadings)}else{vector_split <- sort(rep(pos,pos[2]-1));list_dates <- split(loadings,vector_split)}
  tables=lapply(list_dates,"[",-1,-1)
  cor_ca=x[1:length(tables),1]
  final=Reduce("+",Map(function(x,y) abs(x)*y ,tables,cor_ca))/sum(cor_ca)
  final_vec=as.vector(as.matrix(t(final)))
  
  return(final_vec)
}

files_x=function(raster,cor,na,years){
  
  coor_min=apply(coordinates(raster),2,min) 
  coor_max=apply(coordinates(raster),2,max) 
  coor_all=cbind(coor_min,coor_max)
  
  year_p=paste0("cpt:T=",years)
  
  for(i in seq(0.1,0.9,0.1)){
    
    #pos_data=which(!is.na(values(raster)[,1]))
    pos_selec=which(cor<quantile(cor,i,na.rm=T))
    #pos_final=pos_data*pos_selec
    val=values(raster)
    val[pos_selec,]=NA
    val[which(is.na(val),arr.ind = T)]= -999
    val_l=split(val,col(val))
    
    
    lat=sort(seq(coor_all[2,1],coor_all[2,2]),decreasing = T)
    lon=sort(seq(coor_all[1,1],coor_all[1,2]))
    val_matrix=lapply(val_l,function(x)matrix(x,length(lat),length(lon),byrow=TRUE,dimnames=list(lat,lon)))
    
    
    p="xmlns:cpt=http://iri.columbia.edu/CPT/v10/"
    p1="cpt:nfields=1"
    p2=paste0("cpt:field=ssta, ",year_p[1],", cpt:nrow=",length(lat),", cpt:ncol=",length(lon),", cpt:row=Y, cpt:col=X, cpt:units=Kelvin_scale, cpt:missing=-999")
    
    name_file=paste0(na,"_",i,".txt")
    sink(name_file)
    cat(p)
    cat("\n")
    cat(p1) 
    cat("\n")
    cat(p2) 
    cat("\n")
    u=Map(function(x,y){write.table(t(c(" ",lon)),sep="\t",col.names=F,row.names=F,quote = F);write.table(x,sep="\t",col.names=F,row.names=T,quote = F);cat(y);cat("\n")},val_matrix,c(year_p[-1],""))
    sink()
  }
  
  return("Successful process")   
}

best_GI=function(x){
  
  all_path=paste0(x,"_",seq(0,0.9,0.1),"_GI.txt")
  names=substr(basename(all_path),1,nchar(basename(all_path))-4)
  all_GI=lapply(all_path,function(x)read.table(x,header=T,dec=".",skip=5))  
  best=lapply(all_GI,function(x) x[dim(x)[1],dim(x)[2]] )
  pos=which(unlist(best)==max(unlist(best)))[1]
  output=substr(names[pos],1,nchar(names[pos])-3)
  return(output)
  
}

metricas=function(x){
  
  p=read.table(paste0(x,"_pearson.txt"),header=T,dec=".",skip=2,col.names = c("id","latitud","longitud","pearson"))
  k=read.table(paste0(x,"_2afc.txt"),header=T,dec=".",skip=2,col.names = c("id","latitud","longitud","kendall"))
  g=read.table(paste0(x,"_GI.txt"),header=T,dec=".",skip=5)
  goodness=g[dim(g)[1],dim(g)[2]]
  roc_b=read.table(paste0(x,"_roc_b.txt"),header=T,dec=".",skip=2,col.names = c("id","latitud","longitud","roc_b"))
  roc_a=read.table(paste0(x,"_roc_a.txt"),header=T,dec=".",skip=2,col.names = c("id","latitud","longitud","roc_a"))
  
  hit_s=read.table(paste0(x,"_hit_s.txt"),header=T,dec=".",skip=2,col.names = c("id","latitud","longitud","hit_s"))
  hit_ss=read.table(paste0(x,"_hit_ss.txt"),header=T,dec=".",skip=2,col.names = c("id","latitud","longitud","hit_ss"))
  
  hit=merge(hit_s,hit_ss)
  
  roc=merge(roc_b,roc_a)
  p_k=merge(p,k)
  all=merge(p_k,roc)
  all_final=merge(all,hit)
  metrics=cbind(file=basename(x),all_final,goodness)
  
  below=read.table(paste0(x,"_prob.txt"),header=T,nrow=3,dec=".",fill=T,skip=3,check.names = FALSE)[-1:-2,]
  normal=read.table(paste0(x,"_prob.txt"),header=T,nrow=3,dec=".",fill=T,skip=8,check.names = FALSE)[-1:-2,]
  above=read.table(paste0(x,"_prob.txt"),header=T,nrow=3,dec=".",fill=T,skip=13,check.names = FALSE)[-1:-2,]
  coor=read.table(paste0(x,"_prob.txt"),header=T,nrow=2,dec=".",fill=T,skip=3,check.names = FALSE)
  prob=cbind(id=names(below),below=as.matrix(below)[1,],normal=as.matrix(normal)[1,],above=as.matrix(above)[1,])
  
  all_data=merge(metrics,prob)
  
  return(all_data)
}

save_areas=function(ras,cor,all_name){
  
  name=basename(all_name) 
  dec=substr(name,nchar(name)-2,nchar(name))
  cor_raster=ras[[1]]
  values(cor_raster)=cor
  q1=quantile(cor,as.numeric(dec),na.rm=T)
  jBrewColors <- brewer.pal(n = 9, name = "Reds")
  tiff(paste0(all_name,".tiff"),compression = 'lzw',height = 6.5,width = 5.7,units="in", res=150)
  par(mfrow=c(2,1))
  plot(cor_raster,main="Weighted loadings",col=jBrewColors,colNA="gray",legend.width=1,legend.shrink=1)
  plot(cor_raster>= q1,main="Selected pixels",colNA="gray",legend=F,col=jBrewColors)
  dev.off()
  
  return(print("Área seleccionada guardada en formato Raster"))
  
}

eigen_plot <- function(path_raw, path_output){
  xeigen<- read.csv(paste0(path_raw, "_pca_eigen_x.txt"),skip =2, header=T, sep="")
  yeigen <- read.csv(paste0(path_raw,"_pca_eigen_y.txt"),skip =2, header=T, sep="")
  
  
  datos_e <- data.frame(modes=xeigen$Mode, eigenx=xeigen$variance, eigeny=yeigen$variance, row.names = NULL)
  
  # grÃ¡fico de los modos 
  modosx <-  ggplot(datos_e, aes(modes,group = 1)) +   geom_line(aes(y = eigenx ),  colour="firebrick3" ) + geom_point(aes(y = eigenx ),  colour="firebrick3" ) +
    theme_bw() + theme( title =element_text(size=12, face='bold'),axis.text.y = element_text(size=12),  legend.position = "none", axis.text.x = element_text(angle = 0, hjust = 1, size = 12)) +
    guides(colour = guide_legend(title = " ")) + labs(x="Mode",y="% variance",title = "X Scree Plot") 
  
  
  modosy <-  ggplot(datos_e, aes(modes,group = 1)) +   geom_line(aes(y = eigeny ),  colour="firebrick3" ) + geom_point(aes(y = eigeny ),  colour="firebrick3" ) +
    theme_bw() + theme( title =element_text(size=12, face='bold'),axis.text.y = element_text(size=12),  legend.position = "none", axis.text.x = element_text(angle = 0, hjust = 1, size = 12)) +
    guides(colour = guide_legend(title = " ")) + labs(x="Mode",y="% variance",title = "Y Scree Plot") 
  
  layt<-grid.layout(nrow=1,ncol=2)
  trim_n = unlist(strsplit(path_raw,"/")) 
  trim_n = trim_n[length(trim_n)]
  
  tiff(filename = paste0(path_output,"/",trim_n, "_eigen_plot.tif"), width = 1500, height = 800,res=150,compression = 'lzw')
  grid.newpage()
  pushViewport(viewport(layout=layt))
  print(modosx,vp=viewport(layout.pos.row=1,layout.pos.col=1))
  print(modosy,vp=viewport(layout.pos.row=1,layout.pos.col=2))
  dev.off()
  cat(paste0(trim_n)," Scree plots realizados...\n")
  
}

cca_map <- function(path_raw , path_output,i, coor) {
  
  xserie <- read.csv(paste0(path_raw, "_cca_scores_x.txt"),skip =2, header=T, sep="")
  yserie <- read.csv(paste0(path_raw,"_cca_scores_y.txt"),skip =2, header=T, sep="")
  
  yloadcca <-  read.csv(paste0(path_raw, "_cca_load_y.txt"),skip =2, header=T, sep="")
  
  ## =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= ##
  ## CCA Maps
  ## =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= ##
  x = tables[[i]]
  
  ext = c(min(as.numeric(coor[[2]])),max(as.numeric(coor[[2]])),min(as.numeric(coor[[1]])),max(as.numeric(coor[[1]])))
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
    labs(title=paste0("X Spatial Loadings (Mode ",i,")"),x="",y=" ", fill = " ")  + theme(legend.key.height=unit(0.5,"cm"),legend.key.width=unit(2,"cm"),
                                                                                          legend.text=element_text(size=10),
                                                                                          panel.background=element_rect(fill="white",colour="black"),
                                                                                          axis.text=element_text(colour="black",size=12),
                                                                                          axis.title=element_text(colour="black",size=12,face="bold"),
                                                                                          legend.position = "bottom", 
                                                                                          legend.title = element_text(size = 12.5))  +
    scale_fill_gradientn(colours =myPalette(100), limits=c(-1,1))
  
  
  yloadcca$val = yloadcca[,i+3]/max(abs(yloadcca[,i+3]),na.rm=T)
  y = yloadcca
  ext_y = c(min(as.numeric(y$Latitude))-0.5,max(as.numeric(y$Latitude))+0.5,min(as.numeric(y$Longitude))-0.5,max(as.numeric(y$Longitude))+0.5)
  
  sPDF <<- getMap()  
  new_map = fortify(sPDF)
  sel <<-new_map[new_map$lat>ext_y[1] & new_map$lat<ext_y[2] & new_map$long>ext_y[3] & new_map$long<ext_y[4],]
  
  
  p <- ggplot(sel, aes(x=long,y=lat)) 
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
                  legend.title = element_text(size = 12)) + labs(title=paste0("Y Spatial Loadings (Mode ",i,")"),  x= "", y= "", colour="") 
  
  
  ## GrÃ¡ficos Componentes 
  
  # Se crea una trama de datos con la fecha y las componentes 
  datos <- data.frame(X=xserie[,i], Y=yserie[,i], row.names = NULL)
  datos$X = round(datos$X ,4) # redondee los modos 
  datos$Y = round(datos$Y ,4) # redondee los modos 
  datos$date = as.numeric(substring(rownames(xserie),1,4))
  
  datos[datos$X==-999.0000,1:2]=0
  datos[datos$Y==-999.0000,1:2]=0 # quite los valore NA
  # quite los valore NA
  datos[,1:2]=datos[,1:2]*100 # multipliquelos * 100
  
  modos <-  ggplot(datos, aes(date,group = 1)) +   geom_line(aes(y = X ),  colour="firebrick3" ) + 
    geom_line(aes(y = Y),  colour="forestgreen")  + 
    geom_hline(yintercept = 0, colour="gray") + theme_bw() + 
    theme( title =element_text(size=12, face='bold'),axis.text.y = element_text(size=12),  legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1, size = 12)) +
    guides(colour = guide_legend(title = " ")) + labs(subtitle=paste("Canonical Correlation = ", round(cor(datos$X,datos$Y),3),sep = ""),x="",y="Scores (X red; Y green) (*100)",
                                                      title = paste0("Temporal Scores (Mode " ,i,")")) + scale_x_continuous(breaks = seq(min(datos$date),max(datos$date),3))
  
  
  layt<-grid.layout(nrow=1,ncol=3,widths=c(4/9,2.5/9, 2.5/9),default.units=c('null','null'))
  
  trim_n = unlist(strsplit(path_raw,"/")) 
  trim_n = trim_n[length(trim_n)]
  
  tiff(filename = paste0(path_output, "/",trim_n,"_mode_",i,"_cca_maps.tif"), width = 1700, height = 400,res=100,compression = 'lzw')
  grid.newpage()
  pushViewport(viewport(layout=layt))
  print(Map_x,vp=viewport(layout.pos.row=1,layout.pos.col=1))
  print(modos,vp=viewport(layout.pos.row=1,layout.pos.col=2))
  print(p,vp=viewport(layout.pos.row=1,layout.pos.col=3))
  dev.off()
  cat(paste0(trim_n), " Mapas CCA realizados...\n")
  
  
  
}

cca_map_all <- function(path_raw,path_output){
  
  y <- read.table(paste0(path_raw, "_cca_load_x.txt"),sep="\t",dec=".",skip =2,fill=TRUE,na.strings =-999,stringsAsFactors=FALSE)
  y[1,1]=""
  loadings <- na.omit(y)
  loadings[loadings==-999]=NA
  pos=which(loadings[,1]=="")
  if(length(pos)==1){list_dates=list(loadings)}else{vector_split <- sort(rep(pos,pos[2]-1));list_dates <- split(loadings,vector_split)}
  coor <- list(list_dates[[1]][1,][-1],list_dates[[1]][,1][-1])
  
  tables <<- lapply(list_dates,"[",-1,-1)
  
  for(i in 1:length(tables)) cca_map(path_raw,path_output,i,coor)
}

metric_map <- function(path_metric, path_output,path_raw){
  
  yloadcca <-  read.csv(paste0(path_raw, "_cca_load_y.txt"),skip =2, header=T, sep="")
  
  y = yloadcca
  ext_y = c(min(as.numeric(y$Latitude))-0.5,max(as.numeric(y$Latitude))+0.5,min(as.numeric(y$Longitude))-0.5,max(as.numeric(y$Longitude))+0.5)
  
  sPDF <<- getMap()  
  new_map = fortify(sPDF)
  sel <<-new_map[new_map$lat>ext_y[1] & new_map$lat<ext_y[2] & new_map$long>ext_y[3] & new_map$long<ext_y[4],]
  
  
  all_ind_complete <-read.csv(paste0(path_metric, '/metrics.csv')) 
  
  trimesters <- unique(all_ind_complete$file)
  
  for(i in 1:length(trimesters)){
    
    
    
    all_ind <- all_ind_complete[all_ind_complete$file==trimesters[i],]
    
    myPalette <-  colorRampPalette(c("dodgerblue4", "dodgerblue1","deepskyblue","darkslategray1", "lightcyan",  "lemonchiffon1","khaki","sandybrown", "darkorange2","firebrick2"))
    
    ind <- ggplot(sel, aes(x=long,y=lat)) + 
      geom_polygon(aes(group=group),fill="snow") + 
      geom_point(data=all_ind, aes(x= longitud, y= latitud, group= 1 ,col=kendall),size=3) + 
      geom_path(aes(long,lat,group=group),color="black",size=0.3)  + 
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
    cat(paste0(trimesters[i])," Mapas Metricas realizados...\n")
    
    
    layt<-grid.layout(nrow=2,ncol=1)
    
    tiff(filename = paste0(path_output,"/" ,trimesters[i], "_roc_maps.tif"), width = 800, height = 800,res=130,compression = 'lzw')
    grid.newpage()
    pushViewport(viewport(layout=layt))
    print(ind_2,vp=viewport(layout.pos.row=1,layout.pos.col=1))
    print(ind_3,vp=viewport(layout.pos.row=2,layout.pos.col=1))
    
    dev.off()
    cat(paste0(trimesters[i])," Mapas ROC realizados...\n")
    
    layt<-grid.layout(nrow=2,ncol=1)
    
    tiff(filename = paste0(path_output,"/" ,trimesters[i], "_hit_maps.tif"), width = 800, height = 800,res=130,compression = 'lzw')
    grid.newpage()
    pushViewport(viewport(layout=layt))
    print(ind_4,vp=viewport(layout.pos.row=1,layout.pos.col=1))
    print(ind_5,vp=viewport(layout.pos.row=2,layout.pos.col=1))
    
    dev.off()
    cat(paste0(trimesters[i])," Mapas Hit realizados...\n")
    
    max_C<-apply(all_ind[,c("below","normal","above")], 1, max)
    cat<-ifelse(all_ind[,"below"]==max_C, "Below", ifelse(all_ind[,"normal"]==max_C, "Normal", ifelse(all_ind[,"above"]==max_C, "Above",0)))
    maximos<-cbind.data.frame(all_ind$id, all_ind$longitud, all_ind$latitud, max_C, cat)
    
    # Aqui se ingresan los datos de las estaciones
    maximos$cat = factor(maximos$cat, levels = c("Below", "Normal", "Above"))
    p <- ggplot(sel, aes(x=long,y=lat)) + 
      geom_polygon(aes(fill=hole,group=group),fill="snow")
    
    maxi <- p + geom_point(data=maximos, aes(x=all_ind$longitud, y=all_ind$latitud, colour=cat))+
      geom_path(aes(long,lat,group=group),color="black",size=0.3) + scale_colour_manual(values = c("tomato2","lightgreen","steelblue3")) +
      coord_equal() + theme( legend.key.height=unit(1,"cm"),legend.key.width=unit(0.5,"cm"),
                             legend.text=element_text(size=8),
                             panel.background=element_rect(fill="white",colour="black"),
                             axis.text=element_text(colour="black",size=10),
                             axis.title=element_text(colour="black",size=10,face="bold"),
                             #legend.position = "bottom", 
                             legend.title=element_blank())  + labs(title="Probabilistic Forecast")+labs( x=" ", y=" ", size=" ")
    
    tiff(paste0(path_output,"/" ,trimesters[i], "_probabilistic_maps.tif"),width = 3000, height = 3000, units = "px", res = 400,compression = 'lzw')
    print(maxi)
    dev.off()  
    cat(paste0(trimesters[i]), " Mapas Probabilistico realizados...\n")
    
    
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
    cat(paste0(trimesters[i])," Mapas Probabilidades realizados...\n")
  }
  
  
}

######## run ##############

folders <- list.files(main_dir,full.names = T)
lapply(folders,function(x) dir.create(paste0(x,"/output/raw_output"),recursive = T))
lapply(folders,function(x) dir.create(paste0(x,"/output/all_domain"),recursive = T))
lapply(folders,function(x) dir.create(paste0(x,"/output/opt_domain"),recursive = T))
lapply(folders,function(x) dir.create(paste0(x,"/bat_files")))
path_x <- lapply(folders,function(x)list.files(paste0(x,"/input/sst_cfsv2"),full.names = T))
names_x <- lapply(path_x,function(x) substr(basename(x),1,nchar(basename(x))-4))
path_y <- lapply(folders,function(x)list.files(paste0(x,"/input/stations"),full.names = T))
path_output <- Map(function(x,y) paste0(x,"/output/raw_output/",y,"_0"),folders,names_x)
path_run <- Map(function(x,y) paste0(x,"/bat_files/",y,"_0",".bat"),folders,names_x)

cat("\n Directorios cargados y carpetas creadas \n")

first_run <- Map(function(x,y,z,k,p1,p2,p3,p4,p5)Map(run_cpt,x,y,z,k,p1,p2,p3,p4,p5),path_x,path_y,path_run,path_output,modes_x,modes_y,modes_cca,trans,type_trans)

cat("\n Primera corrida realizada")

tsm_list <- lapply(path_x,function(x)lapply(x,function(x1)read.table(x1,sep="\t",dec=".",skip =2,fill=TRUE,na.strings =-999,stringsAsFactors=FALSE)))
time=lapply(tsm_list,function(x)lapply(x,function(x1) as.character(x1[1,])[-1]))
time_sel=lapply(time,function(x)lapply(x,function(x1)x1[x1!="NA"]))
tsm_raster <- lapply(tsm_list,function(x)lapply(x,data_raster))

cat("\n Datos cargados en formato raster")

path_cc <- lapply(paste0(folders,"/output/raw_output"),function(x)list.files(x,full.names = T,pattern = "cca_cc"))
path_load <- lapply(paste0(folders,"/output/raw_output"),function(x)list.files(x,full.names = T,pattern = "cca_load_x"))
cc <-  lapply(path_cc,function(x)lapply(x,function(x1)read.table(x1,sep="\t",dec=".",header = T,row.names = 1,skip =2,fill=TRUE,na.strings =-999,stringsAsFactors=FALSE)))
load <- lapply(path_load,function(x)lapply(x,function(x1)read.table(x1,sep="\t",dec=".",skip =2,fill=TRUE,na.strings =-999,stringsAsFactors=FALSE)))
cor_tsm <- Map(function(x,y)Map(correl,x,y),cc,load)

cat("\n Correlación calculada")

names_selec <-Map(function(x,y) paste0(x,"/input/sst_cfsv2/",substr(y,1,nchar(y))) ,folders,names_x)
o_empty_1=Map(function(x,y,z,r)Map(files_x,x,y,z,r),tsm_raster,cor_tsm,names_selec,time_sel)

cat("\n Archivos de la TSM construidos por deciles para CPT \n")

path_x_2 <- lapply(folders,function(x)list.files(paste0(x,"/input/sst_cfsv2"),full.names = T,pattern = "0."))
names_x_2 <- lapply(path_x_2,function(x) substr(basename(x),1,nchar(basename(x))-4))
path_run_2 <- Map(function(x,y) paste0(x,"/bat_files/",y,".bat"),folders,names_x_2)
path_output_2 <- Map(function(x,y) paste0(x,"/output/raw_output/",y),folders,names_x_2)
o_empty_2 <-Map(function(x,y,z,k,p1,p2,p3,p4,p5)Map(run_cpt,x,y,z,k,p1,p2,p3,p4,p5),path_x_2,path_y,path_run_2,path_output_2,modes_x,modes_y,modes_cca,trans,type_trans)

cat("\n Segunda corrida realizada\n")

folder_output <- Map(function(x,y) paste0(x,"/output/raw_output/",substr(y,1,nchar(y))),folders,names_x)
best_decil=lapply(folder_output,function(x)lapply(x,best_GI)) %>% lapply(.,unlist)

cat("\n Mejor corrida seleccionada \n")

best_path <- Map(function(x,y) paste0(x,"/output/raw_output/",y),folders,best_decil)
all_metricas<- lapply(best_path,function(x)lapply(x,metricas))
metricas_final=lapply(all_metricas,function(x)do.call("rbind",x))
o_empty_3=Map(function(x,y) write.csv(x,paste0(y,"/output/opt_domain/metrics.csv"),row.names=FALSE),metricas_final,folders)

cat("\n Metricas con dominio optimizado almacenadas \n")

normal_path <- Map(function(x) paste0(x,"_0"),folder_output)
all_metricas_n<- lapply(normal_path,function(x)lapply(x,metricas))
metricas_final_n=lapply(all_metricas_n,function(x)do.call("rbind",x))
o_empty_4=Map(function(x,y) write.csv(x,paste0(y,"/output/all_domain/metrics.csv"),row.names=FALSE),metricas_final_n,folders)

cat("\n Metricas con todo el dominio almacenadas \n")

path_images <- Map(function(x,y) paste0(x,"/output/opt_domain/",y),folders,best_decil)
o_empty_5 <- Map(function(x,y,z)Map(save_areas,x,y,z),tsm_raster,cor_tsm,path_images)

cat("\n Pixeles selecionados almacenados en .tiff  \n")

# Run all_domain
path_metric <-  paste0(folders,"/output/all_domain")
path_out <-  paste0(folders,"/output/all_domain")
path_raw <- normal_path

path_raw_m <- lapply(path_raw, "[[",1)
Map(function(x,y,z)Map(metric_map,x,y,z),path_metric,path_out,path_raw_m)

Map(function(x,y)Map(eigen_plot,x,y),path_raw,path_out)
Map(function(x,y)Map(cca_map_all,x,y),path_raw,path_out)

#Run opt_domain
path_metric <-  paste0(folders,"/output/opt_domain")
path_out <-  paste0(folders,"/output/opt_domain")
path_raw <- best_path

path_raw_m <- lapply(path_raw, "[[",1)
Map(function(x,y,z)Map(metric_map,x,y,z),path_metric,path_out,path_raw_m)

Map(function(x,y)Map(eigen_plot,x,y),path_raw,path_out)
Map(function(x,y)Map(cca_map_all,x,y),path_raw,path_out)

cat("\n Gráficos almacenados en .tiff  \n")





