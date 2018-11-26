# Script to automate runs from CPT(CFSV2) 
#  making predictor area selection 

# Created by: Diego Fernando Agudelo (d.agudelo@cgiar.org)
# Date: July 2018

########### Packages ###############

if(require(stringr)==FALSE){install.packages("stringr",dependencies = TRUE)}
library("stringr")
if(require(corpcor)==FALSE){install.packages("corpcor")}
library("corpcor")
if(require(pcaPP)==FALSE){install.packages("pcaPP")}
library("pcaPP")
if(require(raster)==FALSE){install.packages("raster")}
library("raster")
if(require(RColorBrewer)==FALSE){install.packages("RColorBrewer")}
library("RColorBrewer")
if(require(dplyr)==FALSE){install.packages("dplyr")}
library("dplyr")

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

run_cpt=function(x,y,run,output){
  
  file_y=read.table(y,sep="\t",dec=".",skip =3,fill=TRUE,na.strings =-999,stringsAsFactors=FALSE)
  p=dim(file_y)[2]-1
  mode_y=10
  if(p<10)mode_y=p
  mode_cca=5
  if(p<5)mode_cca=p
  
  GI=paste0(output,"_GI.txt"); pear=paste0(output,"_pearson.txt"); afc=paste0(output,"_2afc.txt")
  prob=paste0(output,"_prob.txt");roc_a=paste0(output,"_roc_a.txt");roc_b=paste0(output ,"_roc_b.txt")
  pca_eigen_x=paste0(output,"_pca_eigen_x.txt"); pca_load_x=paste0(output,"_pca_load_x.txt"); pca_scores_x=paste0(output,"_pca_scores_x.txt")
  pca_eigen_y=paste0(output,"_pca_eigen_y.txt"); pca_load_y=paste0(output,"_pca_load_y.txt"); pca_scores_y=paste0(output,"_pca_scores_y.txt")
  cca_load_x=paste0(output,"_cca_load_x.txt"); cca_cc=paste0(output,"_cca_cc.txt"); cca_scores_x=paste0(output,"_cca_scores_x.txt")
  cca_load_y=paste0(output,"_cca_load_y.txt"); cca_scores_y=paste0(output,"_cca_scores_y.txt")
  
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
  echo 10
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
  echo 2
  echo 541
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
  cmd<-gsub("%modecca%",mode_cca,cmd)
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
  
  roc=merge(roc_b,roc_a)
  p_k=merge(p,k)
  all=merge(p_k,roc)
  metrics=cbind(file=basename(x),all,goodness)
  
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

######## run ##############

main_dir <- "D:/OneDrive - CGIAR/Desktop/Codigos_TNC/cpt_r"
setwd(main_dir)
folders <- list.files(".",full.names = T)
lapply(folders,function(x) dir.create(paste0(x,"/output/raw_output"),recursive = T))
lapply(folders,function(x) dir.create(paste0(x,"/output/all_domain"),recursive = T))
lapply(folders,function(x) dir.create(paste0(x,"/output/opt_domain"),recursive = T))
lapply(folders,function(x) dir.create(paste0(x,"/bat_files")))
path_x <- lapply(folders,function(x)list.files(paste0(x,"/input/cfsv2"),full.names = T))
names_x <- lapply(path_x,function(x) substr(basename(x),1,nchar(basename(x))-4))
path_y <- lapply(folders,function(x)list.files(paste0(x,"/input/stations"),full.names = T))
path_output <- Map(function(x,y) paste0(x,"/output/raw_output/",y,"_0"),folders,names_x)
path_run <- Map(function(x,y) paste0(x,"/bat_files/",y,"_0",".bat"),folders,names_x)

cat("\n Directorios cargados y carpetas creadas \n")

first_run <- Map(function(x,y,z,k)Map(run_cpt,x,y,z,k),path_x,path_y,path_run,path_output)

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

names_selec <-Map(function(x,y) paste0(x,"/input/cfsv2/",substr(y,1,nchar(y))) ,folders,names_x)
o_empty_1=Map(function(x,y,z,r)Map(files_x,x,y,z,r),tsm_raster,cor_tsm,names_selec,time_sel)

cat("\n Archivos de la TSM construidos por deciles para CPT \n")

path_x_2 <- lapply(folders,function(x)list.files(paste0(x,"/input/cfsv2"),full.names = T,pattern = "0."))
names_x_2 <- lapply(path_x_2,function(x) substr(basename(x),1,nchar(basename(x))-4))
path_run_2 <- Map(function(x,y) paste0(x,"/bat_files/",y,".bat"),folders,names_x_2)
path_output_2 <- Map(function(x,y) paste0(x,"/output/raw_output/",y),folders,names_x_2)
o_empty_2 <-Map(function(x,y,z,k)Map(run_cpt,x,y,z,k),path_x_2,path_y,path_run_2,path_output_2)

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

