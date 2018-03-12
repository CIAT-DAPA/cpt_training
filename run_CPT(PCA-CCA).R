# Script to generate a run of CPT since R
# Created by: Diego Fernando Agudelo (d.agudelo@cgiar.org)
# Date: March 2018

########### Packages ###############

library(raster)
library(rgdal)
library(corpcor)

########### Functions ##############

transform_raster=function(x,y){
  mapa_base=raster(ext=y, res=c(1,1))
  val=c(as.matrix(t(x),ncol=1,byrow = T))
  val=as.numeric(val)
  val[val==-999]=NA
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

quarterly_data=function(names,data){
  data=data[-1:-2,]
  trim=substr(names,5,nchar(names)-4)
  months=strsplit(trim,"-")
  number_months=sapply(months[[1]],function(x)which(month.abb==x))
  
  year=as.numeric(unlist(lapply(strsplit(data[-1,1],"-"),"[", 1),recursive=FALSE))
  month=as.numeric(unlist(lapply(strsplit(data[-1,1],"-"),"[", 2),recursive=FALSE))
  new_data=cbind.data.frame(year=year,month=month,sapply(data[-1,-1],as.numeric))
  
  pos_ini=which(new_data$month==number_months[1])
  if(length(number_months)==2){
    pos_data=sort(c(pos_ini,pos_ini+1))
  }else{pos_data=sort(c(pos_ini,pos_ini+1,pos_ini+2))}
  data_out=data.frame(na.omit(aggregate(new_data[pos_data,-1:-2],by=list(sort(rep(1:(length(pos_ini)),length(number_months)))),sum))[,-1])
  year_out=new_data[pos_ini+1,1][1:dim(data_out)[1]]
  rownames(data_out)=year_out
  
  return(data_out)
} 

pca_x_svd=function(x,modos,pon){
  
  n=dim(x)[1]
  p=dim(x)[2]
  #pon_matrix=diag(pon)
  pon_matrix=matrix(rep(pon,n),n,p,byrow = TRUE)
  x0=scale(x)#*(sqrt(n)/sqrt(n-1))
  #x_pon=t(pon_matrix%*%t(x0))
  x_pon=pon_matrix*x0
  svd_o <- fast.svd(x_pon)
  comp  <- svd_o$u[,1:modos,drop=F]%*%diag(svd_o$d[1:modos],length(svd_o$d[1:modos]),length(svd_o$d[1:modos])) 
  vect  <- svd_o$v[,1:modos]
  output <- list(comp,vect)
  return(output)
  
}

pca_y_svd=function(x,modos){
  
  n=dim(x)[1]
  x0=scale(x)#*(sqrt(n)/sqrt(n-1))
  svd_o <- fast.svd(x0)
  comp  <- svd_o$u[,1:modos,drop=F]%*%diag(svd_o$d[1:modos],length(svd_o$d[1:modos]),length(svd_o$d[1:modos])) 
  vect  <- svd_o$v[,1:modos]
  output <- list(comp,vect)
  return(output)
  
}

########## run ############

main_dir <- "C:/Users/dagudelo/Desktop/ppts_honduras/CPT/" #### Cambiar Ruta 
path_tsm <- list.files(paste0(main_dir,"/input/CFSV2"),full.names = T)
names_x <- lapply(path_tsm,basename)
data_tsm <- read.table(path_tsm,sep="\t",dec=".",skip =2,fill=TRUE,na.strings =-999,stringsAsFactors=FALSE)
data_tsm_raster <- data_raster(data_tsm)
data_table=t(getValues(data_tsm_raster))
lat=coordinates(data_tsm_raster)[,2]
ponde <- sqrt(cos((pi/180)*lat))
data_table_na=is.na(data_table)
pos=(colSums(data_table_na)==0)
data_table_final=data_table[,pos]
ponde_final=ponde[pos]

cat("\n Datos de la TSM organizados en formato Años X Pixeles \n")

####### Organizar periodo en data_table_final
pca_x <- pca_x_svd(data_table_final[1:34,],5,ponde_final)[[1]] 
View(pca_x)

cat("\n PCA de los datos x generado \n")

path_y <- list.files(paste0(main_dir,"/input/stations"),full.names = T)
data_y <- read.table(path_y,sep="\t",dec=".",skip =3,fill=TRUE,na.strings =-999,stringsAsFactors=FALSE)
prec_cumu <- quarterly_data(names_x,data_y)
View(prec_cumu)

cat("\n Datos de precipitacion organizados de forma trimestral \n")

x <- data_table_final[1:34,] #### Organizar periodo de entrenamiento
y <- prec_cumu[2:35,] #### Organizar periodo de entrenamiento
pon <- ponde_final
modos_x <- 2 ##### Cambiar modos x
modos_y <- 2 #####  Cambiar modos y
 

n=dim(y)[1]
y_fores_final=matrix(NA,n,dim(y)[2])
y_fores=0

for( i in 1:n){
  
  cross=c(n-1,n,1:n,1,2)
  
  X0=x[-cross[i:(i+4)],]  ;  Y0=y[-cross[i:(i+4)],]  ;  X_for=x[cross[i+2],]
  
  mean_x=colMeans(X0)
  sd_x=apply(X0,2,sd)
  mean_y=colMeans(Y0)
  sd_y=apply(Y0,2,sd)
  
  pca_x=pca_x_svd(X0,modos=2,pon)
  pca_y=pca_y_svd(Y0,modos=2)
  
  x_for_scale=(X_for-mean_x)/sd_x
  x_for_pon=x_for_scale*pon
  x_for_comp=x_for_pon%*%pca_x[[2]]
  
  canonico=cancor(pca_x[[1]],pca_y[[1]])
  x_center=scale(pca_x[[1]],scale=F)
  y_center=scale(pca_y[[1]],scale=F)  
  com_x=x_center%*%canonico$xcoef
  com_y=y_center%*%canonico$ycoef
  R=canonico$cor
  xi=as.matrix(x_for_comp-canonico$xcenter)
  vec_x=canonico$xcoef
  xi_com=as.matrix(xi)%*%vec_x[,1,drop=FALSE]
  y_fores=(R[1]*xi_com)
  y_fores_last=(y_fores%*%solve(canonico$ycoef)[1,])+canonico$ycenter
  
  
  y_fores_comp=y_fores_last%*%t(as.matrix(pca_y[[2]]))
  y_fores_final[i,]=(y_fores_comp*sd_y)+mean_y
  
}

cat("\n Validación cruzada realizada \n")

View(y_fores_final)
mean(diag(cor(y_fores_final,y,method = "kendall")))




