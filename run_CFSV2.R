# Script to automate runs from CPT(CFSv2) 
#  making predictora area selection 

# Created by: Diego Fernando Agudelo (d.agudelo@cgiar.org)
# Date: March 2018

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

quarterly_data=function(names,data){
  data=data[-1:-2,]
  trim=substr(names,5,nchar(names)-4)
  months=strsplit(trim,"-")
  number_months=which(month.abb%in%months[[1]])
  
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

selection_area=function(x,y,ponde){
  
  mode_x=10
  mode_y=10 ; if(dim(y)[2] < mode_y)mode_y=dim(y)[2]
  
  y_pca=pca_y_svd(y,modos=mode_y)[[1]]
  x_pca=pca_x_svd(x,modos=mode_x,ponde)[[1]]
  all_cor=matrix(NA,mode_x*mode_y,dim(x)[2])
  count=0
  
  for(i in 1:mode_x){
    
    for(j in 1:mode_y){
      
      canonico=cancor(x_pca[,1:i],y_pca[,1:j,drop=F])
      x_center=scale(x_pca[,1:i],scale = F)
      y_center=scale(y_pca[,1:j],scale = F)  
      com_x=x_center%*%canonico$xcoef
      com_y=y_center%*%canonico$ycoef
      mode1=cbind(com_x[,1],com_y[,1])
      cor_tsm=cor(x,mode1[,1])
      count=count+1
      all_cor[count,]=cor_tsm[,1]
      
    }
    
  }
  print("Finalizo la selección del area para un mes")
  cor_mean=apply(abs(all_cor),2,mean)
  return(cor_mean)
}

files_x=function(raster,cor,na,years){
  
  coor_min=apply(coordinates(raster),2,min) 
  coor_max=apply(coordinates(raster),2,max) 
  coor_all=cbind(coor_min,coor_max)
  
  year_p=paste0("cpt:T=",years)
  
  for(i in seq(0,0.9,0.1)){
    
    pos_data=which(!is.na(values(raster)[,1]))
    pos_selec=cor<quantile(cor,i)
    pos_final=pos_data*pos_selec
    val=values(raster)
    val[pos_final,]=NA
    val[which(is.na(val),arr.ind = T)]=0
    val_l=split(val,col(val))
    val_matrix=lapply(val_l,function(x)matrix(x,length(coor_all[2,2]:coor_all[2,1]),length(coor_all[1,1]:coor_all[1,2]),byrow=TRUE,dimnames=list(coor_all[2,2]:coor_all[2,1],coor_all[1,1]:coor_all[1,2])))
    
    
    p="xmlns:cpt=http://iri.columbia.edu/CPT/v10/"
    p1="cpt:nfields=1"
    p2=paste0("cpt:field=ssta, ",year_p[1],", cpt:nrow=",length(coor_all[2,2]:coor_all[2,1]),", cpt:ncol=",length(coor_all[1,1]:coor_all[1,2]),", cpt:row=Y, cpt:col=X, cpt:units=Kelvin_scale, cpt:missing=0")
    
    name_file=paste0(na,"_",i,".txt")
    sink(name_file)
    cat(p)
    cat("\n")
    cat(p1) 
    cat("\n")
    cat(p2) 
    cat("\n")
    u=Map(function(x,y){write.table(t(c(" ",coor_all[1,1]:coor_all[1,2])),sep="\t",col.names=F,row.names=F,quote = F);write.table(x,sep="\t",col.names=F,row.names=T,quote = F);cat(y);cat("\n")},val_matrix,c(year_p[-1],""))
    sink()
  }
  
  return("Successful process")   
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
  
  write(cmd,run)
  shell.exec(run)
  #system(run, ignore.stdout = T, show.output.on.console = T)
  
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
  
  p=read.table(paste0(x,"_pearson.txt"),header=T,dec=".",skip=2,col.names = c("id","Latitud","Longitud","pearson"))
  k=read.table(paste0(x,"_2afc.txt"),header=T,dec=".",skip=2,col.names = c("id","Latitud","Longitud","kendall"))
  g=read.table(paste0(x,"_GI.txt"),header=T,dec=".",skip=5)
  goodness=g[dim(g)[1],dim(g)[2]]
  p_k=merge(p,k)
  data=cbind(file=basename(x),p_k,goodness)
  
  return(data)
}

proba=function(x){
  
  below=read.table(paste0(x,"_prob.txt"),header=T,nrow=3,dec=".",fill=T,skip=3,check.names = FALSE)[-1:-2,]
  normal=read.table(paste0(x,"_prob.txt"),header=T,nrow=3,dec=".",fill=T,skip=8,check.names = FALSE)[-1:-2,]
  above=read.table(paste0(x,"_prob.txt"),header=T,nrow=3,dec=".",fill=T,skip=13,check.names = FALSE)[-1:-2,]
  coor=read.table(paste0(x,"_prob.txt"),header=T,nrow=2,dec=".",fill=T,skip=3,check.names = FALSE)
  data=cbind(file=basename(x),id=names(below),latitud=as.matrix(coor)[1,],longitud=as.matrix(coor)[2,],below=as.matrix(below)[1,],normal=as.matrix(normal)[1,],above=as.matrix(above)[1,])
  
  return(data)
}

save_areas=function(ras,cor,all_name){
  
  name=basename(all_name) 
  dec=substr(name,nchar(name)-2,nchar(name))
  cor_raster=ras[[1]]
  pos=!is.na(values(cor_raster))
  values(cor_raster)[pos]=cor
  q1=quantile(cor,as.numeric(dec))
  jBrewColors <- brewer.pal(n = 9, name = "Reds")
  tiff(paste0(all_name,".tiff"),compression = 'lzw',height = 6.5,width = 5.7,units="in", res=150)
  par(mfrow=c(2,1))
  plot(cor_raster,main="Correlación promedio",col=jBrewColors,colNA="gray",legend.width=1,legend.shrink=1)
  plot(cor_raster>= q1,main="Pixeles seleccionados",colNA="gray",legend=F,col=jBrewColors)
  dev.off()
  
  return(print("Área seleccionada guardada en formato Raster"))
  
}

######## run ###########

main_dir <- "C:/Users/dagudelo/Desktop/Codigos_TNC"
folders <- list.files(paste0(main_dir,"/CPT"),full.names = T)
lapply(folders,function(x) dir.create(paste0(x,"/CFSV2_r")))
lapply(folders,function(x) dir.create(paste0(x,"/output")))
lapply(folders,function(x) dir.create(paste0(x,"/run")))
path_x <- lapply(folders,function(x)list.files(paste0(x,"/CFSV2"),full.names = T))
names_x <- lapply(path_x,basename)
path_y <- lapply(folders,function(x)list.files(paste0(x,"/stations"),full.names = T))

cat("\n Directorios cargados y carpetas creadas \n")

prec_list <- lapply(path_y,function(x1)read.table(x1,sep="\t",dec=".",skip =3,fill=TRUE,na.strings =-999,stringsAsFactors=FALSE))

cat("\n Datos de precipitacion cargados \n")

tsm_list <- lapply(path_x,function(x)lapply(x,function(x1)read.table(x1,sep="\t",dec=".",skip =2,fill=TRUE,na.strings =-999,stringsAsFactors=FALSE)))
time=lapply(tsm_list,function(x)lapply(x,function(x1) as.character(x1[1,])[-1]))
time_sel=lapply(time,function(x)lapply(x,function(x1)x1[x1!="NA"]))
tsm_raster <- lapply(tsm_list,function(x)lapply(x,data_raster))
tsm_list_table1 <- lapply(tsm_raster,function(x1)lapply(x1,function(x) t(getValues(x))))
tsm_list_table_na <- lapply(tsm_list_table1,function(x1)lapply(x1,function(x) is.na(x)))
pos <- lapply(tsm_list_table_na,function(x1)lapply(x1,function(x) colSums(x)==0 ))
tsm_list_table <- Map(function(x,y)Map(function(x1,y1)x1[,y1],x,y),tsm_list_table1,pos)
years_predictor <- lapply(tsm_list_table,function(x1)lapply(x1,function(x)as.numeric(substr(rownames(x),2,nchar(rownames(x))))))


cat("\n Datos de la TSM organizados en formato Años X Pixeles \n")

prec_cumu <- Map(function(x,y)lapply(x,function(x)quarterly_data(x,y)),names_x,prec_list)
years_response <- lapply(prec_cumu,function(x1)lapply(x1,function(x)as.numeric(row.names(x))))

cat("\n Datos de precipitacion organizados de forma trimestral \n")

year_model <- Map(function(x,y) Map(function(x1,y1) years_model=intersect(x1,y1),x, y),years_predictor,years_response)
years_final_res <- Map(function(x,y) Map(function(x1,y1) pos_x=x1%in%y1 ,x,y),years_response,year_model)
years_final_prec=Map(function(x,y) Map(function(x1,y1) pos_x=x1%in%y1 ,x,y),years_predictor,year_model)
data_tsm_final=Map(function(x,y) Map(function(x1,y1) x1[y1,] ,x,y),tsm_list_table,years_final_prec)
data_res_final=Map(function(x,y) Map(function(x1,y1) x1[y1,] ,x,y),prec_cumu,years_final_res)

cat("\n Periodo de entrenamiento generado \n")

lat <-lapply(tsm_raster,function(x1)lapply(x1,function(x) coordinates(x)[,2]))
ponderacion <-lapply(lat,function(x1)lapply(x1,function(x) sqrt(cos((pi/180)*x))))
ponde <- Map(function(x,y)Map(function(x1,y1)x1[y1],x,y),ponderacion,pos)

cat("\n Ponderación PCA \n")

cor_tsm=Map(function(x,y,z)Map(selection_area,x,y,z),data_tsm_final,data_res_final,ponde)

cat("\n Correlación de los pixeles calculada \n")

names_selec <-Map(function(x,y) paste0(x,"/CFSV2_r/",substr(y,1,nchar(y)-4)) ,folders,names_x)
o_empty_1=Map(function(x,y,z,r)Map(files_x,x,y,z,r),tsm_raster,cor_tsm,names_selec,time_sel)

cat("\n Archivos de la TSM construidos por deciles para CPT \n")

path_newx <- lapply(folders,function(x)list.files(paste0(x,"/CFSV2_r"),full.names = T))
new_namesx <- lapply(path_newx,function(x) substr(basename(x),1,nchar(basename(x))-4))
path_run <- Map(function(x,y) paste0(x,"/run/",y,".bat"),folders,new_namesx)
path_output <- Map(function(x,y) paste0(x,"/output/",y),folders,new_namesx)
o_empty_2 <-Map(function(x,y,z,k)Map(run_cpt,x,y,z,k),path_newx,path_y,path_run,path_output)

cat("\n Batch CPT realizado \n")

folder_output <- Map(function(x,y) paste0(x,"/output/",substr(y,1,nchar(y)-4)),folders,names_x)
best_decil_l=lapply(folder_output,function(x)lapply(x,best_GI))
best_decil=lapply(best_decil_l,unlist)

cat("\n Mejor corrida seleccionada \n")

best_path <- Map(function(x,y) paste0(x,"/output/",y),folders,best_decil)
all_metricas<- lapply(best_path,function(x)lapply(x,metricas))
metricas_final=lapply(all_metricas,function(x)do.call("rbind",x))
o_empty_3=Map(function(x,y) write.csv(x,paste0(y,"/metrics.csv"),row.names=FALSE),metricas_final,folders)

cat("\n Metricas de validacion almacenadas \n")

all_prob<- lapply(best_path,function(x)lapply(x,proba))
prob_final <- lapply(all_prob,function(x)do.call("rbind",x))
o_empty_4 <- Map(function(x,y) write.csv(x,paste0(y,"/probabilities.csv"),row.names=FALSE),prob_final,folders)

cat("\n Pronosticos probabilisticos almacenados \n")

path_images <- Map(function(x,y) paste0(x,"/",y),folders,best_decil)
o_empty_5 <- Map(function(x,y,z)Map(save_areas,x,y,z),tsm_raster,cor_tsm,path_images)

cat("\n Pixeles selecionados almacenados en .tiff  \n")
