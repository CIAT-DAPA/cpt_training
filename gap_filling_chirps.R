library(raster)
require(ncdf4)

chirps_all = stack("C:\\Users\\lllanos\\Downloads\\data_peru.nc")

setwd("C:/Users/lllanos/Desktop/peru")
dpto="chirps"
dir.create(dpto)

station_data = read.csv(file = paste0("PP_C_N.csv"),header=T,na.strings = "-9999")
station_data = station_data[station_data$year>1980,]
station_coord = read.csv(file = "coordenadas.csv",header=T)

station_chirps.b = raster::extract(x=chirps_all, y=station_coord[,-1], method = 'bilinear')
station_chirps.b = as.data.frame(t(station_chirps.b))[1:nrow(station_data),]
names(station_chirps.b)=names(station_data)[-2:-1]

dates=seq(as.Date("1981/01/01"),as.Date("2016/03/31"),"month")
months=months.Date(dates)
names_st=names(station_chirps.b)


add_legend <- function(...) {
  opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), 
              mar=c(0, 0, 0, 0), new=TRUE)
  on.exit(par(opar))
  plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
  legend(...)
}

setwd("C:/Users/lllanos/Desktop/peru/chirps")

for (i in 1:ncol(station_chirps.b)){
  data.model = as.data.frame(cbind("y"=station_data[,i+2],"x"=station_chirps.b[,i]))
  model = lm(data=data.model,formula = y~x)
  rmse <- round(sqrt(mean(resid(model)^2)), 2)
  coefs <- coef(model)
  b0 <- round(coefs[1], 2)
  b1 <- round(coefs[2],2)
  r2 <- round(summary(model)$r.squared, 2)
  
  eqn <- bquote(italic(y) == .(b0) + .(b1)*italic(x) * "," ~~ 
                  R^2 == .(r2) * "," ~~ RMSE == .(rmse))
  
  to_predict = as.data.frame(station_chirps.b[,i])
  names(to_predict)="x"
  data_model = predict(model,to_predict)
  data_model[data_model<0] = 0
  tiff(paste0(names_st[i],".tiff"),compression = 'lzw',height = 10,width = 10,units="in", res=200)
  par(mfrow=c(2,1))
  
  plot(dates,station_data[,i+2],lwd=1.5,type="l",xlab="",ylab="Precipitation (mm)",main=names_st[i])
  lines(dates,station_chirps.b[,i],col="red",lty=2,lwd=1)
  lines(dates,data_model,col="blue",lty=2)
  
  
  plot(station_data[,i+2],station_chirps.b[,i],xlab="Observed_stations",ylab="CHIRPS")
  abline(model,col="red")
  legend('bottomright', legend = eqn, bty = 'n')
  
  add_legend("topright",c("Observed","CHIRPS","Model"),
             horiz=T, bty='n', cex=0.9,lty=c(1,2,2),lwd=c(1.5,1,1),col=c("black","red","blue")) 
  
  dev.off()
  
  pos.na = which(is.na(station_data[,i+2]))
  station_data[pos.na,i+2] = as.numeric(data_model[pos.na])
  
}

write.csv(station_data,paste0(dpto,"_fill_prec.csv"),row.names = F)


