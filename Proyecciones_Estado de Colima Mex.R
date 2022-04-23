##PROYECCIONES DEL ESTADO DE COLIMA, MEXICO##
##04-12-2019##

rm(list=ls())
gc()

fx5<- read.csv("Fx5_II.CSV",
                   header=TRUE)

fx5q<- read.csv("Fx_quinquenios.CSV",
               header=TRUE)

mx<- read.csv("mx.CSV",
                header=TRUE)

asfr=function(fx5,year){
  
  fx0=fx5[fx5$Year==year,"fx"]
  
  Fx=5*cumsum(fx0)
  TGF=Fx[7]
  FxF=Fx/TGF
  
  x5=seq(17.5,47.5,5)
  Yx= log(-log(FxF))

 
  Yx.lm=lm(Yx[-7] ~ x5[-7])
  a=Yx.lm$coefficients[1]
  b=Yx.lm$coefficients[2]
  
  A=exp(-exp(a))
  B=exp(b)
  x1=c(15:50)
  
  (Fx.estim=TGF*A^(B^(x1)))
  
  fx=Fx.estim[2:36]-Fx.estim[1:35]
  fx=c(Fx.estim[1],fx)
  
  return(fx)
}

fx1=data.frame(matrix(0,36,46))
row.names(fx1)=c(15:50)
names(fx1)=c(1970:2015)


for(i in 1:46){
  
  fx1[,i]=asfr(fx5q,1969+i)
  
}




##Fecundidad##
require(forecast)
require(ggplot2)
require(mvtnorm)


mx<-mx[,-c(1,2)] #tasas centrales de mortalidad
fx <-fx1[1:34,25:45] #fecundidad de los ultimos 30 años 


#migracion
ixt.F<-as.matrix(Inmx[Inmx$Sexo=="Mujeres",-c(1:2)])/
  as.matrix(Nx[Nx$Sexo=="Mujeres",-c(1:2)])
ixt.M<-as.matrix(Inmx[Inmx$Sexo=="Hombres",-c(1:2)])/
  as.matrix(Nx[Nx$Sexo=="Hombres",-c(1:2)])
ext.F<-as.matrix(Emigx[Emigx$Sexo=="Mujeres",-c(1:2)])/
  as.matrix(Nx[Nx$Sexo=="Mujeres",-c(1:2)])
ext.M<-as.matrix(Emigx[Emigx$Sexo=="Hombres",-c(1:2)])/
  as.matrix(Nx[Nx$Sexo=="Hombres",-c(1:2)])

#El 14:66 varía de acuerdo con el numero de años con los que me quiero quedar. Me quedo con el 41:66 porque las demas son edades con puros ceros 
ixt.F<-ixt.F[1:(dim(ixt.F)[1]-20),41:66]
ixt.M<-ixt.M[1:(dim(ixt.M)[1]-20),41:dim(ixt.M)[2]] #41:66 es lo mismo
ext.F<-ext.F[1:(dim(ext.F)[1]-20),41:dim(ext.F)[2]]
ext.M<-ext.M[1:(dim(ext.M)[1]-20),41:dim(ext.M)[2]]



edades <- dim(mx)[1] #dame la dimension y quedate con el numero de renglones = 220
edades.fec <-dim(fx)[1]
tiempo.mort <- dim(mx)[2] #quedo con el número de columnas
añoini.mort <- 1970
añoini.fec <- 1995
tiempo.fec <-dim(fx)[2]
añobase <- 2014
horizonte <- 25
añofin <- añobase+horizonte
tiempo.tot <- tiempo.mort+horizonte

#edades.mig <- dim(ixt.F)[1]
#tiempo.mig <- dim(ixt.F)[2]
#añoini.mig <- 1995


#Inicia la estimacion con el método de Lee Carter (1992)

lc.svd <- function(m,edades,tiempo1, tiempo2,ln){
  if(ln==TRUE){
    lm <- log(m)
  }else{
    lm<-
     m
  }
  
  ax <- rowMeans(lm[,tiempo1:tiempo2])
  lm_a <- lm -ax
  d <- matrix(0,nr = min(edades,tiempo2),
              nc = min(edades,tiempo2))
  
  diag(d) <- svd(lm_a)$d
  
  kt <-(d%*%t(-svd(lm_a)$v))
  bx <- -svd(lm_a)$u
  
 lc.svd=list(ax = ax, bx = bx, kt = kt, D = d)
  
  
}

tabmort <- function(m,edades,sex){
  
  mx <- m
  
  nax <- matrix(0.5,dim(mx)[1],dim(mx)[2])
  ## 1 MUJERES 2 HOMBRES
  if(sex==1){
    for(i in 1:dim(mx)[2]){
      if(mx[1,i]<0.01724){
        nax[1,i] <- 0.14903-2.05527*mx[1,i]
      }else if(mx[1,i]>=0.01724 & mx[1,i]<0.06891){
        nax[1,i] <- 0.04667+3.88089*mx[1,i]
      }else{nax[1,i] <- 0.31411}
    }
  }else{
    for(i in 1:dim(mx)[2]){
      if(mx[1,i]<0.023){
        nax[1,i] <- 0.14929-1.99545*mx[1,i]
      }else if(mx[1,i]>=0.023 & mx[1,i]<0.08307){
        nax[1,i] <- 0.02832+3.26021*mx[1,i]
      }else{nax[1,i] <- 0.29915}
    }
  }
  
  
  nax[edades,] <- 1/mx[edades,]
  
  qx<-matrix(1,dim(mx)[1],dim(mx)[2])
  
  for(i in 1:(dim(mx)[1])){
    qx[i,]<-mx[i,]/(1+(1-nax[i,])*mx[i,])
  }
  
  px <- 1-qx
  
  lx<-matrix(1,dim(mx)[1],dim(mx)[2])
  
  for(i in 2:dim(mx)[1]){
    lx[i,] <- lx[i-1,]*px[i-1,]
  }
  
  dx <- matrix(0,dim(mx)[1],dim(mx)[2])
  dx[dim(mx)[1],] <- lx[dim(mx)[1],]
  for(i in 1:(dim(mx)[1]-1)){
    dx[i,]<-lx[i,]-lx[i+1,]
  }
  
  
  Lx<-matrix(0,dim(mx)[1],dim(mx)[2])
  Lx[1,] <- dx[1,]/mx[1,]
  Lx[edades,] <- dx[edades,]/mx[edades,]
  for(i in 2:(edades-1)){
    Lx[i,]<-(lx[i,]+lx[i+1,])/2
  }
  
  Tx<-matrix(0,dim(mx)[1],dim(mx)[2])
  Tx[edades,]<-Lx[edades,]
  for(i in (edades-1):1){
    Tx[i,]<-Lx[i,]+Tx[i+1,]
  }
  
  ex <- Tx/lx
  
  Sx<-matrix(0,(dim(mx)[1]+1),dim(mx)[2])
  Sx[1,]<-Lx[1,]/lx[1,]
  Sx[(edades+1),] <- Tx[edades,]/Tx[(edades-1),]
  for(i in 2:edades){
    Sx[i,]<-Lx[i,]/Lx[i-1,]
  }
  
  tabmort <- list(Edad=c(0:(edades-1)),mx=mx, nax=nax, qx=qx, 
                  px=px, lx=lx, dx=dx, Lx=Lx, Tx=Tx, ex=ex, Sx=Sx)
}

lc.mort <-lc.svd(mx, edades, tiempo1 = 41, 
                 tiempo2 = tiempo.mort, 
                 ln=TRUE)

lc.fec <-lc.svd(fx, edades=edades.fec, 
                tiempo1 = 10, 
                tiempo2 = tiempo.fec,
                ln=TRUE)

#lc.inmF <-lc.svd(ixt.F, edades=edades.mig, 
                #tiempo1 = 1, 
                #tiempo2 = tiempo.mig,
                #ln=TRUE)

#lc.inmM <-lc.svd(ixt.M, edades=edades.mig, 
                 #tiempo1 = 1, 
                 #tiempo2 = tiempo.mig,
                 #ln=TRUE)

#lc.emigF <-lc.svd(ext.F, edades=edades.mig, 
                 #tiempo1 = 1, 
                 #tiempo2 = tiempo.mig,
                 #ln=TRUE)

#lc.emigM <-lc.svd(ext.M, edades=edades.mig, 
                 #tiempo1 = 1, 
                 #tiempo2 = tiempo.mig,
                 #ln=TRUE)


kt1.fit <-auto.arima(lc.mort$kt[1,], trace=T, d=1)

ft1.fit <-auto.arima(lc.fec$kt[1,], trace=T, d=1
                     ) #lc.fec$kt[1,] para la primera 
it1F.fit <- auto.arima(lc.inmF$kt[1,], trace=T,
                   allowdrift = F)
it1M.fit <- auto.arima(lc.inmM$kt[1,], trace=T,
                   allowdrift = F)
et1F.fit <- auto.arima(lc.emigF$kt[1,], trace=T,
                   allowdrift = F)
et1M.fit <- auto.arima(lc.emigM$kt[1,], trace=T,
                   allowdrift = F)

#proyecciones
h <- 25


kt.for <- forecast(kt1.fit, h = h, c(95))

ft.for <- forecast(ft1.fit, h = h, c(95))

itF.for <- forecast(it1F.fit, h = h, c(95))
itM.for <- forecast(it1M.fit, h = h, c(95))

etF.for <- forecast(et1F.fit, h = h, c(95))
etM.for <- forecast(et1M.fit, h = h, c(95))


#TASAS CENTRALES DE MORTALIDAD
mx.for <- exp(lc.mort$ax + lc.mort$bx[,1]%*%t(kt.for$mean))

#fUNCIONES DE SUPERVIVIENCIA LAS METO A LAS TABLAS DE MORTALIDAD 
SxF.for <- tabmort(mx.for[111:220,], edades = 110, sex=1)$Sx
SxM.for <- tabmort(mx.for[1:110,], edades = 110, sex=2)$Sx

#Tasas especificas de fecundidad
fx.for<-exp(lc.fec$ax + lc.fec$bx[,1]%*%t(ft.for$mean))
colSums(fx.for)           
  
ixF.For <- rbind(exp(lc.inmF$ax + lc.inmF$bx[,1]%*%t(itF.for$mean)),
                 matrix(0,20,35))
         
ixM.For <- rbind(exp(lc.inmM$ax + lc.inmM$bx[,1]%*%t(itM.for$mean)),
                 matrix(0,20,35))

exF.For <- rbind(exp(lc.emigF$ax + lc.emigF$bx[,1]%*%t(etF.for$mean)),
                 matrix(0,20,35))
exM.For <- rbind(exp(lc.emigM$ax + lc.emigM$bx[,1]%*%t(etM.for$mean)),
                 matrix(0,20,35))


#Poblacion 
PxF <- Px[Px$Sexo=="Mujeres", -c(1,2)]
PxM <- Px[Px$Sexo=="Hombres", -c(1,2)]

NxF <- Px[Nx$Sexo=="Mujeres", -c(1,2)]
NxM <- Px[Nx$Sexo=="Hombres", -c(1,2)]

PxF.for <- matrix(0,110,36)
PxM.for <- matrix(0,110,36)

NxF.for <- matrix(0,110,36)
NxM.for <- matrix(0,110,36)

PxF.for[,1] <- PxF[,"2016"]
PxM.for[,1] <- PxM[,"2016"]

NxF.for[,1] <- NxF[,"2015"]
NxM.for[,1] <- NxM[,"2015"]

Bx <- matrix(0,35,35)
BF <- vector(length=35)
BM <- vector(length=35)


#pOBLACION A INICIOS DE AÑO
#PxF.for[2:109,2] <- PxF.for[1:108,1]*(1+0.5*ixF.For[1:108,1])*
 # SxF.for[1:108,1] + PxF.for[2:109,1]*0.5*ixF.For[2:109,1]-
  #PxF.for[1:108,1]*exF.For[1:108,1]


####MUJERES#################
for(i in 2:36){
  
#Con pob a mitad de año 1 a 108
PxF.for[2:109,i] <- (PxF.for[1:108,i-1] + 
  0.5*NxF.for[1:108,i-1]*ixF.For[1:108,i-1])* SxF.for[1:108,i-1]+
  NxF.for[2:109,i-1]*0.5*ixF.For[2:109,i-1]-
  NxF.for[1:108,i-1]*exF.For[1:108,i-1]

#ultimo grupo 109 y mas
PxF.for[110,i]<-(PxF.for[109,i-1] + 
  0.5*NxF.for[109,i-1]*ixF.For[109,i-1])*SxF.for[109,i-1] -
  NxF.for[109, i-1]*exF.For[109,i-1] +
 (PxF.for[110,i-1] + 
  NxF.for[110,i-1]*0.5*ixF.For[110,i-1])*SxF.for[110,i-1]+
    NxF.for[110,i-1]*0.5*ixF.For[110,i-1]-
  NxF.for[110,i-1]*0.5*exF.For[110,i-1]


#Nacimientos 
Bx[,i-1] <- fx.for[,i-1]*(PxF.for[16:50,i-1]+
                        0.5*NxF.for[16:50,i-1]*ixF.For[16:50,i-1] +
                        PxF.for[16:50,i]) * 0.5

BF[i-1] <- (1/2.05)*sum(Bx[,i-1])


#primer grupo de edad

PxF.for[1,i] <- BF[1]*SxF.for[1,i-1] + 
  NxF.for[1,i-1]*0.5*ixF.For[1,i-1]+
  NxF.for[1,i-1]*exF.For[1,i-1]


#POb mitad de año

NxF.for[,i] <- 0.5*(PxF.for[,i-1] + PxF.for[,i])
}
                    
                    matplot(NxF.for, type="l")
                    
    
                    
   #################### hombres #####
            for(i in 2:36){
                      
                     #Con pob a mitad de año 1 a 108
                      PxM.for[2:109,i] <- (PxM.for[1:108,i-1] + 
                                             0.5*NxM.for[1:108,i-1]*ixM.For[1:108,i-1])* SxM.for[1:108,i-1]+
                        NxM.for[2:109,i-1]*0.5*ixM.For[2:109,i-1]-
                        NxM.for[1:108,i-1]*exM.For[1:108,i-1]
                      
                      #ultimo grupo 109 y mas
                      PxM.for[110,i]<-(PxM.for[109,i-1] + 
                                         0.5*NxM.for[109,i-1]*ixM.For[109,i-1])*SxM.for[109,i-1] -
                        NxM.for[109, i-1]*exM.For[109,i-1] +
                        (PxM.for[110,i-1] + 
                           NxM.for[110,i-1]*0.5*ixM.For[110,i-1])*SxM.for[110,i-1]+
                        NxM.for[110,i-1]*0.5*ixM.For[110,i-1]-
                        NxM.for[110,i-1]*0.5*exM.For[110,i-1]
                      
                      
                      #Nacimientos 
                     
                      BM[i-1] <- (1.05/2.05)*sum(Bx[,i-1])
                      
                      
                      #primer grupo de edad
                      
                      PxM.for[1,i] <- BM[1]*SxM.for[1,i-1] + 
                        NxM.for[1,i-1]*0.5*ixM.For[1,i-1]+
                        NxM.for[1,i-1]*exM.For[1,i-1]
                      
                      
                      #POb mitad de año
                      
                      NxM.for[,i] <- 0.5*(PxM.for[,i-1] + PxM.for[,i])
                    }
                    
                    matplot(NxM.for, type="l")
                    
            colSums(NxM.for)
            colSums(NxF.for)       
POB <-   colSums(NxM.for) +      colSums(NxF.for)          

PxF.for <-round(PxF.for,0)
PxM.for <-round(PxM.for,0)
                    