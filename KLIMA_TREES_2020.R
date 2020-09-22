############################################################################################################################################
######## Relationships between tree rings and climate  
######## Based on the Excel routine performed by M. Masiokas (INAIGLA, Argentina) to Excel 
######## This routine was developed by Alvaro Gonzalez Reyes, HEMERA Centro de Observacion de la Tierra, 
######## Universidad Mayor, Chile
######## contact: gonzalezreyesalvaro@gmail.com
############################################################################################################################################


################################# IMPORTANTE ############################################################
############## ANTES DE COMENZAR CARGAR TODO LO QUE CONTIENE ESTE SCRIPT ################################ 
#------------------------------------
# Pegar el siguiente codigo en la consola (sin el signo #). Leera todo el archivo.                                        
# source("KLIMA_TREES_2020.R")
#------------------------------------

##########################################################################################
##### FUNCION PARA TRANSFORMAR LA ESTACION Y CONSIDERAR PERIODOS PREVIOS Y FUTUROS 
##########################################################################################

###### Considerando una estacion meteo que comienza en enero a diciembre. Recordar que la primera columna deben ser los agnos. 

CLIMATE.st <- function(x){
filas = c(nrow(x)+3) ### + 4 porque se va construyendo desde el year prev prev hacia un year futuro, como un escalera hacia arriba. Hacia la izquierda estaran los year prev prev y hacia la derecha estara el year futuro     
year <- as.matrix(c(c(x[1,1]-1):c(x[nrow(x),1]+2)))
data <- x[,1:13]
res <- matrix(NA,ncol=49,nrow=filas)
res[,1]=year[,1]
##### DATA 
for (i in c(2:13)) {res[4:c(filas),i] <- data[,i]} ## ppYear
for (i in c(2:13)) {res[3:c(filas-1),i+12] <- data[,i]} ## pYear 
for (i in c(2:13)) {res[2:c(filas-2),i+24] <- data[,i]} ## Year 
for (i in c(2:13)) {res[1:c(filas-3),i+36] <- data[,i]} ## future Year
res <- data.frame(res)
##### 
m <- c("Ja","Fb","Mr","Ap","My","Jn","Jl","Au","Sp","Oc","Nv","Dc")
meses <- c("year",paste("pp",m,sep=""),paste("p",m,sep=""),m,paste("f",m,sep=""))
colnames(res) <- meses
return(res)
}

############################################################
######## CALCULA LA SIG ESTADISTICA DE UNA CORRELACION 
############################################################

### Se puede utilizar solo, por ejemplo: Sig.correlation(1960,1990)
## t.i = agno inicio
## t.f = agno termino

Sig.correlation <- function(t.i,t.f) {
  t=length(t.i:t.f)
  df=t-2
  t90=qt(1-0.10/2,df)
  t95=qt(1-0.05/2,df)
  t99=qt(1-0.01/2,df)
  t99.9=qt(1-0.001/2,df)
  ### res
  r90 <- t90/sqrt(c(df+t90*t90))
  r95 <- t95/sqrt(c(df+t95*t95))
  r99 <- t99/sqrt(c(df+t99*t99))
  r99.9 <- t99.9/sqrt(c(df+t99.9*t99.9))
  res <- data.frame(critical_r=c("90%","95%","99%","99.9%"),value=c(r90,r95,r99,r99.9))
  return(res)
   }  

ConfLevel <- function(x,n){ # x debe ser la columna de los agnos
CLevel=Sig.correlation(x[1],x[c(length(x) -1)])[,2]
CLevel=CLevel[n]
return(round(CLevel,3))
}

######################################################
######################################################
######################### 
######################### SOLO LEER  
######################### 
######################################################
######################################################

##--------------------------------------------------- 
int <- function(x,delta){
i <- seq(1,c(x-delta))
t <- c(c(1+delta):x)
res <- data.frame(i,t)
return(res)
}


DIS <- function(x,n){
### x vector de datos
#### n es el valor a discriminar
    res <- replace(x, x >= -c(n) & x <= c(n), NA)
    return(x)
}


COMP <- function(chrono,clim){
    x <- chrono ### var 
    y <- clim
    yr.start <- ifelse(x[1,1]>y[1,1],x[1,1],y[1,1]) ### SIEMPRE LA VARIABLE (e.g CRONOLOGIA) DEBE COMENZAR ANTES QUE LA SERIE CLIMATICA 
    yr.end <- ifelse(x[nrow(x),1]< y[nrow(y),1],x[nrow(x),1],y[nrow(y),1])
    r <- c(yr.start,yr.end)
### UBICACION
x <- data.frame((1:nrow(x)),x); xi <- subset(x,x[,2]==r[1])[1] ; xf <- subset(x,x[,2]==r[2])[1]  
y <- data.frame((1:nrow(y)),y); yi <- subset(y,y[,2]==r[1])[1] ; yf <- subset(y,y[,2]==r[2])[1]                     
 return(data.frame(c(xi,xf),c(yi,yf)))
}
##---------------------------------------------------



####### ####### ############## ####### ############## ####### ############## ####### ############## ####### #######
####### ####### ############## ####### ############## ####### ############## ####### ############## ####### #######
####### ####### ############## ####### ############## ####### ############## ####### ############## ####### #######
####### ####### #######     FUNCIONES DE RESPUESTA AL CLIMA 
####### ####### ############## ####### ############## ####### ############## ####### ############## ####### #######
####### ####### ############## ####### ############## ####### ############## ####### ############## ####### #######
####### ####### ############## ####### ############## ####### ############## ####### ############## ####### #######


RESPO.AVG.CLIMATE <- function(var,clim) {
##### Comparacion agnos. DEBE ESTAR CARGADA LA FUNCION COMP 
r <- COMP(var,clim)
var <- var[c(r[[1]]+2):r[[2]],] 
clim <- clim[r[[3]]:r[[4]],]
### Hydro&Climate Station
x=CLIMATE.st(clim)  ### para estacion Enero - Diciembre.
f <- nrow(x)
x=x[4:c(f-2),1:49] ####  CLIMATE DATA CON la col year. en los datos originales, hay que considerar que la estacion tiene los dos year antes, y que llega al mismo year que tiene la variable a correlacionar

##------------------------------------------------------
## Significancia estadistica (por defecto 95%)

time.interval=c(x[1,1]:x[c(nrow(x) -1),1])
CL=ConfLevel(time.interval,2) # 1 = 90%, 2 = 95%, 3 = 99%, 4 = 99.9% 

##------------------------------------------------------

x=x[,2:49] ####  CLIMATE DATA SIN la col year 


y = var[,2] ####  #### CRONOLOGIA o variable (UBICADO EN LA 2DA COLUMNA)


###### ###### ###### ###### ###### ######  IMPORTANTE ###### ###### ###### ###### ###### ###### ###### ######
###### LAS SERIES SERAN CORRELACIONADAS EN un PERIODO COMUN, NORMALMENTE DADO POR LA VARIABLE A CORRELACIONAR 
###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ######

###### DELTA 0
res0 <- cor(x,y,use="pairwise",method="pearson")
res0 <- data.frame(rownames(res0),res0[,1],rep(0,48))


###### DELTA 1 months 
n=1
res1 <- matrix(NA,ncol=nrow(int(48,n)),nrow=1)
res.n1 <- matrix(NA,ncol=nrow(int(48,n)),nrow=1)
for(i in c(1:nrow(int(48,n)))) {res1[,i] <- cor(y,rowMeans(x[,int(48,n)[i,1]:int(48,n)[i,2]],na.rm=TRUE),method="pearson",use="pairwise")}
for(i in c(1:nrow(int(48,n)))) {res.n1[,i] <-paste(colnames(x)[int(48,n)[i,1]],colnames(x)[int(48,n)[i,2]],sep="_")}
res1 <- data.frame(t(res.n1),t(res1),rep(n,nrow(int(48,n))))


###### DELTA 2 months 
n=2
res2 <- matrix(NA,ncol=nrow(int(48,n)),nrow=1)
res.n2 <- matrix(NA,ncol=nrow(int(48,n)),nrow=1)
for(i in c(1:nrow(int(48,n)))) {res2[,i] <- cor(y,rowMeans(x[,int(48,n)[i,1]:int(48,n)[i,2]],na.rm=TRUE),method="pearson",use="pairwise")}
for(i in c(1:nrow(int(48,n)))) {res.n2[,i] <-paste(colnames(x)[int(48,n)[i,1]],colnames(x)[int(48,n)[i,2]],sep="_")}
res2 <- data.frame(t(res.n2),t(res2),rep(n,nrow(int(48,n))))

###### DELTA 3 months 
n=3
res3 <- matrix(NA,ncol=nrow(int(48,n)),nrow=1)
res.n3 <- matrix(NA,ncol=nrow(int(48,n)),nrow=1)
for(i in c(1:nrow(int(48,n)))) {res3[,i] <- cor(y,rowMeans(x[,int(48,n)[i,1]:int(48,n)[i,2]],na.rm=TRUE),method="pearson",use="pairwise")}
for(i in c(1:nrow(int(48,n)))) {res.n3[,i] <-paste(colnames(x)[int(48,n)[i,1]],colnames(x)[int(48,n)[i,2]],sep="_")}
res3 <- data.frame(t(res.n3),t(res3),rep(n,nrow(int(48,n))))

###### DELTA 4 months 
n=4
res4 <- matrix(NA,ncol=nrow(int(48,n)),nrow=1)
res.n4 <- matrix(NA,ncol=nrow(int(48,n)),nrow=1)
for(i in c(1:nrow(int(48,n)))) {res4[,i] <- cor(y,rowMeans(x[,int(48,n)[i,1]:int(48,n)[i,2]],na.rm=TRUE),method="pearson",use="pairwise")}
for(i in c(1:nrow(int(48,n)))) {res.n4[,i] <-paste(colnames(x)[int(48,n)[i,1]],colnames(x)[int(48,n)[i,2]],sep="_")}
res4 <- data.frame(t(res.n4),t(res4),rep(n,nrow(int(48,n))))

###### DELTA 5 months 
n=5
res5 <- matrix(NA,ncol=nrow(int(48,n)),nrow=1)
res.n5 <- matrix(NA,ncol=nrow(int(48,n)),nrow=1)
for(i in c(1:nrow(int(48,n)))) {res5[,i] <- cor(y,rowMeans(x[,int(48,n)[i,1]:int(48,n)[i,2]],na.rm=TRUE),method="pearson",use="pairwise")}
for(i in c(1:nrow(int(48,n)))) {res.n5[,i] <-paste(colnames(x)[int(48,n)[i,1]],colnames(x)[int(48,n)[i,2]],sep="_")}
res5 <- data.frame(t(res.n5),t(res5),rep(n,nrow(int(48,n))))

###### DELTA 6 months 
n=6
res6 <- matrix(NA,ncol=nrow(int(48,n)),nrow=1)
res.n6 <- matrix(NA,ncol=nrow(int(48,n)),nrow=1)
for(i in c(1:nrow(int(48,n)))) {res6[,i] <- cor(y,rowMeans(x[,int(48,n)[i,1]:int(48,n)[i,2]],na.rm=TRUE),method="pearson",use="pairwise")}
for(i in c(1:nrow(int(48,n)))) {res.n6[,i] <-paste(colnames(x)[int(48,n)[i,1]],colnames(x)[int(48,n)[i,2]],sep="_")}
res6 <- data.frame(t(res.n6),t(res6),rep(n,nrow(int(48,n))))

###### DELTA 7 months 
n=7
res7 <- matrix(NA,ncol=nrow(int(48,n)),nrow=1)
res.n7 <- matrix(NA,ncol=nrow(int(48,n)),nrow=1)
for(i in c(1:nrow(int(48,n)))) {res7[,i] <- cor(y,rowMeans(x[,int(48,n)[i,1]:int(48,n)[i,2]],na.rm=TRUE),method="pearson",use="pairwise")}
for(i in c(1:nrow(int(48,n)))) {res.n7[,i] <-paste(colnames(x)[int(48,n)[i,1]],colnames(x)[int(48,n)[i,2]],sep="_")}
res7 <- data.frame(t(res.n7),t(res7),rep(n,nrow(int(48,n))))

###### DELTA 8 months 
n=8
res8 <- matrix(NA,ncol=nrow(int(48,n)),nrow=1)
res.n8 <- matrix(NA,ncol=nrow(int(48,n)),nrow=1)
for(i in c(1:nrow(int(48,n)))) {res8[,i] <- cor(y,rowMeans(x[,int(48,n)[i,1]:int(48,n)[i,2]],na.rm=TRUE),method="pearson",use="pairwise")}
for(i in c(1:nrow(int(48,n)))) {res.n8[,i] <-paste(colnames(x)[int(48,n)[i,1]],colnames(x)[int(48,n)[i,2]],sep="_")}
res8 <- data.frame(t(res.n8),t(res8),rep(n,nrow(int(48,n))))

###### DELTA 9 months 
n=9
res9 <- matrix(NA,ncol=nrow(int(48,n)),nrow=1)
res.n9 <- matrix(NA,ncol=nrow(int(48,n)),nrow=1)
for(i in c(1:nrow(int(48,n)))) {res9[,i] <- cor(y,rowMeans(x[,int(48,n)[i,1]:int(48,n)[i,2]],na.rm=TRUE),method="pearson",use="pairwise")}
for(i in c(1:nrow(int(48,n)))) {res.n9[,i] <-paste(colnames(x)[int(48,n)[i,1]],colnames(x)[int(48,n)[i,2]],sep="_")}
res9 <- data.frame(t(res.n9),t(res9),rep(n,nrow(int(48,n))))

###### DELTA 10 months 
n=10
res10 <- matrix(NA,ncol=nrow(int(48,n)),nrow=1)
res.n10 <- matrix(NA,ncol=nrow(int(48,n)),nrow=1)
for(i in c(1:nrow(int(48,n)))) {res10[,i] <- cor(y,rowMeans(x[,int(48,n)[i,1]:int(48,n)[i,2]],na.rm=TRUE),method="pearson",use="pairwise")}
for(i in c(1:nrow(int(48,n)))) {res.n10[,i] <-paste(colnames(x)[int(48,n)[i,1]],colnames(x)[int(48,n)[i,2]],sep="_")}
res10 <- data.frame(t(res.n10),t(res10),rep(n,nrow(int(48,n))),rep(n,nrow(int(48,n))))


###### DELTA 11 months 
n=11
res11 <- matrix(NA,ncol=nrow(int(48,n)),nrow=1)
res.n11 <- matrix(NA,ncol=nrow(int(48,n)),nrow=1)
for(i in c(1:nrow(int(48,n)))) {res11[,i] <- cor(y,rowMeans(x[,int(48,n)[i,1]:int(48,n)[i,2]],na.rm=TRUE),method="pearson",use="pairwise")}
for(i in c(1:nrow(int(48,n)))) {res.n11[,i] <-paste(colnames(x)[int(48,n)[i,1]],colnames(x)[int(48,n)[i,2]],sep="_")}
res11 <- data.frame(t(res.n11),t(res11),rep(n,nrow(int(48,n))))


###### DELTA 12 months 
n=12
res12 <- matrix(NA,ncol=nrow(int(48,n)),nrow=1)
res.n12 <- matrix(NA,ncol=nrow(int(48,n)),nrow=1)
for(i in c(1:nrow(int(48,n)))) {res12[,i] <- cor(y,rowMeans(x[,int(48,n)[i,1]:int(48,n)[i,2]],na.rm=TRUE),method="pearson",use="pairwise")}
for(i in c(1:nrow(int(48,n)))) {res.n12[,i] <-paste(colnames(x)[int(48,n)[i,1]],colnames(x)[int(48,n)[i,2]],sep="_")}
res12 <- data.frame(t(res.n12),t(res12),rep(n,nrow(int(48,n))))


###### DELTA 13 months 
n=13
res13 <- matrix(NA,ncol=nrow(int(48,n)),nrow=1)
res.n13 <- matrix(NA,ncol=nrow(int(48,n)),nrow=1)
for(i in c(1:nrow(int(48,n)))) {res13[,i] <- cor(y,rowMeans(x[,int(48,n)[i,1]:int(48,n)[i,2]],na.rm=TRUE),method="pearson",use="pairwise")}
for(i in c(1:nrow(int(48,n)))) {res.n13[,i] <-paste(colnames(x)[int(48,n)[i,1]],colnames(x)[int(48,n)[i,2]],sep="_")}
res13 <- data.frame(t(res.n13),t(res13),rep(n,nrow(int(48,n))))

###### DELTA 14 months 
n=14
res14 <- matrix(NA,ncol=nrow(int(48,n)),nrow=1)
res.n14 <- matrix(NA,ncol=nrow(int(48,n)),nrow=1)
for(i in c(1:nrow(int(48,n)))) {res14[,i] <- cor(y,rowMeans(x[,int(48,n)[i,1]:int(48,n)[i,2]],na.rm=TRUE),method="pearson",use="pairwise")}
for(i in c(1:nrow(int(48,n)))) {res.n14[,i] <-paste(colnames(x)[int(48,n)[i,1]],colnames(x)[int(48,n)[i,2]],sep="_")}
res14 <- data.frame(t(res.n14),t(res14),rep(n,nrow(int(48,n))))


###### DELTA 15 months 
n=15
res15 <- matrix(NA,ncol=nrow(int(48,n)),nrow=1)
res.n15 <- matrix(NA,ncol=nrow(int(48,n)),nrow=1)
for(i in c(1:nrow(int(48,n)))) {res15[,i] <- cor(y,rowMeans(x[,int(48,n)[i,1]:int(48,n)[i,2]],na.rm=TRUE),method="pearson",use="pairwise")}
for(i in c(1:nrow(int(48,n)))) {res.n15[,i] <-paste(colnames(x)[int(48,n)[i,1]],colnames(x)[int(48,n)[i,2]],sep="_")}
res15 <- data.frame(t(res.n15),t(res15),rep(n,nrow(int(48,n))))


######## RESULTS MATRIX
values=round(c(res0[,2],res1[,2],res2[,2],res3[,2],res4[,2],res5[,2],res6[,2],res7[,2],res8[,2],res9[,2],res10[,2],res11[,2],res12[,2],res13[,2],res14[,2],res15[,2]),2)

months=c(as.character(res0[,1]),as.character(res1[,1]),as.character(res2[,1]),as.character(res3[,1]),as.character(res4[,1]),as.character(res5[,1]),as.character(res6[,1]),as.character(res7[,1]),as.character(res8[,1]),as.character(res9[,1]),as.character(res10[,1]),as.character(res11[,1]),as.character(res12[,1]),as.character(res13[,1]),as.character(res14[,1]),as.character(res15[,1]))

delta_months=c(res0[,3],res1[,3],res2[,3],res3[,3],res4[,3],res5[,3],res6[,3],res7[,3],res8[,3],res9[,3],res10[,3],res11[,3],res12[,3],res13[,3],res14[,3],res15[,3])

RES <- data.frame(months,values,delta_months)

RES <- data.frame(RES,CL.low.high=rep(CL,nrow(RES)))

##### DISCRIMINA VALORES CON CORRELACIONES CON SIGNIFICANCIA AL P < 0.05

value <- c(RES[1,4]) ### valor de discriminacion 

XX <- replace(RES[,2], RES[,2] >= -c(value) & RES[,2] <= value, NA)

RES <- data.frame(months=RES[,1],values=XX,RES[,3:4])

return(RES)
}


################################################
### ejemplo
################################################

# RESPO.AVG.CLIMATE(tucronologia,tuestacionclimaticaofluviometrica)
# por defecto las correlaciones son discriminadas con un 95% 
# recordar que la estacion climatica o fluviometrica debe considerar siempre dos years antes que la crono 



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
######################## RESPO SUM ##########################
####### UTILIZA LA SUMA ENTRE LOS MESES DE LA ESTACION CLIMATICA.POR EJEMPLO, UTIL PARA PRECIPITACION ACUMULADA. 
#############################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#


RESPO.SUM.CLIMATE <- function(var,clim) {
##### Comparacion agnos. DEBE ESTAR CARGADA LA FUNCION COMP 
r <- COMP(var,clim)
var <- var[c(r[[1]]+2):r[[2]],] 
clim <- clim[r[[3]]:r[[4]],]
### Hydro&Climate Station
x=CLIMATE.st(clim)  ### para estacion Enero - Diciembre.
f <- nrow(x)
x=x[4:c(f-2),1:49] ####  CLIMATE DATA CON la col year. en los datos originales, hay que considerar que la estacion tiene los dos year antes, y que llega al mismo year que tiene la variable a correlacionar

##------------------------------------------------------
## Significancia estadistica (por defecto 95%)

time.interval=c(x[1,1]:x[c(nrow(x) -1),1])
CL=ConfLevel(time.interval,2) # 1 = 90%, 2 = 95%, 3 = 99%, 4 = 99.9% 

##------------------------------------------------------

x=x[,2:49] ####  CLIMATE DATA SIN la col year 


y = var[,2] ####  #### CRONOLOGIA o variable (UBICADO EN LA 2DA COLUMNA)


###### ###### ###### ###### ###### ######  IMPORTANTE ###### ###### ###### ###### ###### ###### ###### ######
###### LAS SERIES SERAN CORRELACIONADAS EN un PERIODO COMUN, NORMALMENTE DADO POR LA VARIABLE A CORRELACIONAR 
###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ######

###### DELTA 0
res0 <- cor(x,y,use="pairwise",method="pearson")
res0 <- data.frame(rownames(res0),res0[,1],rep(0,48))


###### DELTA 1 months 
n=1
res1 <- matrix(NA,ncol=nrow(int(48,n)),nrow=1)
res.n1 <- matrix(NA,ncol=nrow(int(48,n)),nrow=1)
for(i in c(1:nrow(int(48,n)))) {res1[,i] <- cor(y,rowSums(x[,int(48,n)[i,1]:int(48,n)[i,2]],na.rm=TRUE),method="pearson",use="pairwise")}
for(i in c(1:nrow(int(48,n)))) {res.n1[,i] <-paste(colnames(x)[int(48,n)[i,1]],colnames(x)[int(48,n)[i,2]],sep="_")}
res1 <- data.frame(t(res.n1),t(res1),rep(n,nrow(int(48,n))))

###### DELTA 2 months 
n=2
res2 <- matrix(NA,ncol=nrow(int(48,n)),nrow=1)
res.n2 <- matrix(NA,ncol=nrow(int(48,n)),nrow=1)
for(i in c(1:nrow(int(48,n)))) {res2[,i] <- cor(y,rowSums(x[,int(48,n)[i,1]:int(48,n)[i,2]],na.rm=TRUE),method="pearson",use="pairwise")}
for(i in c(1:nrow(int(48,n)))) {res.n2[,i] <-paste(colnames(x)[int(48,n)[i,1]],colnames(x)[int(48,n)[i,2]],sep="_")}
res2 <- data.frame(t(res.n2),t(res2),rep(n,nrow(int(48,n))))

###### DELTA 3 months 
n=3
res3 <- matrix(NA,ncol=nrow(int(48,n)),nrow=1)
res.n3 <- matrix(NA,ncol=nrow(int(48,n)),nrow=1)
for(i in c(1:nrow(int(48,n)))) {res3[,i] <- cor(y,rowSums(x[,int(48,n)[i,1]:int(48,n)[i,2]],na.rm=TRUE),method="pearson",use="pairwise")}
for(i in c(1:nrow(int(48,n)))) {res.n3[,i] <-paste(colnames(x)[int(48,n)[i,1]],colnames(x)[int(48,n)[i,2]],sep="_")}
res3 <- data.frame(t(res.n3),t(res3),rep(n,nrow(int(48,n))))

###### DELTA 4 months 
n=4
res4 <- matrix(NA,ncol=nrow(int(48,n)),nrow=1)
res.n4 <- matrix(NA,ncol=nrow(int(48,n)),nrow=1)
for(i in c(1:nrow(int(48,n)))) {res4[,i] <- cor(y,rowSums(x[,int(48,n)[i,1]:int(48,n)[i,2]],na.rm=TRUE),method="pearson",use="pairwise")}
for(i in c(1:nrow(int(48,n)))) {res.n4[,i] <-paste(colnames(x)[int(48,n)[i,1]],colnames(x)[int(48,n)[i,2]],sep="_")}
res4 <- data.frame(t(res.n4),t(res4),rep(n,nrow(int(48,n))))

###### DELTA 5 months 
n=5
res5 <- matrix(NA,ncol=nrow(int(48,n)),nrow=1)
res.n5 <- matrix(NA,ncol=nrow(int(48,n)),nrow=1)
for(i in c(1:nrow(int(48,n)))) {res5[,i] <- cor(y,rowSums(x[,int(48,n)[i,1]:int(48,n)[i,2]],na.rm=TRUE),method="pearson",use="pairwise")}
for(i in c(1:nrow(int(48,n)))) {res.n5[,i] <-paste(colnames(x)[int(48,n)[i,1]],colnames(x)[int(48,n)[i,2]],sep="_")}
res5 <- data.frame(t(res.n5),t(res5),rep(n,nrow(int(48,n))))

###### DELTA 6 months 
n=6
res6 <- matrix(NA,ncol=nrow(int(48,n)),nrow=1)
res.n6 <- matrix(NA,ncol=nrow(int(48,n)),nrow=1)
for(i in c(1:nrow(int(48,n)))) {res6[,i] <- cor(y,rowSums(x[,int(48,n)[i,1]:int(48,n)[i,2]],na.rm=TRUE),method="pearson",use="pairwise")}
for(i in c(1:nrow(int(48,n)))) {res.n6[,i] <-paste(colnames(x)[int(48,n)[i,1]],colnames(x)[int(48,n)[i,2]],sep="_")}
res6 <- data.frame(t(res.n6),t(res6),rep(n,nrow(int(48,n))))

###### DELTA 7 months 
n=7
res7 <- matrix(NA,ncol=nrow(int(48,n)),nrow=1)
res.n7 <- matrix(NA,ncol=nrow(int(48,n)),nrow=1)
for(i in c(1:nrow(int(48,n)))) {res7[,i] <- cor(y,rowSums(x[,int(48,n)[i,1]:int(48,n)[i,2]],na.rm=TRUE),method="pearson",use="pairwise")}
for(i in c(1:nrow(int(48,n)))) {res.n7[,i] <-paste(colnames(x)[int(48,n)[i,1]],colnames(x)[int(48,n)[i,2]],sep="_")}
res7 <- data.frame(t(res.n7),t(res7),rep(n,nrow(int(48,n))))

###### DELTA 8 months 
n=8
res8 <- matrix(NA,ncol=nrow(int(48,n)),nrow=1)
res.n8 <- matrix(NA,ncol=nrow(int(48,n)),nrow=1)
for(i in c(1:nrow(int(48,n)))) {res8[,i] <- cor(y,rowSums(x[,int(48,n)[i,1]:int(48,n)[i,2]],na.rm=TRUE),method="pearson",use="pairwise")}
for(i in c(1:nrow(int(48,n)))) {res.n8[,i] <-paste(colnames(x)[int(48,n)[i,1]],colnames(x)[int(48,n)[i,2]],sep="_")}
res8 <- data.frame(t(res.n8),t(res8),rep(n,nrow(int(48,n))))

###### DELTA 9 months 
n=9
res9 <- matrix(NA,ncol=nrow(int(48,n)),nrow=1)
res.n9 <- matrix(NA,ncol=nrow(int(48,n)),nrow=1)
for(i in c(1:nrow(int(48,n)))) {res9[,i] <- cor(y,rowSums(x[,int(48,n)[i,1]:int(48,n)[i,2]],na.rm=TRUE),method="pearson",use="pairwise")}
for(i in c(1:nrow(int(48,n)))) {res.n9[,i] <-paste(colnames(x)[int(48,n)[i,1]],colnames(x)[int(48,n)[i,2]],sep="_")}
res9 <- data.frame(t(res.n9),t(res9),rep(n,nrow(int(48,n))))

###### DELTA 10 months 
n=10
res10 <- matrix(NA,ncol=nrow(int(48,n)),nrow=1)
res.n10 <- matrix(NA,ncol=nrow(int(48,n)),nrow=1)
for(i in c(1:nrow(int(48,n)))) {res10[,i] <- cor(y,rowSums(x[,int(48,n)[i,1]:int(48,n)[i,2]],na.rm=TRUE),method="pearson",use="pairwise")}
for(i in c(1:nrow(int(48,n)))) {res.n10[,i] <-paste(colnames(x)[int(48,n)[i,1]],colnames(x)[int(48,n)[i,2]],sep="_")}
res10 <- data.frame(t(res.n10),t(res10),rep(n,nrow(int(48,n))),rep(n,nrow(int(48,n))))


###### DELTA 11 months 
n=11
res11 <- matrix(NA,ncol=nrow(int(48,n)),nrow=1)
res.n11 <- matrix(NA,ncol=nrow(int(48,n)),nrow=1)
for(i in c(1:nrow(int(48,n)))) {res11[,i] <- cor(y,rowSums(x[,int(48,n)[i,1]:int(48,n)[i,2]],na.rm=TRUE),method="pearson",use="pairwise")}
for(i in c(1:nrow(int(48,n)))) {res.n11[,i] <-paste(colnames(x)[int(48,n)[i,1]],colnames(x)[int(48,n)[i,2]],sep="_")}
res11 <- data.frame(t(res.n11),t(res11),rep(n,nrow(int(48,n))))


###### DELTA 12 months 
n=12
res12 <- matrix(NA,ncol=nrow(int(48,n)),nrow=1)
res.n12 <- matrix(NA,ncol=nrow(int(48,n)),nrow=1)
for(i in c(1:nrow(int(48,n)))) {res12[,i] <- cor(y,rowSums(x[,int(48,n)[i,1]:int(48,n)[i,2]],na.rm=TRUE),method="pearson",use="pairwise")}
for(i in c(1:nrow(int(48,n)))) {res.n12[,i] <-paste(colnames(x)[int(48,n)[i,1]],colnames(x)[int(48,n)[i,2]],sep="_")}
res12 <- data.frame(t(res.n12),t(res12),rep(n,nrow(int(48,n))))


###### DELTA 13 months 
n=13
res13 <- matrix(NA,ncol=nrow(int(48,n)),nrow=1)
res.n13 <- matrix(NA,ncol=nrow(int(48,n)),nrow=1)
for(i in c(1:nrow(int(48,n)))) {res13[,i] <- cor(y,rowSums(x[,int(48,n)[i,1]:int(48,n)[i,2]],na.rm=TRUE),method="pearson",use="pairwise")}
for(i in c(1:nrow(int(48,n)))) {res.n13[,i] <-paste(colnames(x)[int(48,n)[i,1]],colnames(x)[int(48,n)[i,2]],sep="_")}
res13 <- data.frame(t(res.n13),t(res13),rep(n,nrow(int(48,n))))

###### DELTA 14 months 
n=14
res14 <- matrix(NA,ncol=nrow(int(48,n)),nrow=1)
res.n14 <- matrix(NA,ncol=nrow(int(48,n)),nrow=1)
for(i in c(1:nrow(int(48,n)))) {res14[,i] <- cor(y,rowSums(x[,int(48,n)[i,1]:int(48,n)[i,2]],na.rm=TRUE),method="pearson",use="pairwise")}
for(i in c(1:nrow(int(48,n)))) {res.n14[,i] <-paste(colnames(x)[int(48,n)[i,1]],colnames(x)[int(48,n)[i,2]],sep="_")}
res14 <- data.frame(t(res.n14),t(res14),rep(n,nrow(int(48,n))))


###### DELTA 15 months 
n=15
res15 <- matrix(NA,ncol=nrow(int(48,n)),nrow=1)
res.n15 <- matrix(NA,ncol=nrow(int(48,n)),nrow=1)
for(i in c(1:nrow(int(48,n)))) {res15[,i] <- cor(y,rowSums(x[,int(48,n)[i,1]:int(48,n)[i,2]],na.rm=TRUE),method="pearson",use="pairwise")}
for(i in c(1:nrow(int(48,n)))) {res.n15[,i] <-paste(colnames(x)[int(48,n)[i,1]],colnames(x)[int(48,n)[i,2]],sep="_")}
res15 <- data.frame(t(res.n15),t(res15),rep(n,nrow(int(48,n))))


######## RESULTS MATRIX

values=round(c(res0[,2],res1[,2],res2[,2],res3[,2],res4[,2],res5[,2],res6[,2],res7[,2],res8[,2],res9[,2],res10[,2],res11[,2],res12[,2],res13[,2],res14[,2],res15[,2]),2)

months=c(as.character(res0[,1]),as.character(res1[,1]),as.character(res2[,1]),as.character(res3[,1]),as.character(res4[,1]),as.character(res5[,1]),as.character(res6[,1]),as.character(res7[,1]),as.character(res8[,1]),as.character(res9[,1]),as.character(res10[,1]),as.character(res11[,1]),as.character(res12[,1]),as.character(res13[,1]),as.character(res14[,1]),as.character(res15[,1]))

delta_months=c(res0[,3],res1[,3],res2[,3],res3[,3],res4[,3],res5[,3],res6[,3],res7[,3],res8[,3],res9[,3],res10[,3],res11[,3],res12[,3],res13[,3],res14[,3],res15[,3])

RES <- data.frame(months,values,delta_months)

RES <- data.frame(RES,CL.low.high=rep(CL,nrow(RES)))

##### DISCRIMINA VALORES CON CORRELACIONES CON SIGNIFICANCIA AL P < 0.05

value <- c(RES[1,4]) ### valor de discriminacion 

XX <- replace(RES[,2], RES[,2] >= -c(value) & RES[,2] <= value, NA)

RES <- data.frame(months=RES[,1],values=XX,RES[,3:4])

return(RES)
 }


##################################################################################################
####### CORRRELACION ENTRE GRUPO DE CRONOLOGIAS Y SERIE CLIMATICA. NO EXISTE LIMITE DE CRONOLOGIAS A PONER 
#### ###########################################################################################
##################################################################################################




CRONOS.CORREL.AVG.CLIM <- function(var,clim){
xx <- clim #### ESTACION HYDROMETEO
yy <- var  #### e.g. Cronologias 
c <- ncol(yy)
yy.n=colnames(yy)[2:c]
M <- matrix(NA,nrow=648,ncol=c)
for (i in c(2:c)) {M[,i] <- RESPO.AVG.CLIMATE(yy[,c(1,i)],xx)[,2]}
M <- data.frame(M)
M[,1] <- as.character(RESPO.AVG.CLIMATE(yy[,c(1:2)],xx)[,1]) ### Agrega nombre de meses primera columna 
data <- RESPO.AVG.CLIMATE(yy[,c(1:2)],xx)[,3:4]
### names
colnames(M)=c("Months",yy.n)
M <- cbind(M,data)
return(M)
 }


CRONOS.CORREL.SUM.CLIM <- function(var,clim){
xx <- clim #### ESTACION HYDROMETEO
yy <- var  #### e.g. Cronologias 
c <- ncol(yy)
yy.n=colnames(yy)[2:c]
M <- matrix(NA,nrow=648,ncol=c)
for (i in c(2:c)) {M[,i] <- RESPO.SUM.CLIMATE(yy[,c(1,i)],xx)[,2]}
M <- data.frame(M)
M[,1] <- as.character(RESPO.SUM.CLIMATE(yy[,c(1:2)],xx)[,1]) ### Agrega nombre de meses primera columna 
data <- RESPO.SUM.CLIMATE(yy[,c(1:2)],xx)[,3:4]
### names
colnames(M)=c("Months",yy.n)
M <- cbind(M,data)
return(M)
 }

#------------------------------------------------------------ IMPORTANTE--------------------------------------------------------------------------------------------------
# por defecto, las correlaciones son significantes al 95%. En caso de que se quiera modificar eso, modificar el CL en las funciones RESPO.SUM.CLIMATE y RESPO.AVG.CLIMATE, luego releer el script completo. 
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------


################################################################
################ EJEMPLO
################ BASADO EN LOS DATOS DE OSCILACION ANTARTICA AAO NCAR NCEP y la PP DE LA CIUDAD DE VALDIVIA ESTACION ISLA TEJA UACh
################ el signo " = " es identico a " <- "

#### PARA LEER ARCHIVO EXCEL

# install.packages("gdata") ## instalar paquete gdata
# library(gdata) ## activar paquete
# EN CASO CONTRARIO SE PUEDEN TRANSFORMAR A .txt o .csv. IMPORTANTE: Hay diferencias entre la creacion de estos datos en UNIX y Windows. 


## DATOS EN .TXT
# AAO_NCAR_NCEP = read.table("AAO_NCAR_NCEP_1948_2011.txt") ### COMIENZA EN 1948 e incluye las medias estaciones de verano, primavera, otogno, invierno 
# PP_VAL = read.table("pp_VAL_Jan_to_Dec.txt") ### COMIENZA EN 1960 


##################################################################
###### CORRE FUNCION ###############################################
##################################################################
### RELACION ENTRE EL AAO Y LA PRECIPITACION de ABRIL DE VALDIVIA. ANALOGO PARA UNA CRONOLOGIA DE ANCHO DE ANILLOS 
###
### IMPORTANTE
### 
########################

##(Est climatica/ Variable)

# RESPO.AVG.CLIMATE(PP_VAL[,c(1,5)],AAO_NCAR_NCEP)

#### SE UTILIZA EL PARENTESIS [ ], que indican [filas,columnas], para seleccionar un intervalo de filas. Se debe ocupar ":" para representar un intervalo continuo. Por ejemplo: 2:10. La variable debe contener en la primera COLUMNA los years.  

### LA FUNCION ESTA CONFIGURADA SOLO PARA ENTREGAR VALORES ESTADISTICAMENTE SIGNIFICANTES A UN 95% DE CONFIANZA



#### PARA UN GRUPO 
## CRONOS.CORREL.AVG.CLIM(PP_VAL,AAO_NCAR_NCEP)
### ESTA FUNCION HACE LO MISMO QUE LA DE ARRIBA, PERO AHORA TOMANDO UN GRUPO DE DATOS (CRONOLOGIAS, MESES DEL AGNO CORRESPONDIENTES A UNA ESTACION METEO, etc.....).#### IMPORTANTE: La primera columna en el grupo de datos deben ser los agnos
### IMPORTANTE: La matriz de cronologias debe considerar el periodo comun para todas las cronologias 



################################################################################################################
################
################ GUARDAR RESULTADOS 
################
################################################################################################################

### ocupar la funcion write.table (por default en R)
### OPCION: instalar paquete para guardar en archivo excel directo ("xlsx")

# install.packages("xlsx")
# write.xlsx(RESPO.AVG.CLIMATE(PP_VAL[,c(1,5)],AAO_NCAR_NCEP),"RES.xlsx",row.names = FALSE)
























