######################################################################################
################ PREPROCESSING & GRADUATION ANALYSES #################################
######################################################################################

library(easypackages)
#Installing packages
#packages("KernSmooth", "splines", "gnm", "mgcv", "demography", "StMoMo", "tseries",
 #       "rainbow", "MortalityLaws") 

packages("KernSmooth", "splines", "plotly", "rainbow", "demography","hts", "fitdistrplus") 
#Loading packages
libraries("KernSmooth", "splines", "plotly", "rainbow", "demography", "hts", "MortalityLaws", "fitdistrplus")

#install.packages("webshot")
#webshot::install_phantomjs()

options(scipen = 999)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path ))

#availableLaws()

#Clearing workspace and loading data
rm(list=ls())
load("../../../R Output/Data Prepping/Data.Data")
standard_table<-read.table("../../../standard_table.txt",header=T)
source("myfunctions.R")


##################################################################
########### SMOOTH SPLINE GRADUATION (AGGREGATE) #################
##################################################################


par(mfrow=c(1,2))
plot(standard_table$age,log((-log(1-(standard_table$female)))),
     main="Kenyan Abridged (1990 - 2017) & Standard (2010) Mortality ",
     ylab="log(mu.hat))",xlab="age",cex.lab=0.7,cex.main=0.7,xlim=c(0,100),type="l",
     ylim=c(-9,0), cex.axis=0.7)
lines(standard_table$age,log((-log(1-(standard_table$male)))),col=2)
for(i in 1:28)
{
  start <- which(AMx$year == (1989+i))[(1)]
  end <- which(AMx$year == (1989+i))[(21)]
  lines(AMx$age[start:end],log(AMx$female[start:end]),col=1)
  lines(AMx$age[start:end],log(AMx$male[start:end]),col=2) 
}

legend(0,0,c("Female","Male"),col=c(1:2),lty=rep(1,2),cex=0.6)
text(40,-3,"Arbridged mu_x Rates",cex=0.6)
text(80,-6,"Standard mu_x Rates",cex=0.6)

plot(standard_table$age,log((-log(1-(standard_table$female)))),
     main="Kenyan Graduated (1990 - 2017) & Standard (2010) Mortality ",
     ylab="log(mu.hat))",xlab="age",cex.lab=0.7,cex.main=0.7,xlim=c(0,100),type="l",
     ylim=c(-9,0), cex.axis=0.7)
lines(standard_table$age,log((-log(1-(standard_table$male)))),col=2)
df1<-df2<-df3<-0
Amort.m <- Amort.f <- matrix(0,nrow=28,ncol=91)
for(i in 1:28)
{
  start <- which(AMx$year == (1989+i))[(1)]
  end <- which(AMx$year == (1989+i))[(21)]
  x<-AMx$age[start:end]
  y<-log(AMx$female[start:end])
  z<-log(AMx$male[start:end])
  fit3<-smooth.spline(x,y,cv = TRUE)
  fit4<-smooth.spline(x,z,cv = TRUE)
  df1[i]<-fit3$df
  df2[i]<-fit4$df
  new1<-predict(fit3,data.frame("x"=seq(0,90,1)))
  new2<-predict(fit4,data.frame("x"=seq(0,90,1)))
  lines(seq(0,90,1),(new1$y)$x,col=1)
  lines(seq(0,90,1),(new2$y)$x,col=2)
  Amort.m[i,]<-(new2$y)$x
  Amort.f[i,]<-(new1$y)$x
}

legend(0,0,c("Female","Male"),col=c(1:2),lty=rep(1,2),cex=0.6)
text(40,-3,"Graduated mu_x Rates",cex=0.6)
text(80,-6,"Standard mu_x Rates",cex=0.6)
df<-data.frame(female=df1,male=df2)
#df
#head(Amort.f)



##################################################################
########### CAUSE OF DEATH SMOOTH SPLINE GRADUATION ##############
##################################################################

CODmort.f<-CODmort.t<-CODmort.m<-matrix(0,nrow=28,ncol=91)
List1<-List2<-List3<-list()
for(j in 3:12)
{
  for (i in 1:28)
  {
    start <- which(CODMx.F$year == (1989+i))[(1)]
    end <- which(CODMx.F$year == (1989+i))[(21)]
    x<-CODMx.F$age[start:end]
    y<-log(CODMx.F[start:end,j]); y[!is.finite(y)] <- 0
    z<-log(CODMx.T[start:end,j]); z[!is.finite(z)] <- 0
    t<-log(CODMx.M[start:end,j]); t[!is.finite(t)] <- 0
    fit6<-smooth.spline(x,y,cv = TRUE)
    fit7<-smooth.spline(x,z,cv = TRUE)
    fit8<-smooth.spline(x,t,cv = TRUE)
    new4<-predict(fit6,data.frame("x"=seq(0,90,1)))
    new5<-predict(fit7,data.frame("x"=seq(0,90,1)))
    new6<-predict(fit8,data.frame("x"=seq(0,90,1)))
    CODmort.f[i,]<-(new4$y)$x
    CODmort.t[i,]<-(new5$y)$x
    CODmort.m[i,]<-(new6$y)$x
  }
  List1[[j-2]]<-CODmort.f
  List2[[j-2]]<-CODmort.t
  List3[[j-2]]<-CODmort.m
}

#head(CODmort.t)

#List1
#List2

##################################################################
#########   PREPPING GRADUATED DATA    ###########################
##################################################################

age.new<-rep(0:90,28)
year.new<-rep(1990:2017,each =91)

CODmux.F<-data.frame(year.new,age.new,exp(c(t(List1[[1]]))),exp(c(t(List1[[2]]))),
                     exp(c(t(List1[[3]]))),exp(c(t(List1[[4]]))),exp(c(t(List1[[5]]))),
                     exp(c(t(List1[[6]]))),exp(c(t(List1[[7]]))),
                     exp(c(t(List1[[8]]))),exp(c(t(List1[[9]]))),exp(c(t(List1[[10]]))))
                     
colnames(CODmux.F)<-names(CODMx.F)

CODmux.T<-data.frame(year.new,age.new,exp(c(t(List2[[1]]))),exp(c(t(List2[[2]]))),
                     exp(c(t(List2[[3]]))),exp(c(t(List2[[4]]))),exp(c(t(List2[[5]]))),
                     exp(c(t(List2[[6]]))),exp(c(t(List2[[7]]))),
                     exp(c(t(List2[[8]]))),exp(c(t(List2[[9]]))),exp(c(t(List2[[10]]))))

colnames(CODmux.T)<-names(CODMx.T)

CODmux.M<-data.frame(year.new,age.new,exp(c(t(List3[[1]]))),exp(c(t(List3[[2]]))),
                     exp(c(t(List3[[3]]))),exp(c(t(List3[[4]]))),exp(c(t(List3[[5]]))),
                     exp(c(t(List3[[6]]))),exp(c(t(List3[[7]]))),
                     exp(c(t(List3[[8]]))),exp(c(t(List3[[9]]))),exp(c(t(List3[[10]]))))

colnames(CODmux.M)<-names(CODMx.M)

Amort.f.new<-c(matrix(c(Amort.f),nrow=length(Amort.f[1,]),byrow=T))
Amort.m.new<-c(matrix(c(Amort.m),nrow=length(Amort.m[1,]),byrow=T))

ex<-rbind(cbind(year=1990,age=0,Ex1[1,3:5]+Ex1[2,3:5]),Ex1[3:92,])
for(i in 2:28)
{
  start <- which(Ex1$year == (1989+i))[(1)]
  end <- which(Ex1$year == (1989+i))[(92)]
  ex<-rbind(ex,rbind(cbind(year=(1989+i),age=0,Ex1[start,3:5]+Ex1[(start+1),3:5]),
                     Ex1[(start+2):end,]))
}
#ex

Adx<-data.frame(year=year.new,age=age.new,female=exp(Amort.f.new)*ex[,3],
                male=exp(Amort.m.new)*ex[,4],total=(exp(Amort.f.new)*ex[,3]) + (exp(Amort.m.new)*ex[,4]) )
Amux<-data.frame(year=year.new,age=age.new,female=Adx[,3]/ex[,3],
                 male=Adx[,4]/ex[,4], total=Adx[,5]/ex[,5])

CODdx.F<-data.frame(year.new,age.new,CODmux.F[,3:12]*ex[,3])
colnames(CODdx.F)<-names(CODmux.F)
head(CODdx.F)

CODdx.M<-data.frame(year.new,age.new,CODmux.M[,3:12]*ex[,4])
colnames(CODdx.M)<-names(CODmux.M)
head(CODdx.M)

CODdx.T<-data.frame(year.new,age.new,CODmux.T[,3:12]*ex[,5])
colnames(CODdx.T)<-names(CODmux.T)
head(CODdx.T)

Adxm <- matrix(Adx$male,ncol=28); exm <- matrix(ex$male, ncol=28)
colnames(Adxm) <- 1990:2017; colnames(exm) <- 1990:2017
#head(Adxm)

#Extending to Old Age using predefined functions

pfcasts.f <- pfcasts.m <- newcast.f <- newcast.m <- list()
ex.new.f <- ex.new.m <- myCOD.f <- myCOD.m <- list()
ex.f <- matrix(ex[,3],ncol=28); ex.m <- matrix(ex[,4],ncol=28)
ex1.f <- ex1.m <- 0

#pfcasts are the forecasts based on the linear model, while newcasts
#are the reconciled forecasts
for(i in 1:28)
{
  mykt.f <- hts(CODdx.F[which(CODdx.F$year==(1989+i)),3:12])
  mykt.m <- hts(CODdx.M[which(CODdx.M$year==(1989+i)),3:12])
  y <- aggts(mykt.f);  z <- aggts(mykt.m)
  y.f <- 81 ; y.m <- 81
  loopout.f <- lapply(y, function (x) me(x)); loopout.m <- lapply(z, function (x) me(x))
  #Calculating the base forecasts using exponential decay
  pfcasts.f [[i]] <- sapply(loopout.f, function(x) x$pfcasts)
  pfcasts.m [[i]] <- sapply(loopout.m, function(x) x$pfcasts)
  #Adjusting the first 3 forecasts >- distributing the under/over estimated deaths 
  #equally over the first 3 forecasted years i.e ages 91, 92 & 93
  pfcasts.f[[i]][,1] <- c(pfcasts.f[[i]][1] + (ex.f[91,i]-sum(pfcasts.f[[i]][2:14,1]) - pfcasts.f[[i]][1])/3, 
                        pfcasts.f[[i]][2] + (ex.f[91,i]-sum(pfcasts.f[[i]][2:14,1]) - pfcasts.f[[i]][1])/3,
                        pfcasts.f[[i]][3] + (ex.f[91,i]-sum(pfcasts.f[[i]][2:14,1]) - pfcasts.f[[i]][1])/3,
                        pfcasts.f[[i]][-c(1:3),1])
  pfcasts.m[[i]][,1] <- c(pfcasts.m[[i]][1] + (ex.m[91,i]-sum(pfcasts.m[[i]][2:14,1]) - pfcasts.m[[i]][1])/3, 
                        pfcasts.m[[i]][2] + (ex.m[91,i]-sum(pfcasts.m[[i]][2:14,1]) - pfcasts.m[[i]][1])/3,
                        pfcasts.m[[i]][3] + (ex.m[91,i]-sum(pfcasts.m[[i]][2:14,1]) - pfcasts.m[[i]][1])/3,
                        pfcasts.m[[i]][-c(1:3),1])
  newcast.f [[i]] <- TdFp(pfcasts.f[[i]], mykt.f$nodes)
  newcast.m [[i]] <- TdFp(pfcasts.m[[i]], mykt.m$nodes)
  
  for (j in 1: 14) 
  {
    ex1.f [j] <- sum(rowSums(newcast.f[[i]])[j:14])
    ex1.m [j] <- sum(rowSums(newcast.m[[i]])[j:14])
  }
  ex.new.f[[i]] <- ex1.f ; ex.new.m[[i]] <- ex1.m
  
  myCOD.f [[i]] <- rbind(CODdx.F[which(CODdx.F$year==(1989+i)),3:12], newcast.f[[i]])
  myCOD.m [[i]] <- rbind(CODdx.M[which(CODdx.M$year==(1989+i)),3:12], newcast.m[[i]])
}

#pfcasts.m[[28]]
#ex.new.m[[28]]
#newcast.m[[28]]
#myCOD[[27]]


#meee.f <- myCOD.f[[1]]; meee.m <- myCOD.m[[1]]; 
ex2.f <- ex2.m <- matrix(0,ncol=28,nrow=105)
ex2.f [1:91,] <- ex.f ; ex2.m [1:91,] <- ex.m ; 
ex2.f[92:105,1] <- ex.new.f[[1]]; ex2.m[92:105,1] <- ex.new.m[[1]]
for (a in 2:28)
{
  #meee.f <- rbind(meee.f, myCOD.f[[a]]); meee.m <- rbind(meee.m, myCOD.m[[a]])
  ex2.f[92:105,a] <- ex.new.f[[a]]; ex2.m[92:105,a] <- ex.new.f[[a]]
}

ex <- data.frame(year = rep(1990:2017,each=105), age = rep(0:104,28), female = c(ex2.f),
                 male = c(ex2.m), total = c(ex2.f) + c(ex2.m))

my_kannisto <- function(x, year, age, age.fit)
{
  mux <- list()
  mymx <- matrix(x, ncol=28)
  colnames(mymx)<-year
  rownames(mymx)<-age
  M2 <- MortalityLaw(x = age[age.fit+1], mx = mymx[aget.fit+1,], law = method[i])
}

mymx<-matrix(Amux[,5],ncol=28)
year<-1990:2017
age<-Amux[which(Amux$year==1990),2]
colnames(mymx)<-year
rownames(mymx)<-age

method<-c("gompertz", "makeham", "perks", "beard_makeham", "ggompertz", "HP4",
          "kostaki", "kannisto", "kannisto_makeham")

age.fit <- 61:76
age.fit+1
mux<-list()
par(mfrow=c(3,3))
for (i in 1:length(method))
{
  M2 <- MortalityLaw(x = age[61:76], mx = mymx[61:76,], law = method[i])
  #M2
  #fitted(M2)
  mux[[i]]<-predict(M2, x = 60:110) 
  
  #Rainbow Plots
  ty<-fts(60:110, log(mux[[i]]),xname="age",yname="log.mux")
  #plot.fds(ty, plot.type = "functions",main=paste(method[i], " Model"),
  #         cex.lab=0.7,cex.main=0.7,cex.axis=0.7)
}





extrapolate <- function(x,y)
{
  agg <- list()
  agg$d95 <- colSums(matrix(x,ncol=28)[91:95,])
  agg$d100 <- colSums(matrix(x,ncol=28)[96:100,])
  agg$d105 <- colSums(matrix(x,ncol=28)[101:105,])
  #agg$comp <- c(rbind(matrix(y,ncol=28), colSums(matrix(x,ncol=28)[91:95,]),
  #colSums(matrix(x,ncol=28)[96:100,]), colSums(matrix(x,ncol=28)[101:105,])))
  return(agg)
}

#tail(matrix( myCOD.f[,3],ncol=28))

myloop.f <- lapply(CODdx.F[,3:12], function (x) extrapolate(x))
myloop.m <- lapply(CODdx.M[,3:12], function (x) extrapolate(x))

#sapply(myloop.f, function(x) x)
d95.f <- sapply(myloop.f, function(x) x$d95); d95.m <- sapply(myloop.m, function(x) x$d95)
d100.f <- sapply(myloop.f, function(x) x$d100); d100.m <- sapply(myloop.m, function(x) x$d100)
d105.f <- sapply(myloop.f, function(x) x$d105); d105.m <- sapply(myloop.m, function(x) x$d105)

#Combining 
#Female
CODDx.F.1 <- matrix(0,nrow=672,ncol=10)
for (i in 3:9)
{
  CODDx.F.1[,i] <- c(rbind(matrix(CODDx.F[,i],ncol=28), d95.f[,(i-2)], d100.f[,(i-2)], d105.f[,(i-2)]))
}

CODDx.F1 <- data.frame(rep(1990:2017,each=24), rep(c(unique(CODDx.F$age), 95,100,105),28), CODDx.F.1[,3:9])
names(CODDx.F1) <- colnames(CODDx.F)


#Male
CODDx.M.1 <- matrix(0,nrow=672,ncol=10)
for (i in 3:10)
{
  CODDx.M.1[,i] <- c(rbind(matrix(CODDx.M[,i],ncol=28), d95.m[,(i-2)], d100.m[,(i-2)], d105.m[,(i-2)]))
}

CODDx.M1 <- data.frame(rep(1990:2017,each=24), rep(c(unique(CODDx.M$age), 95,100,105),28), CODDx.M.1[,3:10])
names(CODDx.M1) <- colnames(CODDx.M)

myEX.f <- myEX.m <- matrix(0,ncol=28,nrow=14)
for(i in 1:28)
{
  myEX.f[,i] <- ex.new.f[[i]]
  myEX.m[,i] <- ex.new.m[[i]]
}

myEX1.f <- rbind(matrix(Ex[,3],ncol=28), ex.f[91,]+colSums(myEX.f[1:4,]), 
                 colSums(myEX.f[5:9,]), colSums(myEX.f[10:14,]))
myEX1.m <- rbind(matrix(Ex[,4],ncol=28), ex.m[91,]+colSums(myEX.m[1:4,]), 
                 colSums(myEX.m[5:9,]), colSums(myEX.m[10:14,]))

Ex.new <- data.frame(rep(1990:2017,each=24), rep(c(unique(CODDx.M$age), 95,100,105),28),
                c(myEX1.f), c(myEX1.m), c(myEX1.f) + c(myEX1.m))
names(Ex.new) <- colnames(Ex)


################################################################################
################   GRADUATION OF THE EXTRAPOLATED RATES   ######################
################################################################################

#Getting complete COD data (i.e without zeros for any age) by going
#through each row and determine if a value is zero

colsub.f1 = apply(CODDx.F1, 2, function(col) all(col !=0 ))
colsub.m1 = apply(CODDx.M1, 2, function(col) all(col !=0 ))
colsub.f1
colsub.m1


#Maternal and neonatal mortality will have to be separately graduated
#Getting complete data
CODDx.Fn<-CODDx.F1[ ,colsub.f1]
CODDx.Mn<-CODDx.M1[,colsub.m1]

head(CODDx.Fn)
head(CODDx.Mn)
#CODDx.Fn[,3]/Ex.new[,3]
#morty <- matrix(0,nrow=28,ncol=105)
#for(i in 1:28)
#{
#  x <-c(7/365,28/365,1,seq(5,105,5))
 # y <- log(matrix(CODDx.Fn[,3]/Ex.new[,3],ncol=28)[,i])
#  fit <- smooth.spline(x,y,cv = TRUE)
#  new <- predict(fit,data.frame("x"=seq(0,104,1)))
#  morty[i,]<-(new$y)$x
#}
#c(matrix(exp(c(morty)),ncol=28,byrow=T))


#matrix(CODDx.Fn[,3]/Ex.new[,3],ncol=28)[,1]
mygrad <- function(x, z)
{
  myList <- list()
  morty <- matrix(0,nrow=28,ncol=105)
  for(i in 1:28)
  {
    k <-c(7/365,28/365,1,seq(5,105,5))
    y <- log(matrix(x/z,ncol=28)[,i])
    fit <- smooth.spline(k,y,cv = TRUE)
    new <- predict(fit,data.frame("k"=seq(0,104,1)))
    morty[i,]<-(new$y)$k
  }
  myList$mort <- c(matrix(exp(c(morty)),ncol=28,byrow=T))
  return(myList)
}

myloopf <- lapply(CODDx.Fn[,3:8], function (x, z=Ex.new[,3]) mygrad(x, z=Ex.new[,3]))
myloopm <- lapply(CODDx.Mn[,3:10], function (x, z=Ex.new[,4]) mygrad(x, z=Ex.new[,4]))
mortf105 <- sapply(myloopf, function(x) x$mort); mortm105 <- sapply(myloopm, function(x) x$mort)

#For maternal and neonatal
head(CODDx.F1,15)
CODmort.fm1<-matrix(0,nrow=28,ncol=51)
for (i in 1:28)
{
  start <- which(CODDx.F1$year == (1989+i))[(1)];
  end <- which(CODDx.F1$year == (1989+i))[(13)]
  x <- CODDx.F1$age[start:end]
  y <- log(CODDx.F1[start:end,8]/Ex.new[start:end,3])
  y[!is.finite(y)] <- 0
  fitme<-smooth.spline(x,y,cv = TRUE)
  newe<-predict(fitme,data.frame("x"=seq(0,50,1)))
  CODmort.fm1[i,]<-(newe$y)$x
}


#PREPPING THE GRADUATED DATA
head(CODDx.F1)
head(mortf105)

age.new1<-rep(0:104,28)
year.new1<-rep(1990:2017,each =105)
CODmux.F <- data.frame(year.new1, age.new1, mortf105[,1:5],
                        c(matrix(c(cbind(exp(CODmort.fm1),matrix(0,nrow=28,ncol=54))),
                                 ncol=28,byrow=T)),mortf105[,6])
colnames(CODmux.F)<-names(CODDx.F1)
CODmux.M <- data.frame(year.new1, age.new1, mortm105)
colnames(CODmux.M)<-names(CODDx.M1)

Adx <- data.frame(year=year.new1,age=age.new1,female=rowSums(CODdx.F1[,3:9]),
                  male=rowSums(CODdx.M1[,3:10]),total=rowSums(CODdx.F1[,3:9]) + rowSums(CODdx.M1[,3:10]) )
ex <- data.frame(year=year.new1,age=age.new1,female=rowSums(CODdx.F1[,3:9])/rowSums(CODmux.F[,3:9]),
                 male=rowSums(CODdx.M1[,3:10])/rowSums(CODmux.M[,3:10]),
                 total=rowSums(CODdx.F1[,3:9])/rowSums(CODmux.F[,3:9]) + rowSums(CODdx.M1[,3:10])/rowSums(CODmux.M[,3:10]))
Amux<-data.frame(year=year.new1,age=age.new1,female=rowSums(CODmux.F[,3:9]),
                 male=rowSums(CODmux.M[,3:10]), total=Adx[,5]/ex[,5])




CODdx.F <- CODdx.F1
head(CODdx.F)

CODdx.M <- CODdx.M1
head(CODdx.M)


#head(Adx,100)
#tail(Adx,100)
#head(Amux,100)
#tail(Amux,100)
#head(ex,100)
#tail(ex,100)




#save(Ex.new, CODDx.F1, CODDx.M1, file="extrapolated.Data")
save(Amux, Adx, ex, CODmux.F, CODmux.M, CODdx.F, CODdx.M, COD.key, file="Graduated.Data")

#chunk End

##################################################
#######   RAINBOW PLOTS ##########################
##################################################

################################
####  AGGREGATED DATA  #########
################################

mort.f.fd<-matrix(log(Amux[,3]),nrow=105)
mort.m.fd<-matrix(log(Amux[,4]),nrow=105)
mort.t.fd<-matrix(log(Amux[,5]),nrow=105)
colnames(mort.f.fd)<-c(1990:2017)
colnames(mort.m.fd)<-c(1990:2017)
colnames(mort.t.fd)<-c(1990:2017)
tyf<-fts(0:104, mort.f.fd,xname="year",yname="log.mux")
tym<-fts(0:104, mort.m.fd,xname="year",yname="log.mux")
tyt<-fts(0:104, mort.t.fd,xname="year",yname="log.mux")
par(mfrow=c(1,3))
plot.fds(tyf, plot.type = "functions",main="Female Mortality",
         cex.lab=0.7,cex.main=0.7,cex.axis=0.7)
legend("topleft",c("1990","2003","2017"),col=c("red","green","purple"),lty=rep(1,3),cex=0.6)
plot.fds(tym, plot.type = "functions",main="Male Mortality",
         cex.lab=0.7,cex.main=0.7,cex.axis=0.7)
legend("topleft",c("1990","2003","2017"),col=c("red","green","purple"),lty=rep(1,3),cex=0.6)
plot.fds(tyt, plot.type = "functions",main="Total Mortality",
         cex.lab=0.7,cex.main=0.7,cex.axis=0.7)
legend("topleft",c("1990","2003","2017"),col=c("red","green","purple"),lty=rep(1,3),cex=0.6)


################################
####  DISAGGREGATED DATA  ######
################################

#Females (on different scales)
head(CODmux.F)
mynamelist<-names(CODmux.F)
par(mfrow=c(3,3))
for(i in 3:9)
{
  CODmux.F.fd<-matrix(log(CODmux.F[,i]),nrow=105)
  colnames(CODmux.F.fd)<-c(1990:2017)
  myty<-fts(0:104, CODmux.F.fd,xname="year",yname=mynamelist[i])
  plot.fds(myty, plot.type = "functions",main=paste("Female, ", mynamelist[i]),
           cex.lab=0.7,cex.main=0.7,cex.axis=0.7,ylab="log.mux")
}
plot.fds(tyf, plot.type = "functions",main="Female, Aggregate Mortality",
         cex.lab=0.7,cex.main=0.7,cex.axis=0.7)
plot.new()
legend(0.5,1.055,c("1990","2003","2017"),col=c("red","green","purple"),lty=rep(1,3),cex=0.7)



#Females (on same scale)
mynamelist<-names(CODmux.F)
par(mfrow=c(3,3))
for(i in 3:9)
{
  CODmux.F.fd<-matrix(log(CODmux.F[,i]),nrow=105)
  colnames(CODmux.F.fd)<-c(1990:2017)
  myty<-fts(0:104, CODmux.F.fd,xname="year",yname=mynamelist[i])
  plot.fds(myty, plot.type = "functions",main=paste("Female, ", mynamelist[i]),
           cex.lab=0.7,cex.main=0.7,cex.axis=0.7,ylab="log.mux",
           ylim=c(min(log(CODmux.F[,4])),max(log(Amux[,3]))))
}
plot.fds(tyf, plot.type = "functions",main="Female, Aggregate Mortality",
         cex.lab=0.7,cex.main=0.7,cex.axis=0.7,
         ylim=c(min(log(CODmux.F[,4])),max(log(Amux[,3]))))
plot.new()
legend(0.5,1.055,c("1990","2003","2017"),col=c("red","green","purple"),lty=rep(1,3),cex=0.7)


#Males (on different scales)
head(CODmux.M)
mynamelist<-names(CODmux.M)
par(mfrow=c(3,3))
for(i in 3:10)
{
  CODmux.M.fd<-matrix(log(CODmux.M[,i]),nrow=105)
  colnames(CODmux.M.fd)<-c(1990:2017)
  myty<-fts(0:104, CODmux.M.fd,xname="year",yname=mynamelist[i])
  plot.fds(myty, plot.type = "functions",main=paste("Male, ", mynamelist[i]),
           cex.lab=0.7,cex.main=0.7,cex.axis=0.7,ylab="log.mux")
  #ylim=c(min(log(CODmux.F[,4])),max(log(CODmux.F[,3]))))
  
}
plot.fds(tym, plot.type = "functions",main="Male, Aggregate Mortality",
         cex.lab=0.7,cex.main=0.7,cex.axis=0.7)
#ylim=c(min(log(CODmux.F[,4])),max(log(Amux[,3]))))
#legend(80,-4,c("1990","2003","2017"),col=c("red","green","purple"),
#   lty=rep(1,3),cex=0.4, text.width = strwidth("000,000"))
#plot.new()
#legend(0.5,1.055,c("1990","2003","2017"),col=c("red","green","purple"),lty=rep(1,3),cex=0.7)

#Males (on same scale)
mynamelist<-names(CODmux.M)
par(mfrow=c(3,3))
for(i in 3:10)
{
  CODmux.M.fd<-matrix(log(CODmux.M[,i]),nrow=105)
  colnames(CODmux.M.fd)<-c(1990:2017)
  myty<-fts(0:104, CODmux.M.fd,xname="year",yname=mynamelist[i])
  plot.fds(myty, plot.type = "functions",main=paste("Male, ", mynamelist[i]),
           cex.lab=0.7,cex.main=0.7,cex.axis=0.7,ylab="log.mux",
           ylim=c(min(log(CODmux.M[,5])),max(log(Amux[,4]))))
}
plot.fds(tym, plot.type = "functions",main="Male, Aggregate Mortality",
         cex.lab=0.7,cex.main=0.7,cex.axis=0.7,
         ylim=c(min(log(CODmux.M[,5])),max(log(Amux[,4]))))
#plot.new()
#legend(0.5,1.055,c("1990","2003","2017"),col=c("red","green","purple"),lty=rep(1,3),cex=0.7)

#Totals (on different scales)
#head(CODmux.T)
#mynamelist<-names(CODmux.T)
#par(mfrow=c(3,3))
# for(i in 3:9)
# {
#   CODmux.T.fd<-matrix(log(CODmux.T[,i]),nrow=91)
#   colnames(CODmux.T.fd)<-c(1990:2017)
#   myty<-fts(0:90, CODmux.T.fd,xname="year",yname=mynamelist[i])
#   plot.fds(myty, plot.type = "functions",main=paste("Total, ", mynamelist[i]),
#            cex.lab=0.7,cex.main=0.7,cex.axis=0.7,ylab="log.mux")
#   #ylim=c(min(log(CODmux.F[,4])),max(log(CODmux.F[,3]))))
#   
# }
# plot.fds(tyt, plot.type = "functions",main="Total, Aggregate Mortality",
#          cex.lab=0.7,cex.main=0.7,cex.axis=0.7)
# #ylim=c(min(log(CODmux.F[,4])),max(log(Amux[,3]))))
# #legend(80,-4,c("1990","2003","2017"),col=c("red","green","purple"),
# #   lty=rep(1,3),cex=0.4, text.width = strwidth("000,000"))
# plot.new()
# legend(0.5,1.055,c("1990","2003","2017"),col=c("red","green","purple"),lty=rep(1,3),cex=0.7)
# 
# #TOtals (on same scale)
# mynamelist<-names(CODmux.T)
# par(mfrow=c(3,3))
# for(i in 3:9)
# {
#   CODmux.T.fd<-matrix(log(CODmux.T[,i]),nrow=91)
#   colnames(CODmux.T.fd)<-c(1990:2017)
#   myty<-fts(0:90, CODmux.T.fd,xname="year",yname=mynamelist[i])
#   plot.fds(myty, plot.type = "functions",main=paste("Total, ", mynamelist[i]),
#            cex.lab=0.7,cex.main=0.7,cex.axis=0.7,ylab="log.mux",
#            ylim=c(min(log(CODmux.T[,4])),max(log(Amux[,5]))))
# }
# plot.fds(tyt, plot.type = "functions",main="Male, Aggregate Mortality",
#          cex.lab=0.7,cex.main=0.7,cex.axis=0.7,
#          ylim=c(min(log(CODmux.T[,5])),max(log(Amux[,5]))))
# plot.new()
# legend(0.5,1.055,c("1990","2003","2017"),col=c("red","green","purple"),lty=rep(1,3),cex=0.7)
# 
# 
# #####################################################################
# ###########  3D VIEW OF THE MORTALITY SURFACE    ####################
# #####################################################################
# 
# #females
# plot_ly(x=seq(0,90,1),y=seq(1990,2017,1),z=Amort.f, type="surface")
# #males
# plot_ly(x=seq(0,90,1),y=seq(1990,2017,1),z=Amort.m, type="surface")
# #total
# #plot_ly(x=seq(0,90,1),y=seq(1990,2017,1),z=Amort.t, type="surface")



#THE END
