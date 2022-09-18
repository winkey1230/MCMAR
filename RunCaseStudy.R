##############################################################################-#
## The following codes replicate our results in the motivating example         #
## mv denotes MMR and smv denotes MCMAR                                        #
## Author: Wei Wang                                                            #
##############################################################################-#


################################################################################
##### load data,library and function ###########################################
###############################################################################
`%+%` <- function(x,y) paste0(x,y)
`%>%` <- magrittr::`%>%` # pipe operator
library(sp) 
library(spdep) 
library(dlnm)
library(tmap)
library(RColorBrewer)
path <- "your path"
source(path%+%"fun-MCMAR.R")
load(path%+%"data\\stage1data.Rdata")
load(path%+%"data\\mapdata_China_province_project.Rdata")
# wald test for MMR and MCMAR
waldtest <- function(model){
  m <- model
  coef <- m$coefficients[2,]
  vcov <- m$vcov[seq(2,10,2),seq(2,10,2)]
  waldstat <- coef%*%solve(vcov)%*%coef
  df <- length(coef)
  pvalue <- 1-pchisq(waldstat,df)
  wald <- c(waldstat,df,pvalue)
  wald
}
voronoipolygons = function(layer) {
  # layer: sp class
  require(deldir)
  require(sp)
  crds = layer@coords
  z = deldir(crds[,1], crds[,2])
  w = tile.list(z)
  polys = vector(mode = 'list', length = length(w))
  for (i in seq(along = polys)) {
    pcrds = cbind(w[[i]]$x, w[[i]]$y)
    pcrds = rbind(pcrds, pcrds[1,])
    polys[[i]] = Polygons(list(Polygon(pcrds)),
                          ID = as.character(i))
  }
  SP = SpatialPolygons(polys)
  voronoi = SpatialPolygonsDataFrame(SP, data = data.frame(x = crds[,1], 
                                                           y = crds[,2], row.names = sapply(slot(SP, 'polygons'), 
                                                                                            function(x) slot(x, 'ID'))))
}
################################################################################
##### Get the spatially adjacent matrix            #############################
################################################################################
## Thiessen polygon-based method-----
Rmethod <- "_TS"
temp <- CINF
coordinates(temp) <- ~ POINT_X + POINT_Y 
tspolygon <- voronoipolygons(temp)
plot(tspolygon)
points(CINF$POINT_X,CINF$POINT_Y,pch = 16)
nb <- poly2nb(tspolygon)
W <- nb2mat(nb)
n <- ncol(W)
C <- W * apply(W, 1, function(x) sum(x!=0))
R <- diag(rowSums(C)) - C
Cmatrix <- diag(n) - R
## k-nearest neighbor method-----
k <- 4
Rmethod <- "_"%+%k%+%"n"
if(k == 4) Rmethod <- ""
temp <- CINF
coordinates(temp) <- ~ POINT_X + POINT_Y 
nb <- knearneigh(temp, k=k) 
nb <- knn2nb(nb)
plot(nb,coordinates(temp),pch = 16)
W <- nb2mat(nb)
n <- ncol(W)
C <- W * apply(W, 1, function(x) sum(x!=0))
# get symatric C
for (i in 1:(n-1)) {
  for (j in (i+1):n) {
     if(C[i,j]!= C[j,i]){
       C[i,j] <- 1
       C[j,i] <- 1
     } 
  }
}
R <- diag(rowSums(C)) - C
Cmatrix <- diag(n) - R

##############################################################################
############# Compare MMR with MCMAR with only intercept#################
##############################################################################
method <- "ml"
xvar <- seq(0,100,by=5)
bvar <- do.call("onebasis",c(list(x=xvar),attr(M.CB,"argvar")))
dfvar <- 5

## MMR
fit0 <- mvmeta(yall,Sall,method = method)
mvmetafit_intercept <- fit0
aic0 <- (fit0$logLik * (-2) + 2 * dfvar + dfvar*(dfvar+1)) %>% round(1)

# MCMAR
system.time(
  # This will cost five minutes. opt.iter can be set as 1 to save time
  fit1 <- smvmeta(yall,S=Sall,Cmatrix = Cmatrix,method=method,
                  control = list(maxiter = 200,factr = 1e7,hessian = T,
                                 opt.iter = 10,opt.iter.show = T))
)
# save(fit1,file = path%+%"data\\fit1_smvmeta_onlyintercept_"%+%method%+%Rmethod%+%".Rdata")
# load(path%+%"data\\fit1_smvmeta_onlyintercept_"%+%method%+%Rmethod%+%".Rdata")
smvmetafit_intercept <- fit1
aic1 <- (fit1$logLik * (-2) + 2 * dfvar + dfvar*(dfvar+1) + 2)  %>% round(1)
rho1 <- fit1$rho; rho1_se <- fit1$se.rho

# plot the average rr curve
tiff(filename=path%+%"Average_ERR.tiff", width=15.9, height=8, units="cm",pointsize = 8,res=300)
cpall0 <- crosspred(bvar, coef= as.numeric(fit0$coefficients), vcov=fit0$vcov, model.link="log",by=5)
x<-seq(0,100,by=5)
par(cex.axis=1,cex.lab=1,cex.main=1.2, mar=c(4,4,3.6,0.8))
plot(cpall0,ylim=c(0.6,1.8),type="n", ci="n",xlab="Temperature on relative scale (%)",
     ylab="Risk Ratio (RR)", main=NULL)
redtrans <- rgb(255, 0, 0, 60, maxColorValue=255) 
points(xvar,cpall0$allRRfit, type="l", lwd=1, col="red")
polygon(c(rev(x),x),c(rev(cpall0$allRRhigh),cpall0$allRRlow),col=redtrans, border = NA)

bluetrans <- rgb(0, 0,255, 120, maxColorValue=255) 
cpall1 <- crosspred(bvar, coef= as.numeric(fit1$coefficients), vcov=fit1$vcov, model.link="log",by=5)
lines(cpall1,col="blue",lty=4,lwd=1,ci="area",ci.arg=list(density=20,col=bluetrans))
legend(x="top",inset =0, legend=c("MMR      (AIC = "%+%aic0%+%")", "MCMAR (AIC = "%+%aic1%+%")"),
       lwd=1.5, lty=1, col=c("red", "blue"), bty="n",ncol=1, cex=1.6)
dev.off()

# plot the map for RR refferring to 50
objtem <- c(10,30,70,90)
predbeta0 <- predmvmeta(fit0)
predbeta1 <- fit1$fitted.values.spatial
objbasis <- do.call("onebasis",c(list(x=objtem),attr(M.CB,"argvar"))) 
refbasis <- do.call("onebasis",c(list(x=50),attr(M.CB,"argvar")))
difx <- t(t(objbasis) - as.numeric(refbasis))
allRRfit0 <- exp(predbeta0 %*% t(difx)) 
allRRfit1 <- exp(predbeta1 %*% t(difx))
!any(rownames(allRRfit0) != CINF$City) # check the order
allRRfit <- cbind(CINF$citycd,allRRfit0,allRRfit1) %>% as.data.frame()
names(allRRfit) <- c("citycd","fit0_"%+%objtem,"fit1_"%+%objtem)

citypoints <- CINF[,1:6]
citypoints <- merge(citypoints,allRRfit,by = "citycd")
sp::coordinates(citypoints) <- c("POINT_X","POINT_Y")
citypoints@proj4string <- mapdata@proj4string
tempmap <- sf::st_as_sf(mapdata) # transform sp to sf class
citypoints <- sf::st_as_sf(citypoints)

colornames <- rep(c("fit0_","fit1_"),length(objtem)) %+% rep(objtem,each = 2)
panelnames <- rep(c("MMR (","MCMAR ("),length(objtem)) %+% rep(objtem,each = 2) %+% "%)"
x <- as.data.frame(citypoints)[,colornames] %>% unlist()
a <- round(range(x) + c(-0.005,0.005),2)
lowcolor <- c(a[1],round(quantile(x[x<1],c(0.05,0.25,0.5,0.75)),2),1)
highcolor <- c(round(quantile(x[x>1],c(0.25,0.5,0.75,0.95)),2),a[2]) 
b <- c(lowcolor,highcolor) # c(x[1],0.75,0.85,0.9,1.0,1.1,1.2,1.5,x[2])
colorlables <- vector(mode = "character",length = length(b)-1)
for (i in 1:(length(b)-1)) {
  colorlables[i] <- "["%+%b[i]%+%", "%+%b[i+1]%+%")"
}
paleta <- brewer.pal(length(colorlables),"RdYlGn")[length(colorlables):1]
pdf(path%+%"RR_Map_for_Casestudy_onlyintercept.pdf")
tmap_mode("plot")
tm_shape(tempmap,ylim = c(1500000,5900000)) +  # limiting the range of coordinate
  tm_polygons(col = "grey100") + 
  tm_shape(citypoints) + 
  tm_bubbles(size = 0.15,col = colornames,palette=paleta,
             style="fixed", breaks=b,interval.closure="left",labels=colorlables) +
  tm_layout(panel.labels=panelnames,legend.outside=T,legend.title.size = 0.8,
            legend.outside.size = 0.2,title = "RR",legend.outside.position = "right") +
  tm_facets(ncol = 2,nrow = 4) 
dev.off()

# tmap_mode("view")
# tm_shape(tempmap) + 
#   tm_polygons(col = "grey100")

###############################################################################
############ Compare MMR with MCMAR with a covariate #########################
##############################################################################
method <- "ml"
xvar <- seq(0,100,by=5)
bvar <- do.call("onebasis",c(list(x=xvar),attr(M.CB,"argvar")))
dfvar <- 5
CINF$`GDP per person` <- CINF$GDP/CINF$Population
allcovariate <- c("Latitude","Longitude","Altitude","Temperature","Relative humidity",
                  "Air pressure","Rainfall","Sunshine","Population increase",
                  "Population density","GDP per person","GDP increase",
                  "Licensed physicians per 1000 persons","Hospital beds per 1000 persons",
                  "Traffic","Students per 1000 persons")
fit0 <- mvmeta(yall,Sall,method = method);mvmetafit_intercept <- fit0
fit1 <- smvmeta(yall,S=Sall,Cmatrix = Cmatrix,method=method,
                control = list(maxiter = 200,factr = 1e7,hessian = T,
                               opt.iter = 10,opt.iter.show = T))
aic0 <- fit0$logLik * (-2) + 2 * dfvar + dfvar*(dfvar+1)
aic1 <- fit1$logLik * (-2) + 2 * dfvar + dfvar*(dfvar+1) + 2
fitres <- c(round(aic0,1),NA,round(GetH(fit0)[4],1),round(aic1,1),
            NA,round(rhotest(fit0,fit1)[-2],3))
smvmetafit_intercept <- fit1
for (i in allcovariate) {
  covariate <- CINF[,i]
  atcovar <- quantile(covariate,0.95,na.rm = T)
  aa <- t(diag(dfvar) %x% c(1,0))
  bb <- t(diag(dfvar) %x% c(0,1))
  cc <- aa + atcovar * bb
  ## MMR
  fit0 <- mvmeta(yall ~ covariate,Sall,method = method)
  mvmetafit_covar <- fit0
  aic0 <- fit0$logLik * (-2) + 2 * dfvar * 2 + dfvar*(dfvar+1)
  mvmetaresi <- c(round(aic0,1),round(waldtest(mvmetafit_covar)[3],3),round(GetH(mvmetafit_covar)[4],1))
  coefs0 <- as.vector(fit0$coefficients)
  beta0 <- as.numeric(cc%*%coefs0) 
  vcov0 <- cc %*% fit0$vcov %*% t(cc)
  # MCMAR
  system.time(
    fit1 <- smvmeta(yall ~ covariate,S=Sall,Cmatrix = Cmatrix,method=method,
                    control = list(maxiter = 200,factr = 1e7,hessian = T,
                                   opt.iter = 10,opt.iter.show = T))
  )
  cat("\n")
  # save(fit1,file = path%+%"data\\fit1_smvmeta_covariate_"%+%i%+%"_"%+%method%+%Rmethod%+%".Rdata")
  # load(path%+%"data\\fit1_smvmeta_covariate_"%+%i%+%"_"%+%method%+%Rmethod%+%".Rdata")
  smvmetafit_covar <- fit1
  aic1 <- fit1$logLik * (-2) + 2 * dfvar * 2 + dfvar*(dfvar+1) + 2
  smvmetaresi <- c(round(aic1,1),round(waldtest(smvmetafit_covar)[3],3),round(rhotest(fit0,fit1)[-2],3))
  rho1 <- fit1$rho; rho1_se <- fit1$se.rho
  coefs1 <- as.vector(fit1$coefficients)
  beta1 <- as.numeric(cc%*%coefs1) 
  vcov1 <- cc %*% fit1$vcov %*% t(cc)
  fitres <- rbind(fitres,c(mvmetaresi,smvmetaresi))
}
rownames(fitres) <- c("Intercept only",allcovariate)
colnames(fitres) <- c("mvmeta_aic","mvmeta_wald_p","I^2","smvmeta_aic",
                      "smvmeta_wald_p","rho","rho_p")
xlsx::write.xlsx(fitres,file = path%+%"case_study_covariate_"%+%method%+%Rmethod%+%".xlsx")




