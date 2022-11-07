##############################################################################-#
## The following codes replicate our results in simulation study               #
## Author: Wei Wang                                                            #
##############################################################################-#


################################################################################
## model the simulation datasets for each scenario                             #
## this step will cost several days due to the large number of simulation times#
## Your may need a super computer to run the codes                             #
################################################################################
### load data and function -----
path <- "your path"
setwd(path)
`%+%` <- function(x,y) paste0(x,y)
source("fun-MCMAR.R") 
load("data\\simdata.Rdata")

### fit MMR and MCMAR with or without covariates -----
method <- "reml" # reml and ml are available
smvcontrol <- list(maxiter = 100,factr = 1e7,hessian = F,opt.iter = 0,opt.iter.show = F)

scennames <- names(simdata)[1:8]
for (scen in scennames ) {
  data <- simdata[[scen]]
  Sall <- data$TrueParamter$Sall
  covariate <- data$TrueParamter$covariate
  Cmatrix <- data$TrueParamter$Cmatrix
  
  nsim <- 1000
  res_mvmeta_intercept <- res_mvmeta_covariate <- vector(mode = "list",length = nsim) 
  res_smvmeta_intercept <- res_smvmeta_covariate <- vector(mode = "list",length = nsim)
  for (i in 1:nsim) {
    yall <- data$Simy[[i]]
    fit0_intercept <- mvmeta(yall,Sall,method = method)
    fit0_covariate <- mvmeta(yall ~ covariate,Sall,method = method)
    fit1_intercept <- smvmeta(yall,S=Sall,Cmatrix = Cmatrix,
                              method=method,control = smvcontrol)
    fit1_covariate <- smvmeta(yall ~ covariate,S=Sall,Cmatrix = Cmatrix,
                              method=method,control = smvcontrol)
    res_mvmeta_intercept[[i]] <- fit0_intercept
    res_mvmeta_covariate[[i]] <- fit0_covariate
    res_smvmeta_intercept[[i]] <- fit1_intercept
    res_smvmeta_covariate[[i]] <- fit1_covariate
    if(i%%100 == 0) save(res_mvmeta_intercept,res_mvmeta_covariate,res_smvmeta_intercept,
                        res_smvmeta_covariate,file = "simres_"%+%scen%+%"_"%+%method%+%".Rdata")
    cat(i," ")
  }
  cat("\n")
}




################################################################################
##### get performance for simulation study######################################
###### mv denotes MMR and smv denote MCMAR ####################################
###############################################################################

########## load function and data -----
`%+%` <- function(x,y) paste0(x,y)
`%>%` <- magrittr::`%>%`
path <- "The path for saving the fit results"
setwd(path)
load("data\\simdata.Rdata")
scennames <- names(simdata)

# get city-specific beta for MMR
Getcsbeta_mv <- function(object){
  fitvalue <- object$fitted.values
  obsy <- object$model[[1]]
  psi <- object$Psi
  SIGMA <- object$S
  res <- sapply(1:nrow(obsy), function(x){
    invs <- chol2inv(chol(mvmeta::xpndMat(SIGMA[x,]) + psi))
    fitvalue[x,] + psi %*% invs %*% (obsy[x,] - fitvalue[x,])
  })
  res <- t(res)
  rownames(res) <- rownames(obsy)
  colnames(res) <- colnames(obsy)
  res
}

GetAIC <- function(model){
  class(model) <- "mvmeta"
  df <- attr(logLik(model),"df")
  -logLik(model) * 2 + 2 * df
}
GetAIC_opt_model <- function(mvmodel,smvmodel){
  class(smvmodel) <- "mvmeta"
  logrr <- (logLik(smvmodel) - logLik(mvmodel)) * 2
  if(logrr < 3.841459) xx <- GetAIC(mvmodel)
  else xx <- GetAIC(smvmodel)
  xx
}

cummean <- function(x){
  cumsum(x)/(1:length(x))
}
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

# get error for intercept
GetErrorBeta0 <- function(model,truevalues){
  sqrt(sum((model$coefficients/truevalues - 1)^2))
}
GetErrorBeta0_opt <- function(mvmodel,smvmodel,truevalues){
  class(smvmodel) <- "mvmeta"
  logrr <- (logLik(smvmodel) - logLik(mvmodel)) * 2
  if(logrr < 3.841459) xx <- GetErrorBeta0(mvmodel,truevalues)
  else xx <- GetErrorBeta0(smvmodel,truevalues)
  xx
}

# get error for predictor
GetErrorBeta1 <- function(model,truevalues){
  if(any(truevalues == 0))
    sqrt(sum((model$coefficients[2,]-truevalues)^2))
  else 
    sqrt(sum((model$coefficients[2,]/truevalues - 1)^2))
}
GetErrorBeta1_opt <- function(mvmodel,smvmodel,truevalues){
  class(smvmodel) <- "mvmeta"
  logrr <- (logLik(smvmodel) - logLik(mvmodel)) * 2
  if(logrr < 3.841459) xx <- GetErrorBeta1(mvmodel,truevalues)
  else xx <- GetErrorBeta1(smvmodel,truevalues)
  xx
}

# get error for predictive ERR
GetErrorERR <- function(model,truevalues,method = "MMR"){
  if(method == "MMR") 
    predvalue <- Getcsbeta_mv(model)
  else predvalue <- model$fitted.values.spatial
  nn <- nrow(predvalue)
  aa <- sapply(1:nn, function(x) sqrt(sum((predvalue[x,]-truevalues[x,])^2)))
  mean(aa)
}
GetErrorERR_opt <- function(mvmodel,smvmodel,truevalues){
  class(smvmodel) <- "mvmeta"
  logrr <- (logLik(smvmodel) - logLik(mvmodel)) * 2
  if(logrr < 3.841459) xx <- GetErrorERR(mvmodel,truevalues,"MMR")
  else xx <- GetErrorERR(smvmodel,truevalues,"MCMAR")
  xx
}

########## get performance for specific scenario ------
AMAE_pbeta <- AMAE_ERR <- aic <- power <- NULL
method <- "reml"
for (scen in scennames) {
  load("simres_"%+%scen%+%"_"%+%method%+%".Rdata")
  truedata <- simdata[[scen]]
  m <- 143
  k <- 5
  ## pooled beta -----------
  pbeta_true <- truedata$TrueParamter$beta
  pbeta_true[1,] <- mean(truedata$TrueParamter$covariate) * pbeta_true[2,] + pbeta_true[1,] # the true average association across all regions
  pbeta_MAE_mv_intercept <- sapply(res_mvmeta_intercept, GetErrorBeta0,
                                   truevalues = pbeta_true[1,]) %>% mean %>% round(3) 
  pbeta_MAE_smv_intercept <- sapply(res_smvmeta_intercept, GetErrorBeta0,
                                    truevalues = pbeta_true[1,]) %>% mean %>% round(3)
  pbeta_MAE_opt_intercept <- sapply(1:1000,function(x) GetErrorBeta0_opt(res_mvmeta_intercept[[x]],
                                                                         res_smvmeta_intercept[[x]],pbeta_true[1,])) %>% mean %>% round(3)
  pbeta_MAE_mv_covariate <- sapply(res_mvmeta_covariate, GetErrorBeta1,
                                   truevalues = pbeta_true[2,]) %>% mean %>% round(3) 
  pbeta_MAE_smv_covariate <- sapply(res_smvmeta_covariate, GetErrorBeta1,
                                    truevalues = pbeta_true[2,]) %>% mean %>% round(3)
  pbeta_MAE_opt_covariate <- sapply(1:1000,function(x) GetErrorBeta1_opt(res_mvmeta_covariate[[x]],
                                                                         res_smvmeta_covariate[[x]],pbeta_true[2,])) %>% mean %>% round(3)
  pbeta_MAE <- c(pbeta_MAE_mv_intercept,pbeta_MAE_smv_intercept,pbeta_MAE_opt_intercept,
                 pbeta_MAE_mv_covariate,pbeta_MAE_smv_covariate,pbeta_MAE_opt_covariate)
  AMAE_pbeta <- rbind(AMAE_pbeta,pbeta_MAE)
  
  ## predictive ERR
  ERR_true <- truedata$TrueParamter$truespatialfitvalue
  ERR_MAE_mv_intercept <- sapply(res_mvmeta_intercept, GetErrorERR,method = "MMR",
                                 truevalues = ERR_true) %>% mean %>% round(3) 
  ERR_MAE_smv_intercept <- sapply(res_smvmeta_intercept, GetErrorERR,method = "MCMAR",
                                  truevalues = ERR_true) %>% mean %>% round(3) 
  ERR_MAE_opt_intercept <- sapply(1:1000,function(x) GetErrorERR_opt(res_mvmeta_intercept[[x]],
                                                                     res_smvmeta_intercept[[x]],ERR_true)) %>% mean %>% round(3)
  ERR_MAE_mv_covariate <- sapply(res_mvmeta_covariate, GetErrorERR,method = "MMR",
                                 truevalues = ERR_true) %>% mean %>% round(3) 
  ERR_MAE_smv_covariate <- sapply(res_smvmeta_covariate, GetErrorERR,method = "MCMAR",
                                  truevalues = ERR_true) %>% mean %>% round(3) 
  ERR_MAE_opt_covariate <- sapply(1:1000,function(x) GetErrorERR_opt(res_mvmeta_covariate[[x]],
                                                                     res_smvmeta_covariate[[x]],ERR_true)) %>% mean %>% round(3)
  ERR_MAE <- c(ERR_MAE_mv_intercept,ERR_MAE_smv_intercept,ERR_MAE_opt_intercept,
               ERR_MAE_mv_covariate,ERR_MAE_smv_covariate,ERR_MAE_opt_covariate)
  AMAE_ERR <- rbind(AMAE_ERR,ERR_MAE)
  
  ## AIC----------
  aic_mv_intercept <- sapply(res_mvmeta_intercept, GetAIC)
  aic_mv_covariate <- sapply(res_mvmeta_covariate, GetAIC)
  aic_smv_intercept <- sapply(res_smvmeta_intercept, GetAIC)
  aic_smv_covariate <- sapply(res_smvmeta_covariate, GetAIC)
  aic_opt_intercept <- sapply(1:1000, function(x) GetAIC_opt_model(res_mvmeta_intercept[[x]],res_smvmeta_intercept[[x]]))
  aic_opt_covariate <- sapply(1:1000, function(x) GetAIC_opt_model(res_mvmeta_covariate[[x]],res_smvmeta_covariate[[x]]))
  
  aici <- c(mean(aic_mv_intercept),mean(aic_smv_intercept),mean(aic_opt_intercept),
            mean(aic_mv_covariate),mean(aic_smv_covariate),mean(aic_opt_covariate))
  aic <- rbind(aic,aici) %>% round(3)
  
  ## wald test for the covariate -----
  nobs <- 1000
  getp <- function(mvmodel,smvmodel){
    class(smvmodel) <- "mvmeta"
    logrr <- (logLik(smvmodel) - logLik(mvmodel)) * 2
    if(logrr < 3.841459) xx <- waldtest(mvmodel)[3]
    else xx <- waldtest(smvmodel)[3]
    xx
  }
  p_mv <- sapply(res_mvmeta_covariate, function(x) waldtest(x)[3])
  p_smv <- sapply(res_smvmeta_covariate, function(x) waldtest(x)[3])
  p_opt <- sapply(1:1000, function(x) getp(res_mvmeta_covariate[[x]],res_smvmeta_covariate[[x]]))
  poweri <- c(sum(p_mv<0.05),sum(p_smv<0.05),sum(p_opt<0.05))
  power <- rbind(power,poweri) %>% round(3)
  
  cat(scen," ")
  rm(res_smvmeta_covariate,res_smvmeta_intercept,res_mvmeta_covariate,res_mvmeta_intercept)
  gc()
}

colnames(AMAE_pbeta) <- c("mv_intercept","smv_intercept","opt_intercept",
                          "mv_covariate","smv_covariate","opt_covariate")
colnames(AMAE_ERR) <- c("mv_intercept","smv_intercept","opt_intercept",
                        "mv_covariate","smv_covariate","opt_covariate")
colnames(aic) <- c("mv_intercept","smv_intercept","opt_intercept",
                   "mv_covariate","smv_covariate","opt_covariate")
colnames(power) <- c("mv_covariate","smv_covariate","opt_covariate")
rownames(AMAE_pbeta) <- rownames(AMAE_ERR) <- rownames(aic) <- rownames(power) <- scennames
xlsx::write.xlsx(AMAE_pbeta,file = "performance-results-"%+%method%+%".xlsx",sheetName = "AMAE_pbeta")
xlsx::write.xlsx(AMAE_ERR,file = "performance-results-"%+%method%+%".xlsx",sheetName = "AMAE_ERR",append = T)
xlsx::write.xlsx(aic,file = "performance-results-"%+%method%+%".xlsx",sheetName = "aic",append = T)
xlsx::write.xlsx(power,file = "performance-results-"%+%method%+%".xlsx",sheetName = "power",append = T)


########## get 95% coverage rates for each scenario============================
Getcoverage <- function(model,mu){
  a <- nrow(model$coefficients)
  if(a == 1) ind <- 1:5 else if(a==2) ind <- seq(2,10,2)
  s2 <- model$vcov[ind,ind]
  est <- model$coefficients[a,] - mu[a,]
  b <- est %*% solve(s2) %*% est
  pchisq(q = b,df = length(est)) < 0.95
}
Getcoverage_opt <- function(mvmodel,smvmodel,mu){
  class(smvmodel) <- "mvmeta"
  logrr <- (logLik(smvmodel) - logLik(mvmodel)) * 2
  if(logrr < 3.841459) xx <- Getcoverage(mvmodel,mu)
  else xx <- Getcoverage(smvmodel,mu)
  xx
}


coverage95 <- NULL
method <- method
for (scen in scennames) {
  # scen <- scennames[4]
  load(path%+%"results/simres_"%+%scen%+%"_"%+%method%+%".Rdata")
  truedata <- simdata[[scen]]
  truebeta <- truedata$TrueParamter$beta
  truebeta[1,] <- mean(truedata$TrueParamter$covariate) * truebeta[2,] + truebeta[1,] # the true average association across all regions
  
  coverate_mvmeta_intercept <- sapply(res_mvmeta_intercept, Getcoverage, mu = truebeta) %>% mean
  coverate_smvmeta_intercept <- sapply(res_smvmeta_intercept, Getcoverage, mu = truebeta) %>% mean
  coverate_opt_intercept <- sapply(1:1000, function(x) Getcoverage_opt(res_mvmeta_intercept[[x]],res_smvmeta_intercept[[x]],mu = truebeta)) %>% mean  
  
  coverate_mvmeta_covariate <- sapply(res_mvmeta_covariate, Getcoverage, mu = truebeta) %>% mean
  coverate_smvmeta_covariate <- sapply(res_smvmeta_covariate, Getcoverage, mu = truebeta) %>% mean
  coverate_opt_covariate <- sapply(1:1000, function(x) Getcoverage_opt(res_mvmeta_covariate[[x]],res_smvmeta_covariate[[x]],mu = truebeta)) %>% mean  
  
  a <- c(coverate_mvmeta_intercept,coverate_smvmeta_intercept,coverate_opt_intercept,
         coverate_mvmeta_covariate,coverate_smvmeta_covariate,coverate_opt_covariate)
  coverage95 <- rbind(coverage95,a)
  colnames(coverage95) <- c("mv_intercept","smv_intercept","opt_intercept",
                            "mv_covariate","smv_covariate","opt_covariate")
  
  rm(res_smvmeta_covariate,res_smvmeta_intercept,res_mvmeta_covariate,res_mvmeta_intercept)
  gc()
  cat(scen," ")
}
rownames(coverage95) <- scennames
xlsx::write.xlsx(coverage95,file = "performance-results-"%+%method%+%".xlsx",sheetName = "coverage95",append = T)











