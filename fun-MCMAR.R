###############################################################################-#
# The following codes include functions implementing MCMAR using two parameter #
# estimation methods: ml and reml                                              #
##############################################################################-#


#'@description smvmeta implement the parameter estimation of MCMAR
#'@author Wang
#'@param formula An object of class "formula", see mvmeta::mvmeta.
#'@param S The covariance for estimated outcomes in a stratified analysis. See \r
#' in mvmeta::mvmeta.
#' @param data,subset,offset,na.action,model,offset,na.action seen mvmeta::mvmeta.
#' @param method Currently, only support "ml" and "reml". 
#' @param Cmatrix A matrix denoting spatial autocorelation pattern. Equal to I-R in MCMAR.
#' @param bscov Define the stucture of V in MCMAR. Currently,only support "unstr".
#' @param control seen the following function smvmeta.control
smvmeta <- function (formula, S, data, subset,Cmatrix = NULL, method = "reml", bscov = "unstr", 
                     model = TRUE, contrasts = NULL, offset, na.action, control = list()){
  call <- match.call()
  mcall <- match.call(expand.dots = FALSE)
  mn <- match(c("formula", "data", "subset", "weights", "na.action", 
                "offset"), names(mcall), 0L)
  mcall <- mcall[c(1L, mn)]
  mcall$drop.unused.levels <- TRUE
  mcall[[1L]] <- as.name("model.frame")
  if (missing(data)) 
    data <- parent.frame()
  if (!inherits(eval(substitute(formula), data), "formula")) {
    formula <- as.formula(paste(deparse(substitute(formula), 
                                        width.cutoff = 499L), "~ 1"), env = parent.frame())
    environment(formula) <- parent.frame()
    call[[mn[1]]] <- mcall[[mn[1]]] <- formula
  }
  if (missing(data)) 
    data <- environment(formula)
  mcall$na.action <- "na.pass"
  mf <- eval(mcall, parent.frame())
  class(mf) <- c("data.frame.mvmeta", class(mf))
  if (missing(na.action)) 
    na.action <- getOption("na.action")
  if (length(na.action)) 
    mf <- do.call(na.action, list(mf))
  if (method == "model.frame") 
    return(mf)
  if (is.empty.model(mf)) 
    stop("empty model not allowed")
  method <- match.arg(method, c("fixed", "ml", "reml", "mm", 
                                "vc"))
  bscov <- match.arg(bscov, c("unstr", "diag", "id", "cs", 
                              "hcs", "ar1", "prop", "cor", "fixed"))
  if (bscov != "unstr" && !method %in% c("ml", "reml")) 
    stop("structured Psi only available for methods 'ml' or 'reml'")
  terms <- attr(mf, "terms")
  y <- as.matrix(model.response(mf, "numeric"))
  X <- model.matrix(terms, mf, contrasts)
  offset <- as.vector(model.offset(mf))
  if (!is.null(offset)) {
    if (length(offset) != NROW(y)) 
      stop("number of offsets should equal number of observations")
  }
  S <- eval(call$S, data, parent.frame())
  S <- mvmeta:::mkS(S, y, attr(mf, "na.action"), if (missing(subset)) 
    NULL
    else eval(call$subset, data, parent.frame()))
  if (nrow(y) < 2L) 
    stop("less than 2 valid studies after exclusion of missing")
  
  if(is.null(control$initial_V)){
    mvcontrol <- mvmeta:::mvmeta.control()
    mvcontrol <- modifyList(mvcontrol,control)
    fit0 <- mvmeta::mvmeta.fit(X, y, S, offset, method, bscov, list())
    control$initial_V <- fit0$Psi
  }
  
  fit <- smvmeta.fit(X, y, S, Cmatrix,offset, method, bscov, control)
  fit$model <- if (model) mf  else NULL
  fit$S <- S
  fit$na.action <- attr(mf, "na.action")
  fit$call <- call
  fit$formula <- formula
  fit$terms <- terms
  fit$contrasts <- attr(X, "contrasts")
  fit$xlevels <- .getXlevels(terms, mf)
  class(fit) <- "smvmeta"
  fit
}




library(mvmeta)
mlprof.fn <- function (par,Xstar, Ystar, Cmatrix, D, k, m){
  U <- chol2inv(chol(diag(m) - par[1] * Cmatrix))
  V <- par2V(par[-1],k)
  gls <- glsfit_s(Xstar, Ystar, D, V, U, onlycoef = FALSE)
  pconst <- -0.5 * m * k * log(2 * pi)
  pres <- -0.5 * (crossprod(gls$invtMY - gls$invtMX %*% gls$coef))
  pdet <- -sum(log(diag(gls$M)))
  as.numeric(pconst + pdet + pres)
}
mlprof.fnU <- function (par,Xstar, Ystar, Cmatrix, V, D, k, m){
  U <- chol2inv(chol(diag(m) - par * Cmatrix))
  gls <- glsfit_s(Xstar, Ystar, D, V, U, onlycoef = FALSE)
  pconst <- -0.5 * m * k * log(2 * pi)
  pres <- -0.5 * (crossprod(gls$invtMY - gls$invtMX %*% gls$coef))
  pdet <- -sum(log(diag(gls$M)))
  as.numeric(pconst + pdet + pres)
}
mlprof.fnV <- function (par,Xstar, Ystar, U, D, k, m){
  V <- par2V(par,k)
  gls <- glsfit_s(Xstar, Ystar, D, V, U, onlycoef = FALSE)
  pconst <- -0.5 * m * k * log(2 * pi)
  pres <- -0.5 * (crossprod(gls$invtMY - gls$invtMX %*% gls$coef))
  pdet <- -sum(log(diag(gls$M)))
  as.numeric(pconst + pdet + pres)
}
remlprof.fn <- function (par,Xstar, Ystar, Cmatrix, D, k, m){
  U <- chol2inv(chol(diag(m) - par[1] * Cmatrix))
  V <- par2V(par[-1],k)
  gls <- glsfit_s(Xstar, Ystar, D, V, U, onlycoef = FALSE)
  pconst <- -0.5 * (m * k - length(gls$coef))* log(2 * pi)
  pres <- -0.5 * (crossprod(gls$invtMY - gls$invtMX %*% gls$coef))
  pdet1 <- -sum(log(diag(gls$M)))
  tXWXtot <- crossprod(gls$invtMX)
  pdet2 <- -sum(log(diag(chol(tXWXtot))))
  as.numeric(pconst + pdet1 + pdet2 + pres)
}
remlprof.fnU <- function (par,Xstar, Ystar, Cmatrix, V, D, k, m){
  U <- chol2inv(chol(diag(m) - par * Cmatrix))
  gls <- glsfit_s(Xstar, Ystar, D, V, U, onlycoef = FALSE)
  pconst <- -0.5 * (m * k - length(gls$coef))* log(2 * pi)
  pres <- -0.5 * (crossprod(gls$invtMY - gls$invtMX %*% gls$coef))
  pdet1 <- -sum(log(diag(gls$M)))
  tXWXtot <- crossprod(gls$invtMX)
  pdet2 <- -sum(log(diag(chol(tXWXtot))))
  as.numeric(pconst + pdet1 + pdet2 + pres)
}
remlprof.fnV <- function (par,Xstar, Ystar, U, D, k, m){
  V <- par2V(par,k)
  gls <- glsfit_s(Xstar, Ystar, D, V, U, onlycoef = FALSE)
  pconst <- -0.5 * (m * k - length(gls$coef))* log(2 * pi)
  pres <- -0.5 * (crossprod(gls$invtMY - gls$invtMX %*% gls$coef))
  pdet1 <- -sum(log(diag(gls$M)))
  tXWXtot <- crossprod(gls$invtMX)
  pdet2 <- -sum(log(diag(chol(tXWXtot))))
  as.numeric(pconst + pdet1 + pdet2 + pres)
}
remlprof.grV <- function (par, Xstar, Ystar, U, D, k, m){
  L <- diag(0, k)
  L[lower.tri(L, diag = TRUE)] <- par
  G <- t(L)
  V <- crossprod(G)
  gls <- glsfit_s(Xstar, Ystar, D, V, U, onlycoef = FALSE)
  tXWXtot <- crossprod(gls$invtMX) 
  invtXWXtot <- chol2inv(chol(tXWXtot))
  invSIGMA <- tcrossprod(gls$invM)
  res <- Ystar - Xstar %*% gls$coef 
  ind1 <- rep(1:k, k:1)
  ind2 <- unlist(sapply(1:k, seq, to = k))
  grad <- sapply(seq(length(par)), function(i) {
    A <- B <- C <- diag(0, k)
    A[ind2[i], ] <- B[, ind2[i]] <- G[ind1[i], ]
    C[ind2[i], ] <- C[, ind2[i]] <- 1
    DD <- C * A + C * B
    D <- U %x% DD
    E <- invSIGMA %*% D %*% invSIGMA
    Fe <- crossprod(res, E) %*% res
    G <- sum(diag(invSIGMA %*% D))
    H <- sum(diag(invtXWXtot %*% crossprod(Xstar, E) %*% Xstar))
    as.numeric(0.5 * (Fe - G + H))
  })
  grad
}
opt.iter.UV <- function(fnU,fnV,V,opt.iter,opt.iter.method,Xstar,Ystar,
                        D,k,m,Cmatrix,lower,upper,opt.iter.show){
  optU <- optimize(f = fnU,lower = lower,upper = upper,maximum = T,Xstar = Xstar, 
                   Ystar = Ystar, D = D, k = k, m = m, V = V, Cmatrix = Cmatrix)
  if(opt.iter.show) cat(optU$objective," ")
  llvalues <- optU$objective
  rho <- optU$maximum
  U <- chol2inv(chol(diag(m)-rho*Cmatrix))
  V <-as.matrix(Matrix::nearPD(V)[[1]])  
  parV <- vechMat(t(chol(V)))
  if(opt.iter < 2) par <- c(rho,parV) else{
    for (i in 1:(opt.iter-1)) {
      # the gr function seems not to elevate the computaion speed, so gr = NULL
      optV <- optim(par = parV, fn = fnV, gr = NULL, Xstar = Xstar, 
                    Ystar = Ystar, D = D, k = k, m = m, U = U,
                    method = opt.iter.method, control = list(fnscale = -1), hessian = F)
      V <- par2V(optV$par,k)
      optU <- optimize(f = fnU,lower = lower,upper = upper,maximum = T,Xstar = Xstar, 
                       Ystar = Ystar, D = D, k = k, m = m, V = V, Cmatrix = Cmatrix)
      rho <- optU$maximum
      U <- chol2inv(chol(diag(m)-rho*Cmatrix))
      parV <- optV$par
      if(opt.iter.show) cat(optV$value,optU$objective," ")
      llvalues <- c(llvalues,optV$value,optU$objective)
    }
    par <- c(rho,parV)
  }
  list(par=par,llvalues=llvalues)
}

par2V <- function(par,k){
  L <- diag(0, k)
  L[lower.tri(L, diag = TRUE)] <- par
  tcrossprod(L)
}

#'@description Define the optimal parameters in smvmeta
#'@param opt.iter.method The optimal iterative method for estimation V
#'@param factr Controls the convergence of the "L-BFGS-B" method.Only for the \r
#'final optimal function in MCMAR, i.e., the function with respect to beta, U, \r
#'and V at the same time.
#'@param hessian Logical. Should a numerically differentiated Hessian matrix be returned?
#'@param opt.iter A integer larger 0. The iterative time for V-U.
#'@param opt.iter.show Logical. show whether the iterative values of AIC are printed.
#'@param lltol judge if the log-likelihood is convergency in the iterative process of V-U
#'@param others seen mvmeta::mvmeta.control 
smvmeta.control <- function (optim = list(), initial_V = NULL, showiter = FALSE, maxiter = 100, initPsi = NULL, 
                             Psifix = NULL, Psicor = 0, Scor = 0, inputna = FALSE, inputvar = 10^4, 
                             hessian = FALSE, vc.adj = TRUE, factr = 1e7, lltol = 1e-5,
                             set.negeigen = sqrt(.Machine$double.eps),opt.iter.method = "BFGS",
                             opt.iter = 1,opt.iter.show = F) {
  optim <- modifyList(list(fnscale = -1, maxit = maxiter, factr = factr), 
                      optim)
  if (showiter) {
    optim$trace <- 6
    optim$REPORT <- 1
  }
  list(optim = optim, initial_V = initial_V, showiter = showiter, maxiter = maxiter, 
       hessian = hessian, initPsi = initPsi, Psifix = Psifix, 
       Psicor = Psicor, Scor = Scor, inputna = inputna, inputvar = inputvar, 
       vc.adj = vc.adj, factr = factr, set.negeigen = set.negeigen,lltol= lltol,
       opt.iter.method = opt.iter.method,opt.iter = opt.iter,opt.iter.show = opt.iter.show)
}

glsfit_s <- function (Xstar, Ystar, D, V, U, onlycoef = TRUE){
  SIGMA <- D + U %x% V
  M <- chol(SIGMA)
  invM <- backsolve(M,diag(ncol(M)))
  invtMX <- crossprod(invM,Xstar)
  invtMY <- crossprod(invM,Ystar)
  coef <- as.numeric(qr.solve(invtMX, invtMY))
  if (onlycoef) 
    return(coef)
  list(coef = coef, SIGMA = SIGMA, M = M, invM = invM,invtMX = invtMX, invtMY = invtMY)
}
smvmeta.ml <- function (Xstar, Ystar, D = D, Cmatrix = Cmatrix, k, m, p, bscov, 
                        nall,control, ...){
  fn <- mlprof.fn
  fnU <- mlprof.fnU
  fnV <- mlprof.fnV
  V <- control$initial_V#as.matrix(Matrix::nearPD(xpndMat(optV$par))$mat) 
  lambda <- 1/eigen(Cmatrix)$values
  lower <- max(lambda[lambda<0]) + 1e-8
  upper <- min(lambda[lambda>0]) - 1e-8
  par <- opt.iter.UV(fnU = fnU,fnV = fnV,V = V,opt.iter = control$opt.iter,
                     opt.iter.method = control$opt.iter.method,Xstar = Xstar,
                     Ystar = Ystar,D = D,k = k,m = m,Cmatrix = Cmatrix,
                     lower = lower,upper = upper,opt.iter.show = control$opt.iter.show)
  llvalues <- par[[2]]
  par <- par[[1]]
  nparV <- length(par) - 1
  opt <- optim(par = par, fn = fn, gr = NULL, Xstar = Xstar, 
               Ystar = Ystar, D = D, k = k, m = m, Cmatrix = Cmatrix,
               lower = c(lower,rep(-1e+8,nparV)), upper = c(upper,rep(1e+8,nparV)),
               method = "L-BFGS-B", control = control$optim, hessian = control$hessian)
  if(control$opt.iter.show) cat(opt$value," ")
  llvalues <- c(llvalues,opt$value)
  V <- par2V(opt$par[-1],k)
  rho <- opt$par[1]
  U <- chol2inv(chol(diag(m)-rho * Cmatrix))
  gls <- glsfit_s(Xstar, Ystar,D, V,U, onlycoef = FALSE)
  qrinvtMX <- qr(gls$invtMX)
  R <- qr.R(qrinvtMX)
  Qty <- qr.qty(qrinvtMX, gls$invtMY)
  vcov <- tcrossprod(backsolve(R, diag(1, ncol(gls$invtMX))))
  res <- NULL
  fitted <- Xstar%*%gls$coef
  rank <- qrinvtMX$rank
  converged <- opt$convergence == 0 | abs(llvalues[length(llvalues)]-llvalues[length(llvalues)-1]) < control$lltol
  
  c(list(coefficients = gls$coef, vcov = vcov, V = V, rho = rho,residuals = res, 
         fitted.values = fitted, df.residual = nall - rank - length(par), iter_llvalues = llvalues,
         rank = rank, logLik = opt$value, converged = converged, par = opt$par), if (!is.null(opt$hessian)) 
             list(hessian = opt$hessian, se.rho = sqrt(-solve(opt$hessian)[1,1])), 
    list(niter = opt$counts[[2]], control = control))
}
smvmeta.reml <- function (Xstar, Ystar, D = D, Cmatrix = Cmatrix, k, m, p, bscov, 
                        nall,control, ...){
  fn <- remlprof.fn
  fnU <- remlprof.fnU
  fnV <- remlprof.fnV
  V <- control$initial_V#as.matrix(Matrix::nearPD(xpndMat(optV$par))$mat) 
  lambda <- 1/eigen(Cmatrix)$values
  lower <- max(lambda[lambda<0]) + 1e-8
  upper <- min(lambda[lambda>0]) - 1e-8
  par <- opt.iter.UV(fnU = fnU,fnV = fnV,V = V,opt.iter = control$opt.iter,
                     opt.iter.method = control$opt.iter.method,Xstar = Xstar,
                     Ystar = Ystar,D = D,k = k,m = m,Cmatrix = Cmatrix,
                     lower = lower,upper = upper,opt.iter.show = control$opt.iter.show)
  llvalues <- par[[2]]
  par <- par[[1]]
  nparV <- length(par) - 1
  opt <- optim(par = par, fn = fn, gr = NULL, Xstar = Xstar, 
               Ystar = Ystar, D = D, k = k, m = m, Cmatrix = Cmatrix,
               lower = c(lower,rep(-1e+8,nparV)), upper = c(upper,rep(1e+8,nparV)),
               method = "L-BFGS-B", control = control$optim, hessian = control$hessian)
  if(control$opt.iter.show) cat(opt$value," ")
  llvalues <- c(llvalues,opt$value)
  V <- par2V(opt$par[-1],k)
  rho <- opt$par[1]
  
  U <- chol2inv(chol(diag(m)-rho * Cmatrix))
  gls <- glsfit_s(Xstar, Ystar,D, V,U, onlycoef = FALSE)
  qrinvtMX <- qr(gls$invtMX)
  R <- qr.R(qrinvtMX)
  Qty <- qr.qty(qrinvtMX, gls$invtMY)
  vcov <- tcrossprod(backsolve(R, diag(1, ncol(gls$invtMX))))
  res <- NULL
  fitted <- Xstar%*%gls$coef
  rank <- qrinvtMX$rank
  converged <- opt$convergence == 0 | abs(llvalues[length(llvalues)]-llvalues[length(llvalues)-1]) < control$lltol
  
  c(list(coefficients = gls$coef, vcov = vcov, V = V, rho = rho,residuals = res, 
         fitted.values = fitted, df.residual = nall - rank - length(par), iter_llvalues = llvalues,
         rank = rank, logLik = opt$value, converged = converged, par = opt$par), if (!is.null(opt$hessian)) 
             list(hessian = opt$hessian, se.rho = sqrt(-solve(opt$hessian)[1,1])), 
    list(niter = opt$counts[[2]], control = control))
}
smvmeta.fit <- function (X, y, S, Cmatrix, offset = NULL, method = "reml", bscov = "unstr", 
                          control = list()){
  control <- do.call("smvmeta.control", control)
  y <- as.matrix(y)
  # nay <- is.na(y)
  k <- ncol(y) # k = q
  m <- nrow(y) # m = n
  p <- ncol(X)
  nall <- length(y)
  nk <- colnames(y)
  if (k > 1L && is.null(nk)) 
    nk <- paste("y", seq(k), sep = "")
  nm <- rownames(y)
  np <- colnames(X)
  if (control$inputna) {
    augdata <- inputna(y, S, inputvar = control$inputvar)
    y <- augdata[, seq(k)]
    S <- augdata[, -seq(k)]
    nay[nay] <- FALSE
  }
  if (!is.null(offset)) 
    y <- y - offset
  Xlist <- lapply(seq(m), function(i) diag(1, k) %x% matrix(X[i,],nrow = 1))
  Xstar <- do.call("rbind",Xlist)
  Ystar <- as.numeric(t(y))
  if (dim(S)[2] == k) 
    S <- inputcov(sqrt(S), control$Scor)
  Slist <- lapply(seq(m), function(i) xpndMat(S[i, ]))
  D <- lapply(seq(m), function(i) t(diag(1, m)[i,] %x% Slist[[i]]))
  D <- do.call("rbind",D)
  fun <- paste("smvmeta", method, sep = ".")
  fit <- do.call(fun, list(Xstar = Xstar, Ystar = Ystar, D = D, Cmatrix = Cmatrix,
                           k = k, m = m, p = p, bscov = bscov,nall = nall, control = control))
  if (!fit$converged) {
    warning("Not convergency")
  }
  iter_llvalues <- fit$iter_llvalues
  names(iter_llvalues) <- c("U",rep(c("V","U"),length(iter_llvalues)/2-1),"UV")
  fit$iter_llvalues <- iter_llvalues
  
  U <- chol2inv(chol(diag(m)-fit$rho*Cmatrix))
  V <- fit$V
  H <- U %x% V
  invSIGMA <- chol2inv(chol(H + D))
  ytilde <- Ystar - fit$fitted.values
  fit$xi <- H %*% invSIGMA %*% ytilde
  fit$epsilon <- ytilde - fit$xi
  
  fit$fitted.values.spatial <- fit$fitted.values + fit$xi
  fit$method <- method
  fit$bscov <- bscov
  fit$offset <- offset
  fit$dim <- list(k = k, m = m, p = p)
  fit$df <- list(nobs = nall - (method == "reml") * 
                   fit$rank, df = nall - fit$df.residual, fixed = fit$rank, 
                 random = ifelse(method == "fixed", 0, nall - fit$rank - 
                                   fit$df.residual))
  fit$D <- D
  fit$lab <- list(k = nk, p = np)
  temp <- as.numeric(fit$fitted.values)
  fit$fitted.values <- matrix(temp, m, k, byrow = TRUE)
  temp <- as.numeric(fit$fitted.values.spatial)
  fit$fitted.values.spatial <- matrix(temp, m, k, byrow = TRUE)
  temp <- as.numeric(fit$epsilon)
  fit$epsilon <- matrix(temp, m, k, byrow = TRUE)
  temp <- as.numeric(fit$xi)
  fit$xi <- matrix(temp, m, k, byrow = TRUE)
  if (!is.null(offset)) {
    y <- y + offset
    fit$fitted.values <- fit$fitted.values + offset
  }
  if (method != "fixed") 
    dimnames(fit$V) <- list(nk, nk)
  if (k == 1L) {
    names(fit$coefficients) <- np
    dimnames(fit$vcov) <- list(np, np)
    fit$fitted.values <- drop(fit$fitted.values)
    fit$residuals <- drop(y - fit$fitted.values)
    names(fit$residuals) <- names(fit$fitted.values) <- nm
    names(fit$fitted.values.spatial) <- names(fit$epsilon) <- names(fit$xi) <- nm
  }
  else {
    fit$coefficients <- matrix(fit$coefficients, p, k, dimnames = list(np,nk))
    rownames(fit$vcov) <- colnames(fit$vcov) <- paste(rep(nk, 
                                                          each = p), rep(np, k), sep = ".")
    fit$residuals <- y - fit$fitted.values
    dimnames(fit$residuals) <- dimnames(fit$fitted.values) <- list(nm,nk)
    dimnames(fit$fitted.values.spatial) <- dimnames(fit$epsilon) <- dimnames(fit$xi) <- list(nm,nk)
  }
  fit
}


# LR test only for ml method
llrtest <- function(model,modelref){
  class(model) <- class(modelref) <- "mvmeta"
  lrstat <- -2*(logLik(modelref)-logLik(model))[1]
  df <- attr(logLik(model),"df")-attr(logLik(modelref),"df")
  pvalue <- 1-pchisq(lrstat,df)
  lr <- c(lrstat,df,pvalue)
  names(lr) <- c("chi.square","df","p_value")
  lr
}



# only for mvmeta
predmvmeta <- function(object){
  fitvalue <- object$fitted.values
  obsy <- object$model[[1]]
  psi <- object$Psi
  SIGMA <- object$S
  res <- sapply(1:nrow(obsy), function(x){
    invs <- chol2inv(chol(xpndMat(SIGMA[x,]) + psi))
    fitvalue[x,] + psi %*% invs %*% (obsy[x,] - fitvalue[x,])
  })
  res <- t(res)
  rownames(res) <- rownames(obsy)
  colnames(res) <- colnames(obsy)
  res
}
# only for mvmeta; mvmeta and smvmeta have the same results
GetH <- function(m,mref=NULL) {
  # HETEROGENEITY AND IC STATS
  q <- qtest(m)
  het <- c(q$Q[1],q$df[1],q$pvalue[1],(q$Q[1]-q$df[1])/q$Q[1]*100,AIC(m),BIC(m))
  # LR TEST (ONLY FOR META-REGRESSION)
  if(!is.null(mref)) {
    lrstat <- -2*(logLik(mref)-logLik(m))[1]
    df <- attr(logLik(m),"df")-attr(logLik(mref),"df")
    pvalue <- 1-pchisq(lrstat,df)
    lr <- c(lrstat,df,pvalue)
  }
  # WALD TEST (ONLY FOR META-REGRESSION)
  if(!is.null(mref)) {
    coef <- coef(m)[-grep("Int",names(coef(m)))]
    vcov <- vcov(m)[-grep("Int",names(coef(m))),-grep("Int",names(coef(m)))]
    waldstat <- coef%*%solve(vcov)%*%coef
    df <- length(coef)
    pvalue <- 1-pchisq(waldstat,df)
    wald <- c(waldstat,df,pvalue)
  }
  # RESULTS
  if(!is.null(mref)) {
    return(c(het,lr,wald))
  } else return(het) 
}

# LR test for rho
rhotest <- function(mvmodel,smvmodel){
  class(smvmodel) <- "mvmeta"
  logrr <- (logLik(smvmodel) - logLik(mvmodel)) * 2
  a <- c(smvmodel$rho,logrr,1-pchisq(logrr,df = 1)) 
  names(a) <- c("rho","chiquare","pvalue")
  a
}

