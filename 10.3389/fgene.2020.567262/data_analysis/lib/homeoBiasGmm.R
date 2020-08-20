library(distr)

# Calculate indicator Z
calcIndicator <- function(y, nc, ratio, param.beta) {
  numerator <- matrix(0, nrow=length(y), ncol=nc)
  for(i in 1:nc) {
    numerator[, i] <- ratio[i] * dbeta(y, param.beta[i, 1], param.beta[i, 2])
  }
  if (sum(is.na(numerator))) {
    cat("Warning: indicator Z includes 'NA'.")
    num.row.na <- apply(is.na(numerator), 1, sum)
    is.row.na <- num.row.na >= 1
    numerator[is.row.na, ] <- 1
  }
  z <- numerator / apply(numerator, 1, sum)
  return(z)
}

# Calculate Log-Likelihood
calcLLObs <- function(y, nc, ratio) {
  return(function(param.beta.vec) {
    param.beta.mat <- matrix(param.beta.vec, nc) 	    
    sumL <- 0
    for(i in 1:nc) {
      sumL <- sumL + ratio[i] * dbeta(y, param.beta.mat[i, 1], param.beta.mat[i, 2])
    }
    logL <- log(sumL)
    return(sum(logL))
  })
}

# Calculate Log-Likelihood for 'complete' data
calcLLComp <- function(y, nc, ratio, z) {
  return(function(param.beta.vec) {
    param.beta.mat <- matrix(param.beta.vec, nc)       
    sumL <- 0
    for(i in 1:nc) {
      sumL <- sumL + sum(z[, i] * log(ratio[i] * dbeta(y, param.beta.mat[i, 1], param.beta.mat[i, 2])))
    }
    return(sumL)
  })
}

# Calculate ratio of components
calcRatio <- function(y, z) {  
  ratio.new <- rep(0, ncol(z))  
  for(i in 1:ncol(z)) {  
    ratio.new[i] <- sum(z[, i]) / length(y)
  }
  return(ratio.new)
}

# Update parameter of beta
updateParam.beta <- function(y, z, param.beta) {  
  param.beta.new <- matrix(0, nrow=nrow(param.beta), ncol=ncol(param.beta))
  for(i in 1:nrow(param.beta)) {    
    param.beta.new[i, 1] <- igamma(digamma(param.beta[i, 1] + param.beta[i, 2]) + sum(z[, i] * log(y)) / sum(z[, i]))
    param.beta.new[i, 2] <- igamma(digamma(param.beta[i, 1] + param.beta[i, 2]) + sum(z[, i] * log(1-y)) / sum(z[, i]))  
  }  
  return(param.beta.new)
}

# Estimate BMM
fitBMM <- function(y, ratio, param.beta, method='optim', niter=200 ,tol=1E-5, display=T) {
  # check initial inputs  
  if (ncol(param.beta) != 2) stop("ncol of param.beta must be two shape paramters of beta distribution (alpha and beta).")  
  if (length(ratio) != nrow(param.beta)) stop("length of ratio must be same with ncol of param.beta.")  
  if (method != 'update' & method != 'optim') stop("method must be 'update' or 'optim'.")
  
  nc <- nrow(param.beta) # number of components  

  if (display) {
    cat('Number of data: ', length(y), '\n')	
    cat('Number of components: ', nc, '\n')	
    cat('Ratio: ', ratio, '\n')
    cat('Alpha: ', param.beta[,1], '\n')
    cat('Beta: ', param.beta[,2], '\n')	
    cat('Method: ', method, '\n')    	
    cat('Max interation: ', niter, '\n')	
    cat('TOL(divergence of log-likelihood): ', tol, '\n\n')
    cat('i,', 'is.convergence,', 'Observed Log-likelihood,', 'Complete Log-likelihood,', 'Ratio,', 'Alpha,', 'Beta', '\n')
  }
  
  # generate history
  LL.obs.history <- rep(NA, niter+1)
  LL.comp.history <- rep(NA, niter+1)
  ratio.history <- matrix(NA, ncol=nc, nrow=niter+1)
  param.history <- matrix(NA, ncol=nc * 2, nrow=niter+1)
  
  # initial values
  param.beta.vec <- c(param.beta[, 1], param.beta[, 2])
  z <- calcIndicator(y, nc, ratio, param.beta)
  LL.obs.init <- (calcLLObs(y, nc, ratio))(param.beta.vec)
  LL.comp.init <- (calcLLComp(y, nc, ratio, z))(param.beta.vec)
  if (display){
    cat(1, ': ', '-', LL.obs.init, LL.comp.init, ratio, param.beta.vec, '\n')
  }
  LL.obs.history[1] <- LL.obs.init
  LL.comp.history[1] <- LL.comp.init
  ratio.history[1, ] <- ratio
  param.history[1, ] <- param.beta.vec
  
  is.conv <- F
  for(i in 1:niter) {	
    
    # M-step
    ratio <- calcRatio(y, z)
    if (method == 'update') {
      param.beta <- updateParam.beta(y, z, param.beta)
      param.beta.vec <- c(param.beta[, 1], param.beta[, 2])
      LL.obs <- (calcLLObs(y, nc, ratio))(param.beta.vec)
      LL.comp <- (calcLLComp(y, nc, ratio, z))(param.beta.vec)
      if (display) {
        cat(i+1, ': ', '-', LL.obs, LL.comp, ratio, param.beta.vec, '\n')
      }
    } else if (method == 'optim') {
      res.opt <- optim(param.beta.vec, calcLLComp(y, nc, ratio, z), control=list(fnscale=-1))
      #res.opt <- optim(param.beta.vec, calcLLObs(y, nc, ratio), control=list(fnscale=-1)) ##
      param.beta.vec <- res.opt$par
      param.beta <- matrix(param.beta.vec, nc)
      LL.comp <- res.opt$value
      LL.obs <- (calcLLObs(y, nc, ratio))(param.beta.vec)
      #LL.obs <- res.opt$value ##
      #LL.comp <- (calcLLComp(y, nc, ratio, z))(param.beta.vec) ##
	    if (display) {
        cat(i+1, ': ', res.opt$convergence, LL.obs, LL.comp, ratio, param.beta.vec, '\n')
	    }
    }
    
    LL.obs.history[i+1] <- LL.obs
    LL.comp.history[i+1] <- LL.comp
    ratio.history[i+1, ] <- ratio
    param.history[i+1, ] <- param.beta.vec
    
    # Terminal condition
    #if (abs(LL.comp.history[i]-LL.comp.history[i+1]) < tol) { 
    if (abs(LL.obs.history[i]-LL.obs.history[i+1]) < tol) { ##
      is.conv <- T
      break
    }
    
    # E-step
    z <- calcIndicator(y, nc, ratio, param.beta)
  }
  
  rownames(param.beta) <- 1:nc  
  colnames(param.beta) <- c('param1', 'param2')  
  res <- list(LL.obs.history[i+1], ratio, param.beta, is.conv, i, LL.obs.history[1:(i+1)], LL.comp.history[1:(i+1)], ratio.history[1:(i+1), ], param.history[(1:i+1), ])
  names(res) <- c('LL.obs', 'ratio', 'param.beta', 'is.conv', 'iter', 'LL.obs.history', 'LL.comp.history', 'ratio.history', 'param.history')
  if (display) {
    cat('\n#Result\n')	
    cat('Log-likelihood: ', res$LL.obs, '\n')	
    cat('Ratio: ', res$ratio, '\n')
    cat('Alpha: ', res$param.beta[, 1], '\n')
    cat('Beta: ', res$param.beta[, 2], '\n')
    cat('is.conv: ', res$is.conv, '\n')	
    cat('Number of iteration: ', res$iter, '\n')  
  }  
  return(res)
}

# draw histogram of input data and estimated BMM model
plotBMM <- function(y, ratio, param.beta, hist.by=0.05, my.xlab='x', my.ylab='Density', my.main='', my.lwd=2, pdf.name='') {
  yhist <- hist(y, xlim=c(0, 1), breaks=seq(0, 1, hist.by), freq=F)
  #ymax <- round(max(yhist$density)+1, 0)
  ymax <- round(max(yhist$density) * 1.2, 0)
  plot.new()
  if(pdf.name != '') {
    pdf(pdf.name)
  }
  hist(y, xlim=c(0, 1), breaks=seq(0, 1, hist.by), freq=F, ylim=c(0, ymax), xlab=my.xlab, ylab=my.ylab, main=my.main)
  par(new=T)
  xp.by <- 5E-3
  xp <- seq(xp.by, 1-xp.by, xp.by)
  yp <- rep(0,length(xp))
  for(i in 1:length(ratio)) {
    yp <- yp + ratio[i] * dbeta(xp, param.beta[i, 1], param.beta[i, 2])
    curve(ratio[i] * dbeta(x, param.beta[i, 1], param.beta[i, 2]), ylim=c(0, ymax), xlab='', ylab='', main='', col='blue', lwd=my.lwd)
    par(new=T)
  }
  plot(xp, yp, type='l', xlim=c(0, 1), ylim=c(0, ymax), xlab='', ylab='', main='', col='red', lty="dashed", lwd=my.lwd)
  if(pdf.name != '') {
    dev.off()
  }
}

# get alpha and beta parameters of beta distribution from mean and variance of input data y
getParamMV <- function(y) {
  m <- mean(y)
  v <- var(y) / length(y) * (length(y) - 1)  
  alpha <- m * (m * (1-m) / v - 1)
  beta <- (1-m) * (m * (1-m) / v - 1)
  return(c(alpha, beta))
}

# calculate difference between likelihood ratio (LR) of two beta-distributions and given LR
calcDiffLR <- function(x, ratio, param.beta, lr=1) {
  if (length(ratio) != 2) stop("length of ratio must be two.")
  if (nrow(param.beta) != 2) stop("nrow of param.beta must be two.")
  if (ncol(param.beta) != 2) stop("ncol of param.beta must be two.")
  likehood1 <- ratio[1] * dbeta(x, param.beta[1, 1], param.beta[1, 2])
  likehood2 <- ratio[2] * dbeta(x, param.beta[2, 1], param.beta[2, 2])
  return(lr - likehood1 / likehood2)
}

# calculate value satisfying given likelihood ratio
getValsGivenLR <- function(ratio, param.beta, lr=1, init.range=c(0, 0.4, 0.6, 1)) {
  x1 <- uniroot(calcDiffLR, c(init.range[1], init.range[2]), ratio, param.beta, lr)$root
  x2 <- uniroot(calcDiffLR, c(init.range[3], init.range[4]), ratio, param.beta, lr)$root
  return(c(x1, x2))
}

# calculate log-likelihood (logL) of two beta-distributions and log-likelihood ratio (logLR)
getLogLR <- function(df.val, ratio, param.beta) {
  if (length(ratio) != 2) stop("length of ratio must be two.")
  if (nrow(param.beta) != 2) stop("nrow of param.beta must be two.")
  if (ncol(param.beta) != 2) stop("ncol of param.beta must be two.")
  df.val.tmp <- as.matrix(df.val)
  logL.mat1 <- matrix(0, nrow=nrow(df.val.tmp), ncol=ncol(df.val.tmp))
  logL.mat2 <- matrix(0, nrow=nrow(df.val.tmp), ncol=ncol(df.val.tmp))
  for(i in 1:ncol(df.val.tmp)) {
    logL.mat1[, i] <- log(ratio[1] * dbeta(df.val.tmp[, i], param.beta[1, 1], param.beta[1, 2]))
    logL.mat2[, i] <- log(ratio[2] * dbeta(df.val.tmp[, i], param.beta[2, 1], param.beta[2, 2]))
  }
  logLR <- apply(logL.mat1, 1, sum) - apply(logL.mat2, 1, sum)
  like.info <- list(logL.mat1, logL.mat2, logLR)
  names(like.info) <- c("logL1", "logL2", "logLR") # log-likelihood (logL), log-likelihood ratio (logLR)
  return(like.info)
}

# get biased genes
getBiasGenes <- function(df, loglr, lr=1) {
  if(nrow(df) != length(loglr)) stop("nrow of df must be same with length of logLR.")
  
  is.bias <- loglr > log(lr)
  df.bias <- subset(df, is.bias)
  
  #is.each.low <- df.bias[, 2:ncol(df.bias)] < 0.5
  #is.each.high <- df.bias[, 2:ncol(df.bias)] > 0.5
  is.each.low <-as.matrix( df.bias[, 1] < 1/3)
  is.each.high <- as.matrix(df.bias[, 1] > 2/3)
  n.low <- apply(is.each.low, 1, sum)
  n.high <- apply(is.each.high, 1, sum)
  is.low.only <- n.low > 0 & n.high == 0
  is.high.only <- n.high > 0 & n.low == 0
  is.ambi <- n.low > 0 & n.high > 0
  df.bias.low <- subset(df.bias, is.low.only)
  df.bias.high <- subset(df.bias, is.high.only)
  df.bias.ambi <- subset(df.bias, is.ambi)
  
  n.all <- nrow(df)
  n.bias <- nrow(df.bias)
  n.bias.low <- nrow(df.bias.low)
  n.bias.high <- nrow(df.bias.high)
  n.bias.ambi <- nrow(df.bias.ambi)
  #stats <- c(n.all, n.bias, n.bias.low, n.bias.high, n.bias.ambi)
  #names(stats) <- c("n.all", "n.bias", "n.bias.low", "n.bias.high", "n.bias.ambi")
  stats <- data.frame('n.all'=n.all, 'n.bias'=n.bias, 'n.bias.low'=n.bias.low, 'n.bias.high'=n.bias.high, 'n.bias.ambi'=n.bias.ambi)
  
  out <- list(df.bias.low, df.bias.high, df.bias.ambi, stats)
  names(out) <- c("df.bias.low", "df.bias.high", "df.bias.ambi", "stats")
  return(out)
}
