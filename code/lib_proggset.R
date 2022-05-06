#   Score statistics for Cox PH model
#   It was adopted from samr's subrutine cox.func
#   values: tt = socre, numer = r, sd = ss
#   tt = r/(ss+s0), when s0=0 it is normal distribution
#   tt^2 is chi-sqare with 1 degree of freedom

samr.scor <- function (x, y, censoring.status, s0 = 0) 
{
  Dn <- sum(censoring.status == 1)
  Dset <- c(1:ncol(x))[censoring.status == 1]
  ind <- matrix(0, ncol(x), Dn)
  for (i in 1:Dn) {
    ind[y > y[Dset[i]] - 1e-08, i] <- 1/sum(y > y[Dset[i]] - 1e-08)
  }
  ind.sums <- rowSums(ind)
  x.ind <- x %*% ind
  numer <- x %*% (censoring.status - ind.sums)
  temp <- (x * x) %*% ind.sums - rowSums(x.ind * x.ind)
  p <- which(temp > 0)
  tt <- matrix(NA, length(numer), 1)
  rownames(tt) <- rownames(numer)
  my.sd <- tt
  my.sd[p] <-  sqrt(temp[p])
  tt[p] <- numer[p]/(my.sd[p] + s0)
  return(list(tt = tt, numer = numer, sd = my.sd))
}


# Composite score calculation
# X = n by p matrix with each row representating a sample and each column, a gene
#     Note X should be base 2 log transformed.
# survdays and censored are survival time in days and censoring status
# Values:
#    score = Risk scores between 0 and 1 with 0.5 being the median
calScore <- function(X, survdays, censored, do.scale=T)
{
  if(do.scale) X <- scale(X)
  fit <- samr.scor(x=t(X), y=survdays, censoring.status=censored)
  col.rm <- which(is.na(fit$tt) | is.nan(fit$tt) | is.infinite(fit$tt))
  if(length(col.rm) > 0)
  {
    # data of columns col.rm cannot be fit. These columns has to be removed
    X <- X[, -col.rm]
    fit$tt <- fit$tt[-col.rm]
  }
  LM <- X %*% fit$tt
  m <- median(LM)
  return(list(score=1/(1 + exp(-(LM-m))), col.rm=col.rm))
}

calScore4cv <- function(idx, X, survdays, censored, do.scale=T, quantileValues=c(0.25, 0.75))
{
  if(do.scale) X <- scale(X)
  test.dat <- X[idx, ]
  X <- X[-idx, ]
  fit <- samr.scor(x=t(X), y=survdays[-idx], censoring.status=censored[-idx])
  col.rm <- which(is.na(fit$tt) | is.nan(fit$tt) | is.infinite(fit$tt))
  if(length(col.rm) > 0)
  {
    # data of columns col.rm cannot be fit. These columns has to be removed
    X <- X[, -col.rm]
    fit$tt <- fit$tt[-col.rm]
  }
  LM <- X %*% fit$tt
  crit <- quantile(LM, probs = quantileValues)
  pred.scor <- sum(test.dat*fit$tt)
  if(pred.scor < crit[1]) rtf <- 0
  else if(pred.scor >= crit[2]) rtf <- 1
  else rtf <- NA
  # cat(idx, "\t=>", rtf1, rtf2, pred.scor, crit, "\n")
  return(rtf)
}

# Calculate groups based on CV
# X = n by p matrix with each row representating a sample and each column, a gene
#     Note X should be base 2 log transformed.
# survdays and censored are survival time in days and censoring status
# quantileValues = a vector of two elements, the quantiles for the
#                  lower and upper thresholds to dichotomize the predicted score.
#                  The default is c(0.25, 0.75)
# Values: res = a vector of n risk scores : 0 = low risk, 1 = high risk
#         HR = Hazard ratio
#         p.value = p-value of log-rank test
calScore.cv <- function(X, survdays, censored, do.scale=T, quantileValues=c(0.25, 0.75))
{
  nP <- dim(X)[1]
  res <- rep(NA, nP)
  for( k in 1:nP) res[k] <- calScore4cv(idx=k, X=X, survdays=survdays, censored=censored, do.scale=do.scale, quantileValues=quantileValues)
  dat <- data.frame(tim=survdays, cens=censored, risk=res)
  fit <- summary(coxph(Surv(tim, cens) ~ risk, data=dat))
  return(list(res=res, HR=fit$coefficients[2], p.value=fit$sctest[3]))
}

# Calculate the scores and risk groups using the same data for determining the model
# X = n by p matrix with each row representating a sample and each column, a gene
#     Note X should be base 2 log transformed.
# survdays and censored are survival time in days and censoring status
# quantileValues = a vector of two elements, the quantiles for the
#                  lower and upper thresholds to dichotomize the predicted score.
#                  The default is c(0.25, 0.75)
# Values: res = a vector of n risk scores : 0 = low risk, 1 = high risk
#         stat = chi-square test statistics of long-rank test
#         p.value = p-value of log-rank test (Biased)
calScore.resub <- function(X, survdays, censored, do.scale=T, quantileValues=c(0.25, 0.75))
{
  temp <- calScore(X=X, survdays=survdays, censored=censored, do.scale=do.scale)
  res <- data.frame(tim=survdays, cens=censored, score=temp$score,
                    risk=rep(NA, length(temp$score)))
  crit <- quantile(res$score, probs = quantileValues, na.rm = T)
  res$risk[res$score < crit[1]] <- 0
  res$risk[res$score >= crit[2]] <- 1
  fit <- survdiff(Surv(tim, cens) ~ risk, data=res)
  return(list(res=res, stat=fit$chisq, p.value=pchisq(fit$chisq, df = 1, lower.tail = F), col.rm=temp$col.rm))
}



# Calculate the scores and risk groups using the same data for determining the model
# But p-value is calculated using permutation
# X = n by p matrix with each row representating a sample and each column, a gene
#     Note X should be base 2 log transformed.
# survdays and censored are survival time in days and censoring status
# quantileValues = a vector of two elements, the quantiles for the
#                  lower and upper thresholds to dichotomize the predicted score.
#                  The default is c(0.25, 0.75)
# Values: res = a vector of n risk scores : 0 = low risk, 1 = high risk
#         stat = chi-square test statistics of long-rank test
#         p.value = p-value of log-rank test (Biased)
#         perm.p.value = permutation based p-value of log-rank test (Unbiased)
calScore.perm <- function(X, survdays, censored, do.scale=T, quantileValues=c(0.25, 0.75), nPerm=100, seed=100){
  if(do.scale) X <- scale(X)
  obs <- calScore.resub(X=X, survdays=survdays, censored=censored, do.scale=F, quantileValues=quantileValues)
  res <- matrix(NA, nPerm, 2)
  colnames(res) <- c("chisq", "p.value")
  nP <- length(survdays)
  if(!is.null(seed)) set.seed(seed)
  for(k in 1:nPerm)
  {
    p <- sample(nP)
    tim <- survdays[p]
    cen <- censored[p]
    r <- calScore.resub(X=X, survdays=tim, censored=cen, do.scale=F, quantileValues=quantileValues)
    # if(length(r$col.rm) > 0) cat(k, "\tstat = ", r$stat, r$col.rm, "\n")
    res[k, 1] <- r$stat
    res[k, 2] <- r$p.value
  }
  return(list(res=res, stat=obs$stat, p.value=obs$p.value, perm.p.value = mean(res[, 1] > obs$stat)))
}

#One survival plot for all genes one cancer


writeLassoOutput <- function(winners, fixlcfit) {
  result <- c()
  result$Var <- c(names(winners)[fixlcfit$vars])
  result$ENSID <- c(names(winners)[fixlcfit$vars])
  result$Coef <- fixlcfit$coef0
  result$'Z-score' <- fixlcfit$zscore0
  result$'P-value' <- fixlcfit$pv
  result$LowConfPt <- fixlcfit$ci[, 1]
  result$UpConfPt <- fixlcfit$ci[, 2]
  result$LowTailArea <- fixlcfit$tailarea[, 1]
  result$UpTailArea <- fixlcfit$tailarea[, 2]
  return(as.data.frame(result))
}

#One survival plot for all genes one cancer
plotCoxSurvival <- function(survexp, survdays, censored, lower_threhold, upper_threhold, ...) {
  r1 <- calScore.perm(survexp, survdays, censored, quantileValues=c(lower_threhold, upper_threhold), nPerm=1000, seed=100)
  r <- calScore.resub(survexp, survdays, censored, quantileValues=c(lower_threhold, upper_threhold))
  grp <- factor(r$res$risk, levels=0:1, labels = c("low", "high"))
  fit <- survfit(Surv(survdays, censored) ~ grp)
  fit1 <- summary(coxph(Surv(survdays, censored) ~ grp))
  
  filename=paste0("Expsurv_", canType,"_Permutation_",  stype, ".png")
  png(filename, width = 2000, height = 2000,units = 'px',  res = 400)
  
  plot(fit, col = 2:3, lty=2:1, xlab = "Day", ylab = paste(canType,": ",ylabname), main=sprintf("Composite Exprs Score cut at %4.2f, %4.2f", lower_threhold, upper_threhold))
  legend("bottomleft", legend = c("Low Risk", "High Risk"), col = 2:3, lty=2:1, inset = 0.02)
  txt1 <- sprintf("HR = %4.2f", fit1$coefficients[2])
  txt2 <- ifelse(fit1$sctest[3]<.001, "p-value (Biased) < 0.001", sprintf("p-value (biased) = %5.3f", fit1$sctest[3]))
  txt3 <- ifelse(r1$perm.p.value<.001, "p-value (Permutation) < 0.001", sprintf("p-value (Permutation) = %5.3f", r1$perm.p.value))
  legend("topright", legend = c(txt1, txt2, txt3), bty = 'n')
  dev.off()
}

#Convert gene symbols to Entrez IDs
SYM2EG <- function(sym) {
	if (!require('org.Hs.eg.db')) {
		warning('org.Hs.eg.db is not pre-installed. Original IDs are returned\n.')
		egs <- sym
	} else {
		if (length(sym)>length(unique(sym))) warning('input list includes duplicates!\n')
		egs <- unlist(mget(intersect(sym,keys(revmap(org.Hs.egSYMBOL))),revmap(org.Hs.egSYMBOL)))
		original.geneIDs <- sym
		if (length(setdiff(sym,keys(revmap(org.Hs.egSYMBOL))))>0) {
			cat('non-mapped gene ids: ',setdiff(sym,keys(revmap(org.Hs.egSYMBOL))),'\n')
			write.table(setdiff(sym,keys(revmap(org.Hs.egSYMBOL))),'non-mapped IDs.txt',row.names=F,col.names=F,quote=F)
		}
		egs <- egs[as.character(sym)]
		names(egs)[is.na(egs)] <- original.geneIDs[is.na(egs)]
		egs[is.na(egs)] <- original.geneIDs[is.na(egs)]	
	}
	return(egs)
} 

