to_numeric <- function(x) {as.numeric(as.character(x))}


matrix.r2s <- function(vars1, vars2, dataf, methodf) {

  k <- 1
  mat1 <- matrix(NA, length(vars1) * length(vars2),8)

  for (j in 1:length(vars2)){
    for (i in 1:length(vars1)){
      n <- dim(dataf[complete.cases(dataf[,c(vars1[i], vars2[j])]),])[1]

      temp <- cor.test(dataf[,vars1[i]], dataf[,vars2[j]], na.rm=T, method=methodf)
      p <- temp$p.value
      rho <- temp$estimate
      stat <- temp$statistic
      if(methodf=="pearson"){
        df <- temp$parameter
      } else {df <- NA}
      row <- c(vars1[i], vars2[j], signif(c(rho,stat, df,p),3), n, methodf)

      mat1[k,] <- row
      k <- k + 1
    }
  }
  out <- data.frame(mat1)
  names(out) <- c("var1", "var2", "rho", "stat", "df", "pval", "n", "method")
  out[,c("rho", "stat", "df", "pval", "n")] <- sapply(out[,c("rho", "stat", "df", "pval", "n")],function(x) as.numeric(as.character(x)))
  out$fdr <- signif(p.adjust(as.numeric(as.character(out$pval)), method="fdr"),3)
  return(out)

}



matrix.regs <- function(vars1, vars2,covariates, dataf, model_type){

  k <- 1

  covs2 <- paste(covariates,collapse=" + ")
  mat1 <- matrix(NA, length(vars1) * length(vars2),9)

  sum.covs <- sum(is.na(covariates))

  if (sum.covs > 0) {withcovs=0} else { withcovs=1}

  for (j in 1:length(vars2)){
    for (i in 1:length(vars1)){

      if(sd(dataf[,vars1[i]], na.rm=T)==0 | sd(dataf[,vars2[j]], na.rm=T)==0) {
        mat1[k,] <- c(vars1[i], vars2[j],rep(NA, 7))
      } else {
        if(withcovs==0) {formula1 <- as.formula(paste(vars1[i],"~", vars2[j],sep=""))
                         n <- dim(dataf[complete.cases(dataf[,c(vars1[i], vars2[j])]),])[1]
                       }
        if(withcovs==1) {formula1 <- as.formula(paste(vars1[i],"~", vars2[j], "+",covs2,sep=""))

                         n <- dim(dataf[complete.cases(dataf[,c(vars1[i], vars2[j], covariates)]),])[1]
                       }
        if (model_type=="linear") {

          row <- summary(lm(formula1, data=dataf))$coefficients[2,]
          df <- paste(summary(lm(formula1, data=dataf))$df[1:2], collapse=";")
        } else if (model_type=="logistic") {row <- summary(glm(formula1, data=dataf, family="binomial"))$coefficients[2,];

                                            df <- summary(glm(formula1, data=dataf, family="binomial"))$df[1]

                                          }

        outcovs <- paste(covariates, collapse=";")

        mat1[k,] <- c(vars1[i], vars2[j],signif(row,3), n, df, outcovs )
      }
      k <- k+1
    }
  }

  out <- data.frame(mat1)

  names(out) <- c("outcome", "ind.var", "beta", "se", "tstat", "pval", "n", "df", "covariates")

  out[,c("beta", "se", "tstat", "pval", "n")] <- sapply(out[,c("beta", "se", "tstat", "pval", "n")],function(x) as.numeric(as.character(x)))
  out$fdr <- p.adjust(as.numeric(as.character(out$pval)), method="fdr")
  if (model_type=="logistic") {
    out$df <- as.numeric(as.character(out$df))
    out$ORlow95 <- exp(out$beta - abs(qt(.05/2, df=out$df))*out$se)
    out$ORhigh95 <- exp(out$beta + abs(qt(.05/2, df=out$df))*out$se)
  } else if (model_type=="linear") {}
  return(out)
}


