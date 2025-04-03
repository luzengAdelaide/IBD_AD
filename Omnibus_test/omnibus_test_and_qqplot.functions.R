
get_05 <- function(vect){quantile(vect, p=.05)}
get_95 <- function(vect){quantile(vect, p=.95)}
get_50 <- function(vect){quantile(vect, p=.50)}
get_02.5 <- function(vect){quantile(vect, p=.025)}
get_97.5 <- function(vect){quantile(vect, p=.975)}
get_00.5 <- function(vect){quantile(vect, p=.005)}
get_99.5 <- function(vect){quantile(vect, p=.995)}

permute <- function(pheno.permute, data.permute, var_list.permute, covs.permute){

  data.permute$pheno.permute.r <- sample(data.permute[,pheno.permute], replace=F)

  ps <- rep(NA, length(var_list.permute))

  for (i in 1:length(var_list.permute)){

    data.permute$temp.var <- data.permute[,var_list.permute[i]]

    f <- as.formula(paste("pheno.permute.r ~ temp.var + " , paste(covs.permute, collapse=" + ")))

    lm1 <- lm(f, data=data.permute)

    ps[i] <- summary(lm1)$coefficients[2,4]

  }

  return(ps)

}


create_plot <- function(variable_group, o.real, nsims, label, pheno.cp, data.cp, var_list.cp, covs.cp, cohort.omni){

  o.mat <- matrix(NA,nsims,length(var_list.cp))

  for(i in 1:nsims){

    o = -log10(sort(permute(pheno.cp, data.cp, var_list.cp, covs.cp),decreasing=F))
    e = -log10( 1:length(o)/length(o) )
    o.mat[i,] <- o

  }

  write.table(o.mat, paste0("sim_results/",variable_group,".",nsims,".", label,".", pheno.cp,".", cohort.omni,".02AUG2019.txt"), quote=F, sep="\t", row.names=F)

  o.05 <- apply(o.mat,2,get_05)
  o.95 <- apply(o.mat,2,get_95)
  o.50 <- apply(o.mat,2,get_50)
  o.02.5 <- apply(o.mat,2,get_02.5)
  o.97.5 <- apply(o.mat,2,get_97.5)
  o.00.5 <- apply(o.mat,2,get_00.5)
  o.99.5 <- apply(o.mat,2,get_99.5)

  o.max <- apply(o.mat,2,max)
  o.min <- apply(o.mat,2,min)

  print(o.05)
  print(0.95)

  print(e)
  print(rev(e))
  print(length(e))
  print(length(o.real))

  print(e)
  print(o.50)

  pdf(paste0("plots/",variable_group,".sim.qqplot.",nsims,".",label,".",pheno.cp,".pdf"))

  plot(e,o.real, type="n", xlim=c(0,max(e)), ylim=c(0,max(o.real)), main=paste0(pheno.cp, " " , cohort.omni," ", variable_group),xlab="Expected", ylab="Observed")

  polygon(c(e, rev(e)), c(o.97.5, rev(o.02.5)), col = "light grey", border = NA)

  polygon(c(e, rev(e)), c(o.95, rev(o.05)),col = "dark grey", border = NA)

  lines(lowess(e,o.50, f=.15))

  points(e, o.real, pch=20)

  legend(legend=c("95% CI","90% CI"),"topleft", col=c("light grey", "dark grey"), pch=15, bty="n", pt.cex=2)

  dev.off()
}

omni.test <- function(realpvals,variable_group, nsims, label, pheno.omni, regs.ot, cohort.ot){

                                        #  realpvals <- regs1$pval
  realpvals <- regs.ot$pval

  q.real <- -2*sum(log(realpvals))

  sims <- read.table(paste0("sim_results/",variable_group,".",nsims,".",label,".", pheno.omni,".", cohort.ot,".02AUG2019.txt"), as.is=T, header=T,sep="\t")

  sims <- as.matrix(sims)

  q.sims <- rep(NA, nsims)

  for(i in 1:nsims){

    q.temp <- -2*sum(log(10^-sims[i,]))
    q.sims[i] <- q.temp

  }

  test <- q.real < q.sims

  p.omni <- sum(test)/nsims

  return(p.omni)

}


omni_testing <- function(pheno.omni, var.omni , covs.omni, data.omni, label.omni, nsims.omni, variable_group.omni, cohort.omni_testing){

  regs1 <- matrix.regs(pheno.omni, var.omni, covs.omni, data.omni ,"linear")
  print(summary(regs1))
  write.table(regs1, paste0("results/all.results.", variable_group.omni,".regs1.", nsims.omni, ".", label.omni,".", pheno.omni,".",cohort.omni_testing,".02AUG2019.txt"), quote=F, row.names=F, sep="\t")
  temp <- regs1[order(regs1$pval),][1,]
  var.top <- as.character(temp$ind.var)
  o.real = -log10(sort(regs1$pval,decreasing=F))
  e = -log10( 1:length(o.real)/length(o.real) )
  create_plot(variable_group.omni, o.real, nsims.omni, label.omni, pheno.omni, data.omni, var.omni, covs.omni, cohort.omni_testing)
  omni.p <- omni.test(regs1$pval,variable_group.omni,nsims.omni, label.omni, pheno.omni, regs1, cohort.omni_testing)

  return(c(variable_group.omni, var.top, temp$pval, temp$fdr, length(var.omni), omni.p))

}
