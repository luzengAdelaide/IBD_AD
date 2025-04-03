rm(list=ls())

setwd("FILEPATH/omnibus-test/")
source("matrix_regressions.R")
source("omnibus_test_and_qqplot.functions.R")
system("mkdir sim_results results plots")

# trait: IBD/CD/UC

# read AD endophenotypes
omni_vars <- read.table("AD_endophenotypes.txt", as.is=T, header=T, sep="\t")$x

d0 <- read.table("risk.trait.profile", as.is=T, header=T, sep="\t")

data.seth <- d0
nsims=100
# read IBD/CD/UC risk score
phenos.seth <- "risk_score"
# read covariates
covs.seth <- c("age_death", "msex", "educ")
pheno.seth <- phenos.seth[1]
label.ax <- "risk.ibd"
variable_group <- "ibd"
cohort <- "ROSMAP"
row <- omni_testing(pheno.seth, omni_vars, covs.seth, data.seth, label.ax, nsims, variable_group, cohort)
out <- data.frame(t(row))
names(out) <- c("variable_group", "top_var", "pval", "fdr", "nVars", "OmniBus.p")
write.table(out, "results/global_ibd.risk.txt", quote=F, row.names=F, sep="\t")

