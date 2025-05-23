setwd("FILEPATH/mr")

############## Reverse Mendelian randomization analysis #############

##################
# Method 1: GSMR #
##################

# trait: ibd/uc/cd

################### AD as exposure, IBD/CD/UC as outcome ###################
rm(list=ls())
library("gsmr")
gsmr_data = read.table("use.ad.trait.gwas" ,sep="\t", header=T)

############# calculate ld
# Estimate LD correlation matrix using R
snp_coeff_id = scan("gsmr_ad_trait.xmat.gz", what="", nlines=1)
snp_coeff = read.table("gsmr_ad_trait.xmat.gz", header=F, skip=2)

# Match the SNP genotype data with the summary data
snp_id = Reduce(intersect, list(gsmr_data$SNP, snp_coeff_id))
gsmr_data = gsmr_data[match(snp_id, gsmr_data$SNP),]
snp_order = match(snp_id, snp_coeff_id)
snp_coeff_id = snp_coeff_id[snp_order]
snp_coeff = snp_coeff[, snp_order]

# Calculate the LD correlation matrix
ldrho = cor(snp_coeff, use="pairwise.complete.obs")

# Check the size of the correlation matrix and double-check if the order of the SNPs in the LD correlation matrix is consistent with that in the GWAS summary data
colnames(ldrho) = rownames(ldrho) = snp_coeff_id

dim(ldrho)

# Show the first 5 rows and columns of the matrix  
ldrho[1:5,1:5]

# GSMR analysis
bzx = gsmr_data$bzx    # SNP effects on the risk factor
bzx_se = gsmr_data$bzx_se    # standard errors of bzx
bzx_pval = gsmr_data$bzx_pval   # p-values for bzx
bzy = gsmr_data$bzy    # SNP effects on the disease
bzy_se = gsmr_data$bzy_se    # standard errors of bzy
bzy_pval = gsmr_data$bzy_pval    # p-values for bzy
n_ref = 63926    # Sample size of the reference sample
gwas_thresh = 5e-6    # GWAS threshold to select SNPs as the instruments for the GSMR analysis
single_snp_heidi_thresh = 0.01    # p-value threshold for single-SNP-based HEIDI-outlier analysis
multi_snp_heidi_thresh = 0.01    # p-value threshold for multi-SNP-based HEIDI-outlier analysis
nsnps_thresh = 10   # the minimum number of instruments required for the GSMR analysis
heidi_outlier_flag = T    # flag for HEIDI-outlier analysis
ld_r2_thresh = 0.05    # LD r2 threshold to remove SNPs in high LD
ld_fdr_thresh = 0.05   # FDR threshold to remove the chance correlations between the SNP instruments
gsmr2_beta = 0     # 0 - the original HEIDI-outlier method; 1 - the new HEIDI-outlier method that is currently under development 
gsmr_results = gsmr(bzx, bzx_se, bzx_pval, bzy, bzy_se, bzy_pval, ldrho, snp_coeff_id, n_ref, heidi_outlier_flag, gwas_thresh, single_snp_heidi_thresh, multi_snp_heidi_thresh, nsnps_thresh, ld_r2_thresh, ld_fdr_thresh, gsmr2_beta)    # GSMR analysis 
filtered_index=gsmr_results$used_index
cat("The estimated effect of the exposure on outcome: ",gsmr_results$bxy)
cat("Standard error of bxy: ",gsmr_results$bxy_se)
cat("P-value for bxy: ", gsmr_results$bxy_pval)
cat("Indexes of the SNPs used in the GSMR analysis: ", gsmr_results$used_index[1:5], "...")
cat("Number of SNPs with missing estimates in the summary data: ", length(gsmr_results$na_snps))
cat("Number of non-significant SNPs: ", length(gsmr_results$weak_snps))
cat("Number of SNPs in high LD ( LD rsq >", ld_r2_thresh, "): ", length(gsmr_results$linkage_snps))
cat("Number of pleiotropic outliers: ", length(gsmr_results$pleio_snps))

##### visulization
effect_col = colors()[75]
vals = c(bzx[filtered_index]-bzx_se[filtered_index], bzx[filtered_index]+bzx_se[filtered_index])
xmin = min(vals); xmax = max(vals)
vals = c(bzy[filtered_index]-bzy_se[filtered_index], bzy[filtered_index]+bzy_se[filtered_index])
ymin = min(vals); ymax = max(vals)
par(mar=c(5,4,2,2))
plot(bzx[filtered_index], bzy[filtered_index], pch=20, cex=0.8, bty="n", cex.axis=1.1, cex.lab=1.2,
     col=effect_col, xlim=c(xmin, xmax), ylim=c(ymin, ymax),
     xlab=expression(AD~(italic(b[zx]))),
     ylab=expression(trait~(italic(b[zy]))),
     main = bquote( atop(b[xy] ~'=-0.060' , 
                         'p-value = 1.84E-05')))
abline(0, gsmr_results$bxy, lwd=1.5, lty=2, col="dim grey")

nsnps = length(bzx[filtered_index])
for( i in 1:nsnps ) {
  # x axis
  xstart = bzx[filtered_index [i]] - bzx_se[filtered_index[i]]; xend = bzx[filtered_index[i]] + bzx_se[filtered_index[i]]
  ystart = bzy[filtered_index[i]]; yend = bzy[filtered_index[i]]
  segments(xstart, ystart, xend, yend, lwd=1.5, col=effect_col)
  # y axis
  xstart = bzx[filtered_index[i]]; xend = bzx[filtered_index[i]] 
  ystart = bzy[filtered_index[i]] - bzy_se[filtered_index[i]]; yend = bzy[filtered_index[i]] + bzy_se[filtered_index[i]]
  segments(xstart, ystart, xend, yend, lwd=1.5, col=effect_col)
}

#######################
##### TwoSampleMR #####
#######################
rm(list=ls())
library(TwoSampleMR)

ibd_file <- "FILEPATH/use.ad.ibd.gwas"
cd_file <- "FILEPATH/use.ad.cd.gwas"
uc_file <- "FILEPATH/use.ad.uc.gwas"

################### AD as exposure, IBD as outcome ###################
ad_exp_dat <- read_exposure_data(
  filename = ibd_file,
  sep = "\t",
  snp_col = "SNP",
  beta_col = "bzx",
  se_col = "bzx_se",
  effect_allele_col = "a1",
  other_allele_col = "a2",
  pval_col = "bzx_pval",
  samplesize_col = "bzx_n"
)
ad_exp_dat$exposure="AD"
ibd_out_dat <- read_outcome_data(
  filename = ibd_file,
  sep = "\t",
  snp_col = "SNP",
  beta_col = "bzy",
  se_col = "bzy_se",
  effect_allele_col = "a1",
  other_allele_col = "a2",
  pval_col = "bzy_pval",
  samplesize_col = "bzy_n"
)
ibd_out_dat$outcome="IBD"
ibd.dat2 <- harmonise_data(
  exposure_dat = ad_exp_dat, 
  outcome_dat = ibd_out_dat
)
ibd.res2 <- mr(ibd.dat2)
print(ibd.res2)
ibd.p2 <- mr_scatter_plot(ibd.res2,ibd.dat2)
ibd.p2[[1]]

################### AD as exposure, CD as outcome ###################
cd_exp_dat <- read_exposure_data(
  filename = cd_file,
  sep = "\t",
  snp_col = "SNP",
  beta_col = "bzx",
  se_col = "bzx_se",
  effect_allele_col = "a1",
  other_allele_col = "a2",
  pval_col = "bzx_pval",
  samplesize_col = "bzx_n"
)
cd_exp_dat$exposure="AD"
ad_out_dat <- read_outcome_data(
  filename = cd_file,
  sep = "\t",
  snp_col = "SNP",
  beta_col = "bzy",
  se_col = "bzy_se",
  effect_allele_col = "a1",
  other_allele_col = "a2",
  pval_col = "bzy_pval",
  samplesize_col = "bzy_n"
)
ad_out_dat$outcome="CD"
cd.dat2 <- harmonise_data(
  exposure_dat = cd_exp_dat, 
  outcome_dat = ad_out_dat
)
cd.res2 <- mr(cd.dat2)
print(cd.res2)
cd.p2 <- mr_scatter_plot(cd.res2, cd.dat2)
cd.p2[[1]]

################### AD as exposure, UC as outcome ###################
ad_exp_dat <- read_exposure_data(
  filename = uc_file,
  sep = "\t",
  snp_col = "SNP",
  beta_col = "bzx",
  se_col = "bzx_se",
  effect_allele_col = "a1",
  other_allele_col = "a2",
  pval_col = "bzx_pval",
  samplesize_col = "bzx_n"
)
ad_exp_dat$exposure="AD"
uc_out_dat <- read_outcome_data(
  filename = uc_file,
  sep = "\t",
  snp_col = "SNP",
  beta_col = "bzy",
  se_col = "bzy_se",
  effect_allele_col = "a1",
  other_allele_col = "a2",
  pval_col = "bzy_pval",
  samplesize_col = "bzy_n"
)
uc_out_dat$outcome="UC"
uc.dat2 <- harmonise_data(
  exposure_dat = ad_exp_dat, 
  outcome_dat = uc_out_dat
)
uc.res2 <- mr(uc.dat2)
print(uc.res2)
uc.p2 <- mr_scatter_plot(uc.res2, uc.dat2)
uc.p2[[1]]

library(patchwork)
combined_plot <- ibd.p2[[1]] + cd.p2[[1]] + uc.p2[[1]] & theme(legend.position = "top")
combined_plot <- combined_plot + plot_layout(guides = "collect")
