library(dplyr)
library(data.table)
library(stringr)

# Load RiskSNPs in GRCH38
risk_snps <- read.csv("FILEPATH/GRCH38_AD_Wightman_RiskSNP.csv", header=TRUE)

risk_snps <- risk_snps %>% 
  dplyr::rename(Chr_37 = mcols.Chr_GRCh37,
                Pos_GRCh37 = mcols.bp_GRCh37,
                SNP_ID = mcols.SNP_ID 
  )

# extract necessary columns
risk_snps <- risk_snps[, c("SNP_ID", "Chr_37", "Pos_GRCh37")]

# load GWAS summary statistics
GWAS <- read.table("FILEPATH/wightman2021.gwas", sep="\t", header=T)
GWAS <- as.data.frame(GWAS)

# output directory
outdir <- "FILEPATH/mic_1MB/AD_1MB"
dir.create(outdir, recursive = T)

################################################################################
## extract GWAS summary stats for each risk SNP in 1MB
window_size <- 1e+6

risk_snps$Position_37 <- risk_snps$Pos_GRCh37

for(i in 1:nrow(risk_snps)){
  risk_snp <- risk_snps[i, ]
  
  risk_snp$Position_37 <-  as.numeric(risk_snp$Position_37)
  risk_snp$Chr_37 <- as.numeric(risk_snp$Chr_37)
  
  start_pos_in_GRCh37 <- risk_snp$Position_37 - window_size
  end_pos_in_GRCh37 <- risk_snp$Position_37 + window_size
  
  # extract GWAS summary stats around the risk SNP
  region <- GWAS %>%
    filter(
      Chr == risk_snp$Chr_37,
      start_pos_in_GRCh37 < Pos,
      Pos < end_pos_in_GRCh37
    )
  
  region <- region %>%
    group_by(SNP) %>%
    filter(p == min(p)) %>%
    ungroup()
  
  # export the extracted stats
  outfile <- sprintf("GWAS_AD_%s.txt", risk_snp$SNP_ID)
  outpath <- file.path(outdir, outfile)
  write.table(region, file = outpath, sep = "\t",
              row.names = FALSE, col.names = TRUE, quote=FALSE)
  
}

