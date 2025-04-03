library(dplyr)
library(data.table)
library(stringr)

# note for all AD/IBD/CD/UC risk variants, we have included the chromosome position for both GRCh37 and GRCh38
# microglia eQTL was calculated based on GRCh38, while the monocyte eQTL was calculated based on GRCh37
# both of them were downloaded from previously published studies

# Load AD risk SNPs
risk_snps <- read.csv("FILEPATH/GRCH38_AD_Wightman_RiskSNP.csv", header=TRUE)

risk_snps <- risk_snps %>% 
  dplyr::rename(Chr_37 = mcols.Chr_GRCh37,
                Pos_GRCh37 = mcols.bp_GRCh37,
                SNP_ID = mcols.SNP_ID 
  )

# extract necessary columns
risk_snps <- risk_snps[, c("SNP_ID", "Chr_37", "Pos_GRCh37")]

# output directory
outdir <- "FILEPATH/AD_1MB"
dir.create(outdir, recursive = T)

################################################################################
eqtl <- fread("FILEPATH/monocyte.blueprint.tsv")
bh_eqtl <- fread("FILEPATH/monocyte.blueprint.BH-signif.txt")

get_eqtls_in_window <- function(eqtls, window_chr, window_center, window_size){
  window_start <- window_center - window_size
  window_end <- window_center + window_size
  
  eqtls_in_window <- eqtls %>%
    filter(
      chr == window_chr,
      window_start < pos,
      pos < window_end
    )
  
  return(eqtls_in_window)
}

create_output_path <- function(risk_snp, egene) {
  risk_snp_dir <- file.path(outdir, risk_snp$SNP_ID)
  dir.create(risk_snp_dir, recursive = T)
  outfile <- str_c("eQTL_", egene, ".txt")
  outpath <- file.path(risk_snp_dir, outfile)
  
  return(outpath)
}

################################################################################
risk_snps$Position_37 <- risk_snps$Pos_GRCh37

for(i in 1:nrow(risk_snps)) {
  risk_snp <- risk_snps[i, ]
  ##
  risk_snp$Position_37 <-  as.numeric(risk_snp$Position_37)
  #risk_snp$Chr_38 <- as.numeric(risk_snp$Chr_38)
  ###
  # get eQTL around the risk SNP
  bh_eqtls_around_risk_snp <- get_eqtls_in_window(bh_eqtl,
                                                  risk_snp$Chr_37,
                                                  risk_snp$Position_37,
                                                  5e+5)
  
  # skip the risk SNP if it has no eQTL around it
  if(nrow(bh_eqtls_around_risk_snp) == 0) {
    sprintf("No eQTL was found around %s\n", risk_snp$SNP_ID) %>% cat
    next 
  }
  
  # list all eGenes whose eSNPs are around the risk SNP
  egenes <- unique(bh_eqtls_around_risk_snp$gene)
  
  for(egene in egenes) {
    
    # get the lead eSNP
    lead_esnp <- bh_eqtls_around_risk_snp %>%
      filter(gene == egene) %>%
      slice_min(pvalue)      
    
    # get all MatrixEQTL output around the lead eSNP
    eqtls_around_lead_esnp <- eqtl %>%
      filter(gene == egene) %>% 
      get_eqtls_in_window(., lead_esnp$chr, lead_esnp$pos, 1e+6)
    
    # export
    outpath <- create_output_path(risk_snp, egene)
    write.table(eqtls_around_lead_esnp, file = outpath, sep = "\t",
                row.names = FALSE, col.names = TRUE, quote=FALSE)
  }
  
}

