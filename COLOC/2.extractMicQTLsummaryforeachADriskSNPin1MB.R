library(dplyr)
library(data.table)
library(stringr)

# eQTL directory of the microglia
eqtl_dir <- "FILEPATH/"
# output directory
outdir <- "FILEPATH/mic_1MB/AD_1MB"
dir.create(outdir, recursive = T)

# read risk snps in GRCH38 coordinates
risk_snps <- read.csv("FILEPATH/GRCH38_AD_Wightman_RiskSNP.csv", header=TRUE)

risk_snps <- risk_snps %>% 
  dplyr::rename(Chr_38 = Chr_GRCh38,
                Pos_GRCh38 = Pos_GRCh38,
                SNP_ID = mcols.SNP_ID 
  )


risk_snps <- risk_snps %>% mutate(Chr_38 = str_c("chr", Chr_38))  # do not forget "chr"

# extract necessary columns
risk_snps <- risk_snps[, c("SNP_ID", "Chr_38", "Pos_GRCh38")]

################################################################################

# load SNP positions
snppos <- fread("/mnt/mfs/ctcn/team/masashi/snuc-eqtl/genotype/get-dosage.ALL.snppos")

# load eQTL and BH significant eQTLs
eqtl_file <- file.path(eqtl_dir, "matrix-eqtl.rds")
eqtl <- readRDS(eqtl_file)$cis$eqtls

bh_eqtl_file <- file.path(eqtl_dir, "matrix-eqtl.BH-signif.rds")
bh_eqtl <- readRDS(bh_eqtl_file)

get_eqtls_in_window <- function(eqtls, window_chr, window_center, window_size){
  window_start <- window_center - window_size
  window_end <- window_center + window_size
  
  snps_in_window <- snppos %>%
    filter(
      chr == window_chr,
      window_start < pos,
      pos < window_end
    )
  
  eqtls_in_window <- eqtls %>%
    inner_join(snps_in_window, by = c("snps" = "snpid"))
  
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
risk_snps$Position_38 <- risk_snps$Pos_GRCh38

for(i in 1:nrow(risk_snps)) {
  risk_snp <- risk_snps[i, ]
  ##
  risk_snp$Position_38 <-  as.numeric(risk_snp$Position_38)
  #risk_snp$Chr_38 <- as.numeric(risk_snp$Chr_38)
  ###
  # get eQTL around the risk SNP
  bh_eqtls_around_risk_snp <- get_eqtls_in_window(bh_eqtl,
                                                  risk_snp$Chr_38,
                                                  risk_snp$Position_38,
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
