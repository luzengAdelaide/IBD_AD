library(dplyr)
library(data.table)
library(stringr)

# trait: IBD/CD/UC

# Load Delanger risk SNPs
risk_snps <- read.table("FILEPATH/reported.delange.trait.txt", header=TRUE, sep="\t")

risk_snps <- risk_snps %>% 
  dplyr::rename(Chr_37 = Chr,
                Pos_GRCh37 = Position,
                SNP_ID = SNP 
  )

# extract necessary columns
risk_snps <- risk_snps[, c("SNP_ID", "Chr_37", "Pos_GRCh37")]

# load GWAS summary statistics
GWAS <- fread("FILEPATH/delange.ytsiy.gwas.txt", sep="\t", header=T)
df_gwas <- as.data.frame(GWAS)

# IBD samples size
df_gwas$N = rep("59957", nrow(df_gwas))
# CD samples size
#df_gwas$N = rep("40266", nrow(df_gwas))
# UC samples size
#df_gwas$N = rep("45975", nrow(df_gwas))

## select the useful colnames for GWAS summary
 df_gwas <- df_gwas %>%
   dplyr::select(Chr, Position, MarkerName, TestedAllele, OtherAllele, Beta, StdErr, P.value, N) %>%
   dplyr::rename(
     Chr_GRCh37 = Chr,
     Pos_GRCh37 = Position,
     SNP = MarkerName,
     P = P.value
   )

# output directory
outdir <- "FILEPATH/IBD_1MB"
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
  region <- df_gwas %>%
    filter(
      Chr_GRCh37 == risk_snp$Chr_37,
      start_pos_in_GRCh37 < Pos_GRCh37,
      Pos_GRCh37 < end_pos_in_GRCh37
    )
  
  # export the extracted stats
  outfile <- sprintf("GWAS_IBD_%s.txt", risk_snp$SNP_ID)
  outpath <- file.path(outdir, outfile)
  write.table(region, file = outpath, sep = "\t",
              row.names = FALSE, col.names = TRUE, quote=FALSE)
  
}

