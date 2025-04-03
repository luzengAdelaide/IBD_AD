library(dplyr)
library(stringr)
library(coloc)
library(data.table)

# trait: AD/IBD/CD/UC

# load GWAS summary stats of each risk SNP
gwas_window_dir <-  "FILEPATH/mic_1MB/trait_1MB/"

# Load eQTL for each risk snp for each cell type
eqtl_window_dir <-  "FILEPATH/mic_1MB/trait_1MB/"
pattern <- file.path(eqtl_window_dir, "rs*", "eQTL_*.txt")
eqtl_windows <- Sys.glob(pattern)

# apply COLOC function for eqtl and gwas, can use MAF or varbeta depends what provide from GWAS 
for(eqtl_window_file in eqtl_windows) {
  rsid = eqtl_window_file %>% dirname %>% basename
  eqtl_window <- fread(eqtl_window_file)
  
  gwas_window_file <- file.path(gwas_window_dir, str_c("GWAS_trait_", rsid, ".txt"))
  gwas_window <- fread (gwas_window_file)
  
  df <- gwas_window %>%
    inner_join(eqtl_window, by = c("SNP" = "snps"))
  
  df$varbeta <- (df$beta / df$statistic)^2
  df$varbeta2 <- (df$se)^2
  df <- df %>% filter(!is.na(varbeta))
  df <- df %>% filter(!is.na(b))
  dataset2 <- list(beta=df$beta, varbeta=df$varbeta, sdY=1, snp=df$SNP, type="quant")
  
  dataset1 <- list(beta=df$b, varbeta=df$varbeta2, snp=df$SNP, type="cc", s=0.08018904)
  
  res <- coloc.abf(dataset1, dataset2)
  
  # regular expression
  gene <- eqtl_window_file %>% basename %>% str_match("eQTL_(.+)\\.txt") %>% .[2]
  
  dir = "FILEPATH/mic_1MB/trait_1MB"

  outfile <- sprintf("%s/%s-%s.coloc.tsv", dir, rsid, gene)
  
  summary.df <- data.frame(
    key = names(res$summary),
    value = res$summary,
    row.names = NULL) %>%
    mutate(
      snp = rsid,
      gene = gene
    )
  
  fwrite(summary.df, file=outfile, sep = "\t")
}
