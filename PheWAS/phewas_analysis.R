# PheWAS analysis
library(PheWAS)
library(reshape2)
library(ggpubr)
library(ggrepel)

setwd("FILEPATH/")
# read phenotype data
pheno = read.table("FILEPATH/phenos_AD_ROSMAP.txt", sep="\t", header=T)
# 1 is male, 0 is female
pheno$msex = gsub('0', 'F',
                  gsub('1', 'M', pheno$msex))

covariate = pheno[,c(1:4)]
use.pheno = pheno[,-c(2:4)]

# read genotype data
geno = read.csv("FILEPATH/risk.uc.profile", sep="\t", header=T)
use.geno = geno[,c(1,6)]

joind.data = inner_join(use.pheno,use.geno, by="projid")
joind.data2 = inner_join(joind.data,covariate, by="projid")

results=phewas_ext(phenotypes=names(use.pheno)[-1],
                   genotypes=names(use.geno)[-1],
                   covariates = names(covariate)[-1],
                   data=joind.data2,cores=4)

results$fdr = p.adjust(results$p)

cbPalette <- c( "purple", "navyblue", "#2077B4", "darkgreen", "lightgreen", "red", "hotpink","pink", "plum", "skyblue")

p1=ggplot(results, aes(x=phenotype, y=-log10(p))) + 
  geom_point(aes(col=phenotype, size=beta)) + 
  theme_classic() + 
  scale_colour_manual(values=cbPalette) +
  theme(axis.text.x = element_blank(), 
        panel.grid.minor=element_line(colour = "grey", linetype="dashed"), axis.ticks=element_blank()) + 
  labs(color="Category", size="Effect size", x="AD endophenotypes", y="-log10(FDR)") +
  theme(text = element_text(size=12)) +
  geom_text_repel(data=. %>% mutate(label = ifelse(p < 0.05, as.character(phenotype), "")), aes(label=label), size=4, box.padding = unit(0.7, "lines")) + 
  geom_hline(yintercept=-log10(0.05), color="red", size=1, alpha=0.5)

### CD
cd.geno = read.csv("FILEPATH/risk.cd.profile", sep="\t", header=T)
use.cd.geno = cd.geno[,c(1,2)]

cd.joind.data = inner_join(use.pheno,cd.geno, by="projid")
cd.joind.data2 = inner_join(cd.joind.data,covariate, by="projid")

cd.results=phewas_ext(phenotypes=names(use.pheno)[-1],
                   genotypes=names(use.cd.geno)[-1],
                   covariates = names(covariate)[-1],
                   data=cd.joind.data2,cores=4)

cd.results$fdr = p.adjust(cd.results$p)

cbPalette <- c( "purple", "navyblue", "#2077B4", "darkgreen", "lightgreen", "red", "hotpink","pink", "plum", "skyblue")

p2 = ggplot(cd.results, aes(x=phenotype, y=-log10(p))) + 
  geom_point(aes(col=phenotype, size=beta)) + 
  theme_classic() + 
  scale_colour_manual(values=cbPalette) +
  theme(axis.text.x = element_blank(), 
        panel.grid.minor=element_line(colour = "grey", linetype="dashed"), axis.ticks=element_blank()) + 
  labs(color="Category", size="Effect size", x="AD endophenotypes", y="-log10(p-value)") +
  theme(text = element_text(size=12)) +
  geom_text_repel(data=. %>% mutate(label = ifelse(p < 0.05, as.character(phenotype), "")), aes(label=label), size=4, box.padding = unit(0.7, "lines")) + 
  geom_hline(yintercept=-log10(0.05), color="red", size=1, alpha=0.5)

### IBD
ibd.geno = read.csv("FILEPATH/risk.ibd.profile", sep="\t", header=T)
use.ibd.geno = ibd.geno[,c(1,6)]

ibd.joind.data = inner_join(use.pheno,ibd.geno, by="projid")
ibd.joind.data2 = inner_join(ibd.joind.data,covariate, by="projid")

ibd.results=phewas_ext(phenotypes=names(use.pheno)[-1],
                      genotypes=names(use.ibd.geno)[-1],
                      covariates = names(covariate)[-1],
                      data=ibd.joind.data2,cores=4)

ibd.results$fdr = p.adjust(ibd.results$p)

cbPalette <- c( "purple", "navyblue", "#2077B4", "darkgreen", "lightgreen", "red", "hotpink","pink", "plum", "skyblue")

p3 = ggplot(ibd.results, aes(x=phenotype, y=-log10(p))) + 
  geom_point(aes(col=phenotype, size=beta)) + 
  theme_classic() + 
  scale_colour_manual(values=cbPalette) +
  theme(axis.text.x = element_blank(), 
        panel.grid.minor=element_line(colour = "grey", linetype="dashed"), axis.ticks=element_blank()) + 
  labs(color="Category", size="Effect size", x="AD endophenotypes", y="-log10(p-value)") +
  theme(text = element_text(size=12)) +
  geom_text_repel(data=. %>% mutate(label = ifelse(p < 0.05, as.character(phenotype), "")), aes(label=label), size=4, box.padding = unit(0.7, "lines")) + 
  geom_hline(yintercept=-log10(0.05), color="red", size=1, alpha=0.5)

combined_plot <- p3 + p2 + p1
print(combined_plot)
