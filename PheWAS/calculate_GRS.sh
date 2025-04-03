# calculate genetic risk score of IBD/CD/UC
module load Plink/1.9.10

plink \
    --vcf ibd.rosmap.vcf  \
    --score ../gwas/reported.delange.ibd_ibd.txt 1 4 6 header \
    --out risk.ibd

plink \
    --vcf cd.rosmap.vcf  \
    --score ../gwas/reported.delange.ibd_cd.txt 1 4 6 header \
    --out risk.cd

plink \
    --vcf uc.rosmap.vcf  \
    --score ../gwas/reported.delange.ibd_uc.txt 1 4 6 header \
    --out risk.uc

sed -i 's/ \+/\t/g; s/^\t//; s/FID/projid/g' risk.*.profile 
