# calculate pairwise LD scores between AD and IBD/CD/UC

# trait: IBD/CD/UC

# extract snps based on chromosome and convert it to bed format
awk '{print $1"\t"$2"\t"$2+1"\t"$3}' ad_snps.txt > ad_snps.bed
awk '{print $1"\t"$2"\t"$2+1"\t"$3}' ibd_snps.txt > ibd_snps.bed
awk '{print $1"\t"$2"\t"$2+1"\t"$3}' cd_snps.txt > cd_snps.bed
awk '{print $1"\t"$2"\t"$2+1"\t"$3}' uc_snps.txt > uc_snps.bed

# separate bed files to chromosome 1-22
for trait in *bed; do for chr in {1..22}; do grep -w $chr $trait > $chr.$trait; done; done

# intersect SNPs within a 500kb window of upstream and downstream of AD variants:
module load Bedtools/2.31 
for i in {1..22}; do bedtools window -a $i.ad_snps.bed -b $i.$trait\_snps.bed -w 5000000 > overlap.$i.ad.$trait.txt; done

# remove empty file
find . -type f -empty -delete

# extract overlapped SNPs to test
cat overlap.*$trait* |sort -nk1,1 -nk2,2 |cut -f1,4,8 > test.ad.$trait.txt

# calculate ld scores
module load Plink/1.9.10
while read -r value1 value2 value3 remainder; 
do
    echo "$value1 $value2 $value3";   
    plink --bfile ~/1000_genome/eur/eur.$value1 --ld $value2 $value3
done < test.ad.$trait.txt
