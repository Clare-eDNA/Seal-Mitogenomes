module load  VCFtools/0.1.14-gimkl-2018b-Perl-5.28.1

vcftools --vcf furseal.vcf --min-alleles 2 --max-alleles 2 --remove-indels --recode-INFO-all --recode --out ./FilteredSNPs/SNPset1

vcftools --vcf ./FilteredSNPs/SNPset1.recode.vcf --missing-indv

awk '!/IN/' out.imiss | cut -f5 > totalmissing
cat out.imiss

gnuplot << \EOF
set terminal dumb size 120, 30
set autoscale 
unset label
set title "Histogram of % missing data per individual"
set ylabel "Number of Occurrences"
set xlabel "% of missing data"
set yr [0:100000]
binwidth=0.01
bin(x,width)=width*floor(x/width) + binwidth/2.0
plot 'totalmissing' using (bin($1,binwidth)):(1.0) smooth freq with boxes
pause -1
EOF

awk '$5 > 0.40' out.imiss | cut -f1 > lowDP.indv

vcftools --vcf ./FilteredSNPs/SNPset1.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out ./FilteredSNPs/RemovedIndivid

vcftools --vcf ./FilteredSNPs/RemovedIndivid.recode.vcf --hardy --max-missing 1.0 --out HWEs


vcftools --vcf ./FilteredSNPs/RemovedIndivid.recode.vcf  --max-missing 0.98 --maf 0.05 --min-meanDP 3 --hwe 0.05 --remove-filtered-all --recode-INFO-all --recode --out ./FilteredSNPs/Test1


vcftools --vcf ./FilteredSNPs/Test1.recode.vcf --plink


#you should have .ped and .map files now

module load  R/3.6.2-gimkl-2020a


# See file "How I filtered SNPs for NZ Fur Seals" in documents in NZ Fur Seals folder
