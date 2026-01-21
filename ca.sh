#!/bin/bash
#SBATCH -A hugen2072-2026s
#SBATCH --mem-per-cpu=60G 
#SBATCH --output=Log_%x.out
#SBATCH --time=3:00:00 
#SBATCH --mail-user=tsb31@pitt.edu
#SBATCH --mail-type=BEGIN,END,FAIL,TIMEOUT
#SBATCH --cpus-per-task=2

# load our required tools
module load bcftools/1.22
module load vcftools
module load plink/1.90b7.7
# used the officia PLINK docs to help generate "olink" commands below
# https://www.cog-genomics.org/plink/1.9/

# source our convenience variables
source vars.sh

# 1. Create and run a script on the CRC cluster that 
# Oh would you look at that...we're in said script...

# move down into our working directory
# (and make the directory if it doesn't exist)
mkdir -p $WD
cd $WD

# 2. Copies data.vcf to your working directory
cp $DATA/data.vcf .

# 4. Is filtered to include only positions on chromosome 4; 
# NOTE: running this before part 3...
# using this reference from biostars: https://www.biostars.org/p/201603/
# went with this methid because I wanted to retain the headers
vcftools --vcf data.vcf --out chr4 --chr 4 --recode

# 3. Then uses bcftools to create a data4.bcf.gz file that is sorted and indexed, and
# create and sort the file
bcftools sort chr4.recode.vcf -Ob -o data4.bcf.gz 
# index the file
bcftools index data4.bcf.gz

# 5. Then uses PLINK to create a PLINK binary file set version of data4.bcf.gz called data4.{fam,bim,bed};
# --AND--
# 6. Then uses PLINK to update sex variable in the data4.fam file using the sex variable in sex.txt 
# (without copying sex.txt to your own directory),
plink --bcf data4.bcf.gz --update-sex $SEX --make-bed --out data4

# 7. Then uses PLINK, coupled with the phenotypes in phenotype.txt 
# (without copying phenotype.txt to your own directory—and you should use phenotype.txt as an auxiliary file 
# for the next few tasks without altering data4.fam to include phenotype data), to
# NOTE: going to save this to a var in vars.sh that is sourced for this script.

# can calculate the allele frequencies of both and then filer to individual files
# used this resource for the correct way to use --freq in this case: https://www.cog-genomics.org/plink/1.9/basic_stats#freq
# used this resource to know how to split to individual case and control files: https://www.cog-genomics.org/plink/1.9/formats#frq
# NOTE: I am just running the combined case-control option given by plink and manually splitting the file
plink --bfile data4 --pheno $PHE --freq case-control --out data4
# create a variable to load the heading for the freq files to be made
HEADERS="CHR SNP A1 A2 MAF NCHROBS"

# 8. Calculate the allele frequencies of the markers in data4.{fam,bim,bed} in the cases only and 
# calculated above, going to filter to ONLY cases
CTRLS="data4.ctrls.freq"
echo $HEADERS > $CTRLS
tail -n+2 data4.frq.cc | awk '{ print $1, $2, $3, $4, $6, $8 }' >> $CTRLS

# 9. Calculate the allele frequencies of the markers in data4.{fam,bim,bed} in the controls only; then 
# calculated above, going to filter to ONLY controls
CASES="data4.cases.freq"
echo $HEADERS > $CASES
tail -n+2 data4.frq.cc | awk '{ print $1, $2, $3, $4, $5, $7 }' >> $CASES

# 10. Performs a GWAS of the phenotype using logistic regression with no covariates. 
# NOTE: used this resource as a reminder (because Ryan was right and I forgot how to do a GWAS XD):
# https://www.cog-genomics.org/plink/2.0/tutorials/qwas1
# which also led me to this github repo to reference
# https://www.cog-genomics.org/plink/2.0/tutorials/qwas1
plink --vcf data.vcf --update-sex $SEX --pheno $PHE --logistic --out data

# Some clean up because I am going to use RStudio on my laptop instead of on the CRC
mkdir -p ../data
cp $CTRLS ../data
cp $CASES ../data
cp data.assoc.logistic ../data

