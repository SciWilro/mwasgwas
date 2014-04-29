# run with sh extract-and-rename-snps-from-ped.sh plink-base prefix
# e.g. sh extract-and-rename-snps-from-ped.sh ~/Desktop/archive/IBDPRISM_newdata/subset prism

SNPTABLE="$HOME/drive/research/prism/data/genetic/2013-01-09-generate-snps/all-snps-metadata-with-pos.txt"
RISKALLELES="$HOME/drive/research/prism/data/genetic/2013-01-09-generate-snps/risk_alleles.txt"
PLINKBASE=$1
PREFIX=$2

# for all ichip variants, find all variants within 10 bases using dbsnp data
# this file will be called ichipmap-full.txt.
# use this file to add coords to snp metadata table

# using snp metadata and an R script and plink map,
# extract all identical matches from plink map first by rsID, then by coords
# output is two columns
# ichipID rsID
Rscript ../../../code/extract_best_snp_matches_from_plink_map.r ${SNPTABLE} ${PLINKBASE}.map > ic-rs-mapping-${PREFIX}.txt

# use this to 
# 1. extract subset from plink/ped file
cut -f 1 ic-rs-mapping-${PREFIX}.txt > ic-id-subset-${PREFIX}.txt
plink --file ${PLINKBASE} --extract ic-id-subset-${PREFIX}.txt --out subset-${PREFIX} --tab --recode
# dedup subjects
sort -u -k 2,2 subset-${PREFIX}.ped > tmp; mv tmp subset-${PREFIX}.ped

# 2. rename snps in map file with identical rsID matches
time plink --file subset-${PREFIX} --update-map ic-rs-mapping-${PREFIX}.txt --update-name --out subset-${PREFIX}-rename --tab --recode

# convert to raw format using risk allele table
time plink --file subset-${PREFIX}-rename --recodeA --recode-allele ${RISKALLELES} --tab --out subset-${PREFIX}-rename

