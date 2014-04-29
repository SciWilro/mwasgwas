# usage:
# flip_strands inbasename

# get freq table
plink --noweb --file ${1} --freq --out ${1}
# get list of snps to flip
grep " C [ ]*T \|  G [ ]*T \|  T [ ]*C \|  T [ ]*G " ${1}.frq | sed -e 's: [ ]*: :g' | sed -e 's:^ ::' | cut -f 2 -d ' ' > ${1}_fliplist.txt
# flip
plink --noweb --file ${1} --flip ${1}_fliplist.txt --recode --out ${1}-flip2
