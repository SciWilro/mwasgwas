# usage
# Rscript deduplicate_plink_data.r -m plink.map -f plink.frq -o duplicate_snps_list.txt
library('optparse',quietly=TRUE)

# set up and parse command-line options
option_list <- list(
    make_option(c("-f", "--freq_table"), type="character", default=NULL,
        help="Plink frequency table (output from plink --freq) [default %default]",
        metavar="path"),
    make_option(c("-m", "--plink_map"), type="character", default=NULL,
        help="Path to plink map file [default %default]",
        metavar="path"),
    make_option(c("-o", "--outfile"), type="character", default=NULL,
        help="path for output list of snps to drop [default %default]",
        metavar="path"),
    make_option(c("-a", "--include_alleles"), action="store_true", default=FALSE,
        help="include alleles in uids [default %default]",
        metavar="path")  

)
args <- parse_args(OptionParser(option_list = option_list))

# load plink map
mp <- read.table(args$plink_map,sep='\t',head=F)

# load plink frequency table
pfreq <- read.table(args$freq_table,head=T,check=F)

# get overlapping snps
ix <- intersect(mp[,2], pfreq[,'SNP'])
pfreq <- pfreq[match(ix,pfreq[,'SNP']),]
mp <- mp[match(ix,mp[,2]),]

# create new snp IDs using chr, pos, minor allele, major allele
mp.ids <- sprintf('uid_%d_%d',mp[,1],mp[,4])

if(args$include_alleles){
    alleles <- cbind(pfreq[,'A1'],pfreq[,'A2'])
    for(i in 1:nrow(alleles)){alleles.i <- sort(alleles[i,]); mp.ids[i] <- sprintf('%s_%s',mp.ids[i], paste(alleles.i,collapse='_'))}
}

# check for duplicates
# dups <- which(table(mp.ids) > 1)
# if(length(dups)>0){
#     stop('There are duplicate snp IDs.\n')
# }

# write output
sink(args$outfile)
cat(sprintf('%s\t%s\n',mp[,2], mp.ids),sep='')
sink(NULL)