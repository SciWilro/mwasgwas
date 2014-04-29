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
        metavar="path")  
)
args <- parse_args(OptionParser(option_list = option_list))

# load plink map
mp <- read.table(args$plink_map,sep='\t',head=F)

# load plink frequency table
pfreq <- read.table(args$freq_table,head=T,check=F)

# create new snp IDs using chr, pos, minor allele, major allele
mp.ids <- sprintf('uid_%d_%d',mp[,1],mp[,4])
# alleles <- cbind(pfreq[,'A1'],pfreq[,'A2'])
# for(i in 1:nrow(alleles)){alleles.i <- sort(alleles[i,]); mp.ids[i] <- sprintf('%s_%s',mp.ids[i], paste(alleles.i,collapse='_'))}

# get list of duplicates
dups <- which(table(mp.ids) > 1)
if(length(dups) > 1){
    bad.snps <- NULL
    for(i in 1:length(dups)){
        snp.names <- pfreq[mp.ids==names(dups)[i],'SNP']
        mafs <- pfreq[mp.ids==names(dups)[i],'MAF']
        
        # all NA - remove all
        if(all(is.na(mafs))){
            bad.snps <- c(bad.snps, snp.names)
        } else {
            # keep the snp with highest MAF
            mafs[is.na(mafs)] <- -1
            bad.snps <- c(bad.snps, snp.names[-which.max(mafs)])
        }
    }
    sink(args$outfile)
    cat(bad.snps,sep='\n')
    sink(NULL)
} else {
    sink(args$outfile)
    cat('XYZXYZXYZXYZ\n')
    sink(NULL)
}
