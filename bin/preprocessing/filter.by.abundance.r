# filters OTU/taxon/pathway table by abundance
# usage:
# filter.by.abundance.py -i taxa.txt -m .001 -p .01 -o taxa-filtered.txt

source(sprintf('%s/src/lib/rpackage/util.load.r',Sys.getenv('MWAS_GWAS_DIR')))
source(sprintf('%s/src/lib/rpackage/general-utilities.r',Sys.getenv('MWAS_GWAS_DIR')))
source(sprintf('%s/src/lib/rpackage/util.r',Sys.getenv('MWAS_GWAS_DIR')))
library('optparse',quietly=TRUE)

# set up and parse command-line options
option_list <- list(
    make_option(c("-i", "--input_table"), type="character", default=NULL,
        help="Path to tab-delimited input table, see input_type [default %default]",
        metavar="path"),
    make_option(c("-t", "--input_type"), type="character", default='taxa',
        help="Input table type, can be 'taxa' or 'functions' or 'otutable' [default %default]",
        metavar="type"),
    make_option(c("-m", "--min.fraction"), type="numeric", default=0.001,
        help="Minimum mean fraction of communities composed of this feature [default %default]",
        metavar="path"),
    make_option(c("-p", "--min.prevalence"), type="numeric", default=0.01,
        help="Minimum fraction of communities in which feature is present [default %default]",
        metavar="path"),
    make_option(c("-o", "--outpath"), type="character", default=NULL,
        help=" [default %default]",
        metavar="path")  
)
args <- parse_args(OptionParser(option_list = option_list))

if(args$input_type == 'taxa'){
    x <- load.qiime.taxon.table(args$input_table)
} else if(args$input_type == 'otutable'){
    x <- load.qiime.otu.table(args$input_table)
} else {
    x <- load.qiime.mapping.file(args$input_table)
}

nx <- sweep(x,1,rowSums(x),'/')
x <- x[,colMeans(nx) > args$min.fraction] 
x <- x[,colMeans(x > 0) > args$min.prevalence]

sink(args$outpath)
if(args$input_type == 'taxa'){
    cat('Taxon\t')
    write.table(t(x),sep='\t',quote=F)
} else if(args$input_type == 'otutable'){
    cat('# Obligatory comment line\n#OTU ID\t')
    write.table(t(x),sep='\t',quote=F)
} else {
    cat('Function\t')
    write.table(x,sep='\t',quote=F)
}
sink(NULL)