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
        help="Input table type, can be 'metadata' or 'taxa' or 'otutable' [default %default]",
        metavar="type"),
    make_option(c("-m", "--metadata"), type="character", default=NULL,
        help="Path to metadata file [default %default]",
        metavar="path"),
    make_option(c("-c", "--category"), type="character", default=NULL,
        help="Collapse according to this category/column in the metadata [default %default]",
        metavar="category"),
    make_option(c("-s", "--states"), type="character", default=NULL,
        help="Only keep samples with these values in metadata, e.g. 'GENDER:Male;COUNTRY:*,!USA'; if NULL, includes only samples in metadata file [default %default]",
        metavar="state_list"),
    make_option(c("-f", "--collapse_function"), type="character", default='mean',
        help="Function with which to collapse ('mean', 'sum', 'mode', 'cat'). 'mode' automatically used for non-numeric data. [default %default]",
        metavar="function"),
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

metadata <- load.qiime.mapping.file(args$metadata)
# filter metadata if requested
if(is.null(args$states)) {
    sample.ids <- rownames(metadata)
} else {
    sample.ids <- extract.samples.by.metadata(metadata, args$states)
    metadata <- droplevels(metadata[sample.ids,,drop=F])
}

common <- intersect(rownames(x), rownames(metadata))
cat(sprintf('%d of %d samples retained from data file\n',length(common), nrow(x)))
metadata <- droplevels(metadata[common,,drop=F])
x <- x[common,,drop=F]
category <- metadata[[args$category]]
cat(sprintf('%d samples after collapsing\n',length(unique(category))))


# do the collapse
# f is the collapsing function
f <- function(xx) {
    if(any(is.na(xx))){
        return(NA)
    }
    if(class(xx)=='factor' || class(xx)=='character' || args$collapse_function=='mode'){
        if(args$collapse_function=='cat'){
            # concatenate all types
            return(paste(unique(xx),collapse='__'))
        } else {
            # return the mode
            return(names(sort(-table(xx)))[1])
        }
    } else if(args$collapse_function == 'mean'){
        return(mean(xx))
    } else {
        return(sum(xx))
    }
}

xnew <- sapply(1:ncol(x), function(ix) {sapply(split(x[,ix],category),f)})
colnames(xnew) <- colnames(x)
# xnew <- apply(x, 2, function(xx) sapply(split(xx,category),f))
sink(args$outpath)
if(args$input_type == 'taxa'){
    cat('Taxon\t')
    write.table(t(xnew),sep='\t',quote=F)
} else if(args$input_type == 'otutable'){
    cat('# Obligatory comment line\n#OTU ID\t')
    write.table(t(xnew),sep='\t',quote=F)
} else {
    cat('#SampleID\t')
    write.table(xnew,sep='\t',quote=F)
}
sink(NULL)