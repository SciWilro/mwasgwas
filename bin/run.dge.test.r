basedir <- Sys.getenv('MWAS_GWAS_DIR')
source(sprintf('%s/lib/rpackage/load2.R',basedir))
library('optparse',quietly=TRUE)

# set up and parse command-line options
option_list <- list(
  make_option(c("-i", "--microbiome_table"), type="character", default=NULL,
              help="Path to tab-delimited input table, see input_type [default %default]",
              metavar="path"),
  make_option(c("-t", "--microbiome_type"), type="character", default='taxa',
              help="Input table type, can be 'taxa', 'functions', or 'otutable' [default %default]",
              metavar="type"),
  make_option(c("-m", "--metadata"), type="character", default=NULL,
              help="Path to metadata file [default %default]"),
  make_option(c("-c", "--covariates"), type="character", default=NULL,
              help="Use these covariates in tests. Comma-separated list of metadata column headers. [default %default]",
              metavar="covariates"),
  make_option(c("-x", "--x_table"), type="character", default=NULL,
              help="Path to 'x' data file (e.g. snps counts) [default %default]",
              metavar="path"),
  make_option(c("-a", "--x_annotations"), type="character", default=NULL,
              help="Path to 'x' annotations file (gene labels for snps) [default %default]",
              metavar="path"),
  make_option(c("-X", "--x_feature"), type="character", default=NULL,
              help="Test only this feature from x table [default %default]",
              metavar="feature"),
  make_option(c("-T", "--trended_disp"), action="store_true",
              help="Estimate per-sample trended dispersion (slow) [default %default]"),
  make_option(c("-n", "--norm_factor_method"), type="character",default='none',
              help="edgeR calcNormFactors method (none, RLE, or upperquartile) [default %default]"),
  make_option(c('-p','--make_parallel_commands'), action='store_true',
              help="Make parallel commands for testing one SNP at a time [default %default]"),
  make_option(c('--use_bsub'), action='store_true',
              help="Use bsub when making parallel commands [default %default]"),
  make_option(c("-o", "--outpath"), type="character", default=NULL,
              help="path for saving output [default %default]",
              metavar="path")  
)

args <- parse_args(OptionParser(option_list = option_list))
if(is.null(args$x_feature) & is.null(args$make_parallel_commands)){
  stop('x_feature or make_parallel_commands required')
} 

if(args$make_parallel_commands){
  gx <- load.genotype.data(args$x_table, annotations.fp=args$x_annotations)
  for(snp in colnames(gx$x)){
    bsub <- sprintf('bsub -o maketable-%s.lsf -q hour -R "rusage[mem=32]" "',snp)
    basecmd <- sprintf('Rscript $MWAS_GWAS_DIR/bin/run.dge.test.r -i %s -t %s -m %s -x %s -n %s ',
                       args$microbiome_table, args$microbiome_type, args$metadata,
                       args$x_table, args$norm_factor_method)
    if(!is.null(args$trended_disp)){
      basecmd <- paste(basecmd,' -T', sep='')
    }
    if(!is.null(args$x_annotations)){
      basecmd <- paste(basecmd,' -a ',args$x_annotations, sep='')
    }
    basecmd <- paste(basecmd,' -X ',snp, sep='')
    basecmd <- paste(basecmd,' -o ',args$outpath,'-',snp,sep='')
    cmd <- paste(bsub, basecmd, '"', sep='')
    cat(cmd,'\n',sep='')
  }
} else {
  
  mb <- load.microbiome.data(args$microbiome_table,type=args$microbiome_type,collapse.at=.95)
  gx <- load.genotype.data(args$x_table, min.maf=0.1, annotations.fp=args$x_annotations)
  md <- load.metadata(args$metadata)
  
  # create mxwas obj
  mxd <- get.MX.dataset(mb=mb, gx=gx, md=md)
  
  # clean up covariates
  # This is a hack
  # need to generalize
  covariate.names=c('Antibiotics','Immuno','Inflamed','Age','Gender',
                    'Biopsy_Location_General','Disease','Disease_Location',
                    'Years_since_diagnosis','Collection_Center')
  mm <- covariate.model.matrix(mxd$md$x[,covariate.names])
  mm <- mm[,!grepl('Disease_LocationLx',colnames(mm))]
  mm <- mm[,!grepl('Disease_LocationL3',colnames(mm))]
  
  # run tests 
  cat('Running tests for SNP > 0...\n')
  res0 <- exact.test.edgeR.covariates(x=mxd$mb$x, y=mxd$gx$x[,args$x_feature]  > 0, covariates=mm, verbose=TRUE)
  cat('Running tests for SNP == 2...\n')
  res2 <- exact.test.edgeR.covariates(x=mxd$mb$x, y=mxd$gx$x[,args$x_feature] == 2, covariates=mm, verbose=TRUE)
  
  # extract and collate results
  tt1 <- topTags(res0,n=NULL)$table
  tt2 <- topTags(res2,n=NULL)$table
  tt1 <- cbind(rownames(tt1),tt1); colnames(tt1)[1] <- 'mbID'
  tt2 <- cbind(rownames(tt2),tt2); colnames(tt2)[1] <- 'mbID'
  tt1 <- cbind(rep(args$x_feature,nrow(tt1)),tt1); colnames(tt1)[1] <- 'xID'
  tt2 <- cbind(rep(args$x_feature,nrow(tt2)),tt2); colnames(tt2)[1] <- 'xID'
  rownames(tt1) <- NULL
  rownames(tt2) <- NULL
  res <- rbind(tt1,tt2)
  
  # save results
  write.table(res,file=args$outpath,col.names=FALSE,row.names=FALSE,sep='\t',quote=FALSE)
}