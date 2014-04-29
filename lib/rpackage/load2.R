srcdir=Sys.getenv('MWAS_GWAS_DIR')
source(sprintf('%s/src/lib/rpackage/util.r',srcdir))
source(sprintf('%s/src/lib/rpackage/util.load.r',srcdir))
source(sprintf('%s/src/lib/rpackage/general-utilities.r',srcdir))
source(sprintf('%s/src/lib/rpackage/plink.r',srcdir))
source(sprintf('%s/src/lib/rpackage/MX.dataset.functions.r',srcdir))
source(sprintf('%s/src/lib/rpackage/MX.statistical.tests.r',srcdir))
#source('~/drive/enterotypes/scripts/package/clustering.r',srcdir))
library('vegan')

"get.MX.dataset" <- function(
  mb, # microbiome table (taxa, otus, functions, bdiv pc's, alpha diversities)
  gx, # other table (e.g. genetics)
  md # metadata
){
  
  # get common rows of gx and mb and md
  common.rows <- intersect(rownames(mb$x), rownames(gx$x))
  common.rows <- intersect(common.rows, rownames(md$x))
  if(length(common.rows) < 2) stop(sprintf('%d rows in common between microbiome, genotype, metadata.\n',length(common.rows)))

  mb$x <- mb$x[common.rows,]
  gx$x <- gx$x[common.rows,]
  md$x <- md$x[common.rows,]  
  
  md$x <- droplevels(md$x)
  res <- list(mb=mb, gx=gx, md=md)  
  class(res) <- 'MX.dataset'
  return(res)
}