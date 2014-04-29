"print.MX.dataset" <- function(mxd){
  # Summary
  cat(sprintf('Microbiome+X dataset with %d subjects.\n',nrow(mxd$md$x)))
  
  # Metadata
  cat(sprintf('Metadata: %d columns\n',ncol(mxd$md$x)))
  print(mxd$md$x[1:5,1:5])
  cat('...\n\n')

  # Microbiome
  cat(sprintf('Microbiome: %d features, type "%s"\n',ncol(mxd$mb$x),mxd$mb$description))
  print(mxd$mb$x[1:5,1:5])
  cat('...\n\n')

  # X/genotype
  cat(sprintf('X-ome: %d features, type %s\n',ncol(mxd$gx$x),mxd$mb$description))
  print(mxd$gx$x[1:5,1:5])
  cat('...\n\n')
}