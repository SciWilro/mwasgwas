mbdir <- Sys.getenv('MWAS_DIR')
source(sprintf('%s/lib/wrap.edgeR.r',mbdir))

# tests every microbiome feature (taxon/otu/function) 
# against a single binary or numeric outcome
# uses a linear test, so most appropriate for numeric outcome
#
# returns a matrix of taxon, coef, pval, fdr for test
#
# mb is the microbiome table (taxa/otus/functions/pcoa/alpha div)
# x is the outcome being tested
# drop.outliers.range is the number of quartiles above 3rd
"association.tests" <- function(mb, x, covariates=NULL,
                                verbose=FALSE,
                                drop.outliers.range=NULL
){
  # convenience: reshape, provide names, etc.
  if(is.null(dim(mb))) mb <- matrix(mb, ncol=1)
  if(!is.null(covariates) & is.null(dim(covariates))) covariates <- data.frame(covariate1=covariates)    
  if(is.null(colnames(mb))) colnames(mb) <- sprintf('col%d',1:ncol(mb))
  
  # validate inputs
  association.tests.validation(mb, x, covariates)
  
  N <- ncol(mb)
  res <- data.frame(mbID=colnames(mb),
                    direction=numeric(N),
                    statistic=numeric(N),
                    stringsAsFactors=FALSE,
                    pvalue=numeric(N),
                    qvalue=numeric(N))
  res[,c('pvalue','qvalue','direction','statistic')] <- NA
  
  for(i in 1:N){
    if(verbose & i %% 1 == 0) cat('',i,'')
    # drop outliers
    ix <- !is.na(x)
    if(!is.null(drop.outliers.range)) ix <- ix & !is.outlier(mb[,i],range=drop.outliers.range) 
    
    # precompute model matrix
    mm <- covariates[ix,,drop=F]
    
    if(is.null(covariates)){
      mlm <- lm(mb[ix,i] ~ x[ix])
    } else {
      mlm <- lm(mb[ix,i] ~ x[ix] + mm)
    }
    res$pvalue[i] <- summary(mlm)[[4]][2,'Pr(>|t|)']
    res$direction[i] <- mlm$coefficients[2]
    res$statistic[i] <- summary(mlm)[[4]][2,'t value']
  }		
  
  qvals <- p.adjust(drop.na(res[,'pvalue']),method='fdr')
  res[,'qvalue'] <- qvals
  
  res <- res[order(res[,'pvalue']),]
  res <- res[order(res[,'qvalue']),]
  
  return(res)
}

# tests every microbiome feature (taxon/otu/function) 
# against every column in matrix X (e.g. genotype allele counts)
# uses a linear test, so most appropriate for numeric outcome
#
# returns a matrix of taxon, coef, pval, fdr for test
#
# mb is the microbiome table (taxa/otus/functions/pcoa/alpha div)
# x is the outcome being tested
# drop.outliers.range is the number of quartiles above 3rd
"batch.association.tests" <- function(mb, X, covariates=NULL,
                                verbose=TRUE,
                                drop.outliers.range=NULL
){
  # convenience: reshape, provide names, etc.
  if(is.null(dim(mb))) mb <- matrix(mb, ncol=1)
  if(is.null(dim(X))) X <- matrix(X, ncol=1)
  if(!is.null(covariates) & is.null(dim(covariates))) covariates <- data.frame(covariate1=covariates)    
  if(is.null(colnames(mb))) colnames(mb) <- sprintf('col%d',1:ncol(mb))
  if(is.null(colnames(X))) colnames(X) <- sprintf('col%d',1:ncol(X))
  
  N <- ncol(mb)
  M <- ncol(X)
  res <- data.frame(mbID=rep(colnames(mb),M),
                    xID=rep(colnames(X),each=N),
                    direction=numeric(N*M),
                    statistic=numeric(N*M),
                    pvalue=numeric(N*M),
                    qvalue=numeric(N*M),
                    stringsAsFactors=FALSE)
  res[,c('pvalue','qvalue','direction','statistic')] <- NA

  for(j in 1:ncol(X)){
    if(verbose) cat(j,'')
    res.j <- association.tests(mb, X[,j], covariates=covariates,
                               drop.outliers.range=drop.outliers.range)
    for(i in 1:nrow(res.j)){
      ix <- res$mbID == res.j$mbID[i]  & res$xID == colnames(X)[j]
      headers <- c('pvalue','direction','statistic')
      res[ix,headers] <- res.j[i,headers]
    }
  }
  if(verbose) cat('\n')
  
  
  qvals <- p.adjust(drop.na(res[,'pvalue']),method='fdr')
  res[!is.na(res[,'pvalue']),'qvalue'] <- qvals
  
  res <- res[order(res[,'pvalue']),]
  res <- res[order(res[,'qvalue']),]
  
  return(res)
}


"association.tests.validation" <- function(mb, x, covariates){
  # ensure no constant-valued columns in mb
  for(i in 1:ncol(mb)){
    if(var(mb[,i]) == 0) {
      stop(sprintf('Column %d in microbiome data is constant.\n',i))
    }
  }
  
  
  # TO DO
  # check x, covariates
  
}

"load.dge.results" <- function(fp){
  res <- read.table(fp,sep='\t',head=F,row=NULL)
  dim(res)
  colnames(res) <- c('gwasID','mwasID','logFC','logCPM','LR','pvalue','qvalue')
  res$qvalue <- p.adjust(res$pvalue)
  res <- res[order(res$pvalue),]
  return(res)
}