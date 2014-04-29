# returns a matrix of snp, taxon, pval, qval for each pair 
# tests either 0 vs 1/2 or 0/1 vs 2
"np.tests" <- function(gwas, mwas, covariates=NULL, max.p=1e-10, test.against=c(0,2)[1],
            drop.outliers=TRUE,verbose=FALSE,print.results=TRUE,
            collapse.at=0.75){
    if(nrow(gwas) != nrow(mwas)){
        stop(sprintf('Error: gx has %d rows, mb has %d rows, covariates has %d rows.\n',
            nrow(gwas), nrow(mwas), nrow(covariates)))
    }
    if(!is.null(covariates)){
        if(nrow(gwas) != nrow(covariates)){
            stop(sprintf('Error: gx has %d rows, mb has %d rows, covariates has %d rows.\n',
                nrow(gwas), nrow(mwas), nrow(covariates)))
        }
        mm <- covariate.model.matrix(covariates)
    }
    # drop columns with only 1 group of 2 or more
    drop.ix <- apply(gwas,2,function(xx) {counts <- table(xx==test.against); sum(counts >= 2) < 2})
    gwas <- gwas[,!drop.ix]
    
    if(is.null(colnames(mwas))) colnames(mwas) <- sprintf('col%d',1:ncol(mwas))

    # collapse if requested
    # cluster taxa and get one rep
    if(!is.null(collapse.at)){
        # collapse mwas
        mwas.groups <- cluster.by.correlation(mwas, collapse.at)
        names(mwas.groups) <- colnames(mwas)
        reps <- collapse.by.correlation(mwas, collapse.at, "mean", verbose=TRUE)
        mwas <- mwas[,reps]
        # collapse gwas
        gwas.groups <- cluster.by.correlation(gwas, collapse.at)
        names(gwas.groups) <- colnames(gwas)
        reps <- collapse.by.correlation(gwas, collapse.at, "mean", verbose=TRUE)
        gwas <- gwas[,reps]
    }

    N <- ncol(gwas)
    M <- ncol(mwas)
    res <- data.frame(gwasID=rep(colnames(gwas),each=M),
                      mwasID=rep(colnames(mwas),N),
                      pvalue=numeric(N*M),
                      qvalue=numeric(N*M),
                      stringsAsFactors=FALSE)
    
     for(i in 1:N){
        if(verbose & i %% 100 == 0) cat(i,' ',sep='')
        non.na <- !is.na(gwas[,i])
        group1 <- which(gwas[non.na,i] == test.against)
        group2 <- which(gwas[non.na,i] != test.against)
        start.index <- (i - 1) * M
        for(j in 1:M){
            y <- mwas[non.na,j]
            if(!is.null(covariates)){
                mlm <-lm(y ~ mm[non.na,,drop=F])
                y <- mlm$residuals;
            }
            res[start.index + j,'pvalue'] <- wilcox.test(y[group1],y[group2],exact=FALSE)$p.value
        }
    }
    if(verbose) cat('\n')
    res[,'qvalue'] <- p.adjust(res[,'pvalue'],method='fdr')
    res <- res[order(res[,'pvalue']),]
    res <- res[order(res[,'qvalue']),]
    
    if(print.results){
        sig <- res[,'qvalue'] < .1
        if(any(sig)){
            print(res[sig,,drop=F])
        } else {
            print(res[1:5,,drop=F])
        }
    }
    if(is.null(collapse.at)){
        return(res)
    } else {
        return(list(assoc=res,mb.groups=mwas.groups,gx.groups=gwas.groups))
    }
}