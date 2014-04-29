# returns a matrix of snp, taxon, pval, qval for each pair 
# tests either 0 vs 1/2 or 0/1 vs 2
# covariates must be a data frame
#
# if collapse.at is not NULL, collapses gwas + mwas 
# by clustering at collapse.at spearman's correlation
# return value has association matrix with "direction" column,
# this shows either the slope if linear, or the difference in means
# from group 2 to group 1
# perm.x means to permute the genomics
# perm.y means to permute the qtr
# max.pct.zero = 0.1 means taxa with > 0.1 pct zero will be tested with nonparametric
"association.tests" <- function(gwas, mwas, covariates=NULL, top.n.results=NULL,
            verbose=FALSE,print.results=TRUE,warn=TRUE,
            collapse.at=NULL, test.type=c('np','linear','cor','ttest','linear-nonzero','logistic','chisq','perm.x','perm.y','negbin','poisson')[1], test.against=c(0,2)[1],
            drop.outliers.range=NULL,
            do.power.transform=FALSE,
            nperm=999, use.qvalue=FALSE, zi.basis=2000,
            max.pct.zero=NULL,
            ncores=1){

    if(ncores > 1) library('multicore')
    if(is.null(dim(gwas))) gwas <- matrix(gwas, ncol=1)
    if(is.null(dim(mwas))) mwas <- matrix(mwas, ncol=1)
    if(!is.null(covariates) & is.null(dim(covariates))) covariates <- data.frame(covariate1=covariates)
    
    if(nrow(gwas) != nrow(mwas)){
        stop(sprintf('Error: gx has %d rows, mb has %d rows.\n',
            nrow(gwas), nrow(mwas)))
    }
    if(!is.null(covariates)){
        na.ix <- colMeans(is.na(covariates)) > .5
        # drop any metadata columns with > .5 NA
        if(any(na.ix)){
            
            cat(sprintf('Warning: dropping metadata column %s due to >50 %% NA.\n',colnames(covariates)[na.ix]))
            covariates <- covariates[,!na.ix,drop=F]
        }
        na.ix <- rowSums(is.na(covariates)) > 0
        if(any(na.ix)){
         if(warn || mean(na.ix) > .05) cat(sprintf('Warning: dropping %d samples due to NA values in covariates.\n',sum(na.ix)))
            gwas <- gwas[!na.ix,,drop=F]
            mwas <- mwas[!na.ix,,drop=F]
        }
        covariates <- droplevels(covariates[!na.ix,,drop=F])
        if(nrow(gwas) != nrow(covariates)){
            stop(sprintf('Error: gx has %d rows, mb has %d rows, covariates has %d rows.\n',
                nrow(gwas), nrow(mwas), nrow(covariates)))
        }
        mm <- covariate.model.matrix(covariates)
    }

    # drop columns with only 1 group of 2 or more
    # unless test.against is null (i.e. kruskal-wallis)
    min.n.per.group <- 3
    if(test.type == 'np' | test.type == 'chisq'){
        if(test.type == 'np' & !is.null(test.against)){
            drop.ix <- apply(gwas,2,function(xx) {counts <- table(xx==test.against); any(counts < min.n.per.group)})
        } else {
            drop.ix <- rep(FALSE,ncol(gwas))
        }
        
    } else {
        if(verbose) cat('Finding gx columns with fewer than 3 subjects with 0, 1, or 2 alleles...\n')
        drop.ix <- apply(gwas,2,function(xx) {
                # only test for dropping if is a 0,1,2 column
                if(mean(xx %in% 0:2) == 1){
                    counts <- table(xx); 
                    any(counts < min.n.per.group)
                } else {
                    FALSE
                }
            })
        if(verbose) cat('Setting columns with too few "2" entries to 0=0, 1=1 or 2\n')
        for(i in which(drop.ix)){
            if(sum(gwas[,i] == 2) > 0 & sum(gwas[,i] == 2) < min.n.per.group){
                gwas[gwas[,i]==2,i] <- 1
                drop.ix[i] <- FALSE
            }
        }
    }
    
    if(!any(!drop.ix)) stop('Error: no snps have at least 3 subjects in each group.\n')
    if(verbose) cat(sprintf('Dropping %d columns, %d remain.\n',sum(drop.ix), sum(!drop.ix)))
    gwas <- gwas[,!drop.ix,drop=F]

    if(test.type == 'chisq'){
        # ensure >= 10 samples in each 'bin' of mwas
        drop.ix <- apply(mwas,2, function(xx) {counts <- table(xx); any(counts < min.n.per.group) || length(counts) < 2})
        if(!any(!drop.ix)) stop('Error: no qtrs have at least 3 subjects in each group.\n')
        if(verbose || warn) cat(sprintf('Dropping %d qtr columns, %d remain.\n',sum(drop.ix), sum(!drop.ix)))
        mwas <- mwas[,!drop.ix,drop=F]
    }

    if(test.type == 'zi-negbin'){
        # ensure pct nonzero >= .1
        drop.ix <- apply(mwas,2, function(xx) {mean(xx>0) < .1})
        if(!any(!drop.ix)) stop('Error: no qtrs have at least 3 subjects in each group.\n')
        if(verbose || warn) cat(sprintf('Dropping %d qtr columns, %d remain.\n',sum(drop.ix), sum(!drop.ix)))
        mwas <- mwas[,!drop.ix,drop=F]
    }

    if(is.null(colnames(mwas))) colnames(mwas) <- sprintf('col%d',1:ncol(mwas))
    if(is.null(colnames(gwas))) colnames(gwas) <- sprintf('col%d',1:ncol(gwas))

    # collapse if requested
    # cluster taxa and get one rep
    if(!is.null(collapse.at)){
        # collapse mwas
        mwas.groups <- cluster.by.correlation(mwas, collapse.at)
        names(mwas.groups) <- colnames(mwas)
        reps <- collapse.by.correlation(mwas, collapse.at, "mean", verbose=TRUE)
        mwas <- mwas[,reps,drop=F]
        # collapse gwas
        gwas.groups <- cluster.by.correlation(gwas, collapse.at)
        names(gwas.groups) <- colnames(gwas)
        reps <- collapse.by.correlation(gwas, collapse.at, "mean", verbose=TRUE)
        gwas <- gwas[,reps,drop=F]
    }

    N <- ncol(gwas)
    M <- ncol(mwas)
    if(M == 0) stop(sprintf('Error: no samples left in quant trait matrix.\n'))
    if(N == 0) stop(sprintf('Error: no samples left in gx matrix.\n'))
    if(is.null(top.n.results)){
        if(verbose) cat(sprintf('Creating results matrix, %d x %d = %d entries...\n',N,M,N*M))
        res <- data.frame(gwasID=rep(colnames(gwas),M),
                          mwasID=rep(colnames(mwas),each=N),
                          pvalue=numeric(N*M),
                          qvalue=numeric(N*M),
                          direction=numeric(N*M),
                          statistic=numeric(N*M),
                          stringsAsFactors=FALSE)
        res[,c('pvalue','qvalue')] <- NA
    } else {
        if(top.n.results > N * M) top.n.results <- N * M
        if(verbose) cat(sprintf('Creating results matrix, %d entries...\n',top.n.results))
        res <- list(gwasID=character(top.n.results),
                    mwasID=character(top.n.results),
                    pvalue=rep(2,top.n.results),
                    qvalue=numeric(top.n.results),
                    direction=numeric(top.n.results),
                    statistic=numeric(top.n.results)
                )
        res <- as.data.frame(res)        
        res[,'qvalue'] <- NA
    }
    
    if(verbose) cat('Of ',M,': ',sep='')
    
    outlier <- rep(FALSE,nrow(mwas))
    total.zi.tests <- 0 # keep track of total tests performed for zi-negbin
    total.failures <- 0
    mcres <- mclapply(1:M,function(j){
        if(verbose && j %% 1 == 0) cat('***',j,'***')

        mcres.j <- data.frame(
                    pvalue=numeric(N),
                    direction=numeric(N),
                    statistic=numeric(N),
                    stringsAsFactors=FALSE
                )
        mcres.j[,c('pvalue')] <- NA

        
        if(!is.null(drop.outliers.range)) outlier <- is.outlier(mwas[,j],range=drop.outliers.range)
        for(i in 1:N){
			if(verbose & i %% 1 == 0) cat('',i,'')
            # drop subjects with na in this gene (or outliers)
            if(!is.null(drop.outliers.range)){
                if(test.type == 'linear-nonzero' | test.type == 'np-nonzero' | test.type == 'zi-negbin'){            
                    outlier <- rep(FALSE,nrow(mwas))
                    nz <- mwas[,j] > 0
                    outlier[nz] <- is.outlier(mwas[nz,j],range=drop.outliers.range)
                } else {
                    outlier <- is.outlier(mwas[,j],range=drop.outliers.range)
                }            
                # only measure outliers within nonzero portions if linear-nonzero
            }
            non.na <- !is.na(gwas[,i]) & !outlier
            if(any(is.na(non.na))) next
            if(sum(non.na) < 3) next
            if(var(mwas[non.na,j]) == 0) mcres.j[i,'pvalue'] <- NA
            direction <- NA
            statistic <- NA
            if(test.type == 'np' || test.type == 'np-nonzero'){
                nz <- rep(TRUE,nrow(mwas))
                if(test.type == 'np-nonzero') nz <- mwas[,j] > 0
                if(!is.null(test.against)){
                    group1 <- which(gwas[non.na & nz,i] == test.against)
                    group2 <- which(gwas[non.na & nz,i] != test.against)
                }
                y <- mwas[non.na & nz,j]
                if(!is.null(covariates)){
                    if(do.power.transform){
                        require('car')
                        p1 <- powerTransform(y + min(y[y>0])/2 ~ mm[non.na & nz,,drop=F])
                        # fit linear model with transformed response:
                        mlm <- lm(bcPower(y + min(y[y>0])/2, p1$roundlam) ~ mm[non.na & nz,,drop=F])
                    } else {
                        mlm <-lm(y ~ mm[non.na & nz,,drop=F])
                    }
                    y <- mlm$residuals;
                }
                if(!is.null(test.against)){
                    pval <- wilcox.test(y[group1],y[group2],exact=FALSE)$p.value
                    direction <- ifelse(median(y[group2]) > median(y[group1]),1,-1)
                } else {
                    pval <- kruskal.test(y ~ as.factor(gwas[non.na & nz,i]))$p.value
                }
            } else if(test.type == 'perm.x' | test.type == 'perm.y' | test.type == 'perm.xy'){
                y <- mwas[non.na,j]
                
                if(!is.null(test.against)) {
                    # if this is binary, create a grouping indicator
                    x <- gwas[non.na,i] == test.against
                    direction <- ifelse(median(y[!x]) > median(y[x]),1,-1)
                    max.W <- prod(table(x)) # maximum statistic
                } else {
                    # else keep the grouping as a factor for kruskal's
                    x <- gwas[non.na,i]
                }
                W <- numeric(nperm + 1) # will hold all perm test statistics + actual
                # run permutations
                for(ii in 1:(nperm+1)){
                    if(verbose && ii %% 10 == 0) cat('.')
                    if(!is.null(covariates)){
                        if(do.power.transform){
                            require('car')
                            p1 <- powerTransform(y + min(y[y>0])/2 ~ mm[non.na,,drop=F])
                            if(p1$convergence > 0){
                                # if power transform failed, use pseudolog
                                mlm <- lm(pseudolog(y) ~ mm[non.na,,drop=F])
                            } else {
                                # fit linear model with transformed response:
                                mlm <- lm(bcPower(y + min(y[y>0])/2, p1$roundlam) ~ mm[non.na,,drop=F])
                            }
                        } else {
                            mlm <-lm(y ~ mm[non.na,,drop=F])
                        }
                        y.i <- mlm$residuals;
                    } else {
                        y.i <- power.transform(y)
                    }

                    if(!is.null(test.against)){
                        W[ii] <- wilcox.test(y.i[x],y.i[!x],exact=FALSE)$statistic
                        W[ii] <- max(W[ii],max.W - W[ii]) # make this two-sided since W is 0-max.W
                    } else {
                        W[ii] <- kruskal.test(y.i ~ as.factor(x))$statistic
                    }
                    if(test.type %in% c('perm.x','perm.xy')) x <- sample(x)
                    if(test.type %in% c('perm.y','perm.xy')) y <- sample(y)
                }
                pval <- mean(W >= W[1])
            } else if(test.type == 'ttest'){
                if(!is.null(test.against)){
                    group1 <- which(gwas[non.na,i] == test.against)
                    group2 <- which(gwas[non.na,i] != test.against)
                }
                y <- mwas[non.na,j]
                if(!is.null(covariates)){
                    mlm <-lm(y ~ mm[non.na,,drop=F])
                    y <- mlm$residuals;
                }
                if(!is.null(test.against) || length(unique(gwas[non.na,i]))){
                    if(var(y[group1]) == 0 && var(y[group2]) == 0){
                        pval <- 1
                    } else {
                        pval <- t.test(y[group1],y[group2],exact=FALSE)$p.value
                    }
                } else {
                    fit <- aov(y ~ as.factor(gwas[non.na,i]))
                    pval <- summary(fit)[[1]][1,5]
                }
            } else if((test.type == 'hurdle' | test.type == 'hurdlezero') & min(mwas[non.na,j]) == 0){
                require('pscl')
                # if one class has non-zeros but not the other
                # or if metadata prevents testing
                if(any(sapply(split(mwas[non.na,j],gwas[non.na,i,drop=F]),function(xx) sum(xx>0) == 0)) ||
                    any(apply(mm[non.na,,drop=F],2,function(xx) sum(table(mwas[non.na,j]==0,xx==0)==0)) > 0)){
                    
#                     pval <- prop.test(table(gwas[non.na,i,drop=F] == 0, mwas[non.na,j]==0))$p.value
                    pval <- NA
                } else {
                    if(is.null(covariates)){
                        mlm <- hurdle(mwas[non.na,j] ~ gwas[non.na,i,drop=F],dist='negbin')
                    } else {
                        mlm <- hurdle(mwas[non.na,j] ~ gwas[non.na,i,drop=F] + mm[non.na,,drop=F] | gwas[non.na,i,drop=F],dist='negbin')
                    }
                    # note: should eventually return test for zero portion
                    if(test.type == 'hurdle'){
                        pval <- summary(mlm)[[1]]$count[2,4]
                    } else {
                        pval <- summary(mlm)[[1]]$zero[2,4]
                    }
                }
            } else if(test.type == 'chisq'){
                # assumes there are two classes in mwas (e.g. zero/nonzero)
                pval <- chisq.test(table(gwas[non.na,i], mwas[non.na,j]), simulate.p.value=TRUE)$p.value
                direction <- NA
            } else if(test.type == 'logistic'){
				y <- as.numeric(mwas[non.na,j] > 0)
                if(is.null(covariates)){
                    mlm <- glm(y ~ gwas[non.na,i,drop=F],family='binomial')
                } else {
                    mlm <- glm(y ~ matrix(gwas[non.na,i,drop=F],ncol=1) + mm[non.na,,drop=F],family='binomial')
                }
                if(!any(grepl('gwas',rownames(summary(mlm)[[12]])))){
                    pval <- NA
                } else {
                    pval <- summary(mlm)[[12]][2,4]
                }
            } else if(test.type == 'linear-nonzero'){
                # recalculate outliers based on nonzero values
                outlier <- rep(FALSE, nrow(mwas))
                if(!is.null(drop.outliers.range)) outlier[nz] <- is.outlier(mwas[nz,j],range=drop.outliers.range)
                non.na <- !is.na(gwas[,i]) & !outlier
                nz <- mwas[non.na,j] > 0

                if(is.null(covariates)){
                    if(do.power.transform){
                        require('car')
                        p1 <- powerTransform(mwas[non.na,j][nz] ~ gwas[non.na,i,drop=F][nz,,drop=F])
                        if(p1$convergence > 0){
                            # if power transform failed, use pseudolog
                            mlm <- lm(pseudolog(mwas[non.na,j][nz]) ~ gwas[non.na,i,drop=F][nz,,drop=F])
                        } else {
                            # fit linear model with transformed response:
                            mlm <- lm(bcPower(mwas[non.na,j][nz], p1$roundlam) ~ gwas[non.na,i,drop=F][nz,,drop=F])
                        }
                    } else {
                        mlm <- lm(mwas[non.na,j][nz] ~ gwas[non.na,i,drop=F][nz,,drop=F])
                    }
                } else {

                    if(do.power.transform){
                        require('car')
                        p1 <- powerTransform(mwas[non.na,j][nz] ~ gwas[non.na,i,drop=F][nz,,drop=F] + mm[non.na,,drop=F][nz,,drop=F])
                        if(p1$convergence > 0){
                            # if power transform failed, use pseudolog
                            mlm <- lm(pseudolog(mwas[non.na,j][nz]) ~ gwas[non.na,i,drop=F][nz,,drop=F] + mm[non.na,,drop=F][nz,,drop=F])
                        } else {
                            # fit linear model with transformed response:
                            mlm <- lm(bcPower(mwas[non.na,j][nz], p1$roundlam) ~ gwas[non.na,i,drop=F][nz,,drop=F] + mm[non.na,,drop=F][nz,,drop=F])
                        }
                    } else {
                        mlm <- lm(mwas[non.na,j][nz] ~ gwas[non.na,i,drop=F][nz,,drop=F] + mm[non.na,,drop=F][nz,,drop=F])
                    }
                }
                if(!any(grepl('gwas',rownames(summary(mlm)[[4]])))){
                    pval <- NA
                } else {
                    pval <- summary(mlm)[[4]][2,4]
                    direction <- mlm$coefficients[grep('gwas',names(mlm$coefficients))]
                    statistic <- summary(mlm)[[4]][2,'t value']
                }
            } else if(test.type == 'linear' | (test.type == 'hurdle' & min(mwas[non.na,j]) > 0)){
                y <- mwas[non.na,j]

                if(is.null(covariates)){
                    if(do.power.transform){
                        require('car')
                        p1 <- powerTransform(y + min(y[y>0])/2 ~ gwas[non.na,i,drop=F])
                        if(p1$convergence > 0){
                            # if power transform failed, use pseudolog
                            mlm <- lm(pseudolog(y) ~ gwas[non.na,i,drop=F])
                        } else {
                            # fit linear model with transformed response:
                            mlm <- lm(bcPower(y + min(y[y>0])/2, p1$roundlam) ~ gwas[non.na,i,drop=F])
                        }
                    } else {
                        mlm <- lm(y ~ gwas[non.na,i,drop=F])
                    }
                } else {
                    if(do.power.transform){
                        require('car')
						if(any(y==0)){
							y <- y + min(y[y>0])/2 
						}
						# force optim to throw error on warnings:
						# hack to avoid "Warning in sqrt(diag(solve(res$hessian))) : NaNs produced"
						warn.value <- options()[['warn']]
						options(warn=2)
                        p1 <- try(powerTransform(y ~ gwas[non.na,i,drop=F] + mm[non.na,,drop=F], method='BFGS'),silent=TRUE)
                        if(class(p1)=='try-error'){
                        	if(verbose) cat('Error running BFGS optimizer, trying CG...\n')
                        	p1 <- try(powerTransform(y ~ gwas[non.na,i,drop=F] + mm[non.na,,drop=F], method='CG'),silent=TRUE)
                        }
                        if(class(p1)=='try-error'){
                        	if(verbose) cat('Error running CG optimizer, trying SANN...\n')
                        	p1 <- try(powerTransform(y ~ gwas[non.na,i,drop=F] + mm[non.na,,drop=F], method='SANN'),silent=TRUE)
                        }
						options(warn=warn.value)
                        if(class(p1)=='try-error'){
                        	if(verbose) cat('Error running all optimizers, defaulting to pseudolog...\n')
                        }
                        
                        if(class(p1)=='try-error' || p1$convergence != 0){
                            # if power transform failed, use pseudolog
                            mlm <- lm(pseudolog(y) ~ gwas[non.na,i,drop=F] + mm[non.na,,drop=F])
                        } else {
                            # fit linear model with transformed response:
                            mlm <- lm(bcPower(y, p1$roundlam) ~ gwas[non.na,i,drop=F] + mm[non.na,,drop=F])
                        }
                    } else {
                        if(!is.null(max.pct.zero) && mean(y==min(y)) > max.pct.zero){
                            if(all(table(gwas[non.na,i]) >= min.n.per.group)){
                                cat('Column',j,'gene',i,'>',max.pct.zero,'zeros; running np. ')
                                # assumes there are two classes in mwas (e.g. zero/nonzero)
                                wilcox.res <- wilcox.test(y[gwas[non.na,i] != min(gwas[non.na,i])], y[gwas[non.na,i] == min(gwas[non.na,i])],exact=F,conf.int=TRUE)
                                pval <- wilcox.res$p.value
                                statistic <- wilcox.res$statistic
                                direction <- wilcox.res$estimate
                            } else {
                                pval <- NA
                                direction <- NA
                                statistic <- NA
                            }
                        } else {
                            mlm <- lm(y ~ gwas[non.na,i,drop=F] + mm[non.na,,drop=F])
                        }
                    }
                }
                if(is.null(max.pct.zero) || mean(y==min(y)) <= max.pct.zero){
                    if(!any(grepl('gwas',rownames(summary(mlm)[[4]])))){
                        pval <- NA
                    } else {
                        pval <- summary(mlm)[[4]][2,'Pr(>|t|)']
                        direction <- mlm$coefficients[grep('gwas',names(mlm$coefficients))]
                        statistic <- summary(mlm)[[4]][2,'t value']
                    }
                }
            } else if(test.type == 'negbin'){
                require('MASS')
                y <- mwas[non.na,j]
                basis <- ceiling(1/min(mwas[mwas > 0],na.rm=TRUE))
                y <- ceiling(y * basis)
                
                if(is.null(covariates)){
                    # just in case, initialize with poisson fits
#                     mp <- try(glm(y ~ gwas[non.na,i,drop=F], family='poisson'))
#                     mlm <- try(glm.nb(y ~ gwas[non.na,i,drop=F],start=coef(mp)))
                    mlm <- try(glm.nb(y ~ gwas[non.na,i,drop=F]))
                } else {
#                     mp <- try(glm(y ~ gwas[non.na,i,drop=F] + mm[non.na,,drop=F], family='poisson'))
#                     mlm <- try(glm.nb(y ~ gwas[non.na,i,drop=F] + mm[non.na,,drop=F], start=coef(mp)))
                    mlm <- try(glm.nb(y ~ gwas[non.na,i,drop=F] + mm[non.na,,drop=F]))
                }
                # if regression failed, return NA
                if(class(mlm)=='try-error' || !mlm$converged){
                    if(all(table(gwas[non.na,i]) >= min.n.per.group) &&
                        all(table(y > 0) >= min.n.per.group)){
                        cat('Column',j,'failed; running chisq\n')
                        # assumes there are two classes in mwas (e.g. zero/nonzero)
                        pval <- chisq.test(table(gwas[non.na,i], mwas[non.na,j]), simulate.p.value=TRUE)$p.value
                        direction <- NA
                    } else {
                        cat('Column',j,'failed and chisq failed\n')
                        pval <- NA
                        direction <- NA
                    }
                } else {
                    coefs <- summary(mlm)[['coefficients']]
                    if(!any(grepl('gwas',rownames(coefs)))){
                        pval <- NA
                    } else {
                        pval <- coefs[2,'Pr(>|z|)']
                        direction <- mlm$coefficients[grep('gwas',names(mlm$coefficients))]
                        statistic <- coefs[2,'z value']
                    }
                }
            } else if(test.type == 'zi-negbin'){
                pval <- NA
                statistic <- NA
                direction <- NA

                require('MASS')
                require('pscl')
                y <- mwas[non.na,j]
                if(is.null(zi.basis)){
                    basis <- ceiling(1/min(mwas[mwas > 0],na.rm=TRUE))
                } else {
                    basis <- zi.basis
                }
                y <- ceiling(y * basis)
                if(is.null(covariates)){
                    stop('zi-negbin only implemented with covariates\n')
                } else {
                    # if > 10% zeros, do zeroinfl
                    if(mean(y==0) < .1){
#                         if(verbose) cat("ZI not applicable for gwas",i,"mwas", j,", too few zeros, trying linear\n")
                        mlm <- lm(asin(sqrt(mwas[non.na,j])) ~ gwas[non.na,i,drop=F]  + mm[non.na,,drop=F])
                        mlm$converged <- TRUE
#                         mlm <- try(glm.nb(y ~ gwas[non.na,i,drop=F]  + mm[non.na,,drop=F]))
                        is.zi <- FALSE
                    } else if(mean(y == 0) < 0.9){
                        # just in case, initialize with poisson and logistic
#                         mp <- try(glm(y ~ gwas[non.na,i,drop=F] + mm[non.na,,drop=F], family='poisson'))
#                         ml <- try(glm(y > 0 ~ matrix(gwas[non.na,i,drop=F],ncol=1),family='binomial'))
#                         mlm <- try(zeroinfl(y ~ gwas[non.na,i,drop=F] + mm[non.na,,drop=F] | gwas[non.na,i,drop=F], dist='negbin',EM=FALSE,
#                                     start=list(count=coef(mp), zero=coef(ml))),silent=TRUE)
                        mlm <- try(zeroinfl(y ~ gwas[non.na,i,drop=F] + mm[non.na,,drop=F] | gwas[non.na,i,drop=F], dist='negbin',EM=FALSE),silent=TRUE)
                        # make certain we can access p-values
                        # sometimes zeroinfl will succeed but is unusable
                        if(class(mlm) != "try-error"){
                            if(class(try(summary(mlm),silent=TRUE))[1] == 'try-error'){
                                class(mlm) <- 'try-error'
                            } else {
                                coefs <- summary(mlm)[['coefficients']]
                                pval.count <- coefs$count[grep('gwas',rownames(coefs$count)),'Pr(>|z|)']
                                pval.zero <- coefs$zero[grep('gwas',rownames(coefs$zero)),'Pr(>|z|)']
                                if(is.nan(pval.count) || is.nan(pval.zero)) class(mlm) <- 'try-error'
                            }
                        }                                                
                        is.zi <- TRUE
                    }

                    # run non-parametric if test failed or > 90% zeros                    
                    if(mean(y == 0) > 0.9 || class(mlm)[1]=='try-error' || !mlm$converged){
                        if(class(mlm)[1]=='try-error' || !mlm$converged) {
                            total.failures <- total.failures + 1
                            cat('Column',j,'failed; running np...')
                        } else {
                            cat('Column',j,'mean(y==0) > 0.9; running np...')
                        }
                        if(all(table(gwas[non.na,i]) >= min.n.per.group)){
                            # assumes there are two classes in mwas (e.g. zero/nonzero)
                            wilcox.res <- wilcox.test(mwas[non.na,j][gwas[non.na,i] != min(gwas[non.na,i])], mwas[non.na,j][gwas[non.na,i] == min(gwas[non.na,i])],exact=F,conf.int=TRUE)
                            pval <- wilcox.res$p.value
                            statistic <- wilcox.res$statistic
                            direction <- wilcox.res$estimate
                        } else {
                            cat('Column',j,'not testable\n')
                        }
                    } else {
                        coefs <- summary(mlm)[['coefficients']]
                        # if ran zero-infl
                        if(is.zi){
                            total.zi.tests <- total.zi.tests + 2
                            if(any(grepl('gwas',rownames(coefs$count)))){
                                pval.count <- coefs$count[grep('gwas',rownames(coefs$count)),'Pr(>|z|)']
                                pval.zero <- coefs$zero[grep('gwas',rownames(coefs$zero)),'Pr(>|z|)']
                                direction.count <- coefs$count[grep('gwas',rownames(coefs$count)),'Estimate']
                                direction.zero <- coefs$zero[grep('gwas',rownames(coefs$zero)),'Estimate']
                                statistic.count <- coefs$count[grep('gwas',rownames(coefs$count)),'z value']
                                statistic.zero <- coefs$zero[grep('gwas',rownames(coefs$zero)),'z value']
                                best.ix <- which.min(c(pval.count, pval.zero))
                                cat(c('c','z')[best.ix],'')
                                pval <- min(c(pval.count, pval.zero))
                                # for convenience, statistic holds 1 if count-based pvalue, 2 if zero-based pvalue
                                statistic <- best.ix #c(statistic.count, statistic.zero)[best.ix]
                                direction <- c(direction.count, direction.zero)[best.ix]
                            }
                        } else {
                            total.zi.tests <- total.zi.tests + 1
                            if(any(grepl('gwas',rownames(coefs)))){
                                pval <- coefs[grep('gwas',rownames(coefs)),ncol(coefs)]
                                direction <- mlm$coefficients[grep('gwas',names(mlm$coefficients))]
                                statistic <- coefs[grep('gwas',rownames(coefs)),ncol(coefs)-1]
                            }
                        }
                    }
                }
            } else if(test.type == 'poisson'){
                y <- mwas[non.na,j]
                basis <- ceiling(1/min(mwas[non.na,j][mwas[non.na,j] > 0]))
                y <- round(y * basis)
                if(is.null(covariates)){
                    mlm <- glm(y ~ gwas[non.na,i,drop=F], family ='poisson')
                } else {
                    mlm <- glm(y ~ gwas[non.na,i,drop=F] + mm[non.na,,drop=F], family ='poisson')
                }
#                 cat(j,'\n\n\n')
                
                coefs <- summary(mlm)[['coefficients']]
                if(!any(grepl('gwas',rownames(coefs)))){
                    pval <- NA
                } else {
                    pval <- coefs[2,'Pr(>|z|)']
                    direction <- mlm$coefficients[grep('gwas',names(mlm$coefficients))]
                    statistic <- coefs[2,'z value']
                }
            } else if(test.type == 'cor'){
                y <- mwas[non.na,j]
                if(!is.null(covariates)){
                    mlm <-lm(y ~ mm[non.na,,drop=F])
                    y <- mlm$residuals;
                }
                pval <- cor.test(gwas[non.na,i], y)$p.value
            }
#             if(verbose){ cat('TOTAL FAILURES:',total.failures,'\n')}
            mcres.j[i,'pvalue'] <- pval
            mcres.j[i,'direction'] <- direction
            mcres.j[i,'statistic'] <- statistic
            if(verbose && i == N) cat(sprintf('j=%d, i=%d, min adjusted p-value so far: %f\n',j, i, min(drop.na(mcres.j[,'pvalue'])) * N * M))
        }
        return(mcres.j)
    }, mc.cores=ncores)

    # process results
    if(is.null(top.n.results)){
        for(j in 1:M){
            start.index <- (j - 1) * N
            res[start.index + (1:N),'pvalue'] <- mcres[[j]][,'pvalue']
            res[start.index + (1:N),'direction'] <- mcres[[j]][,'direction']
            res[start.index + (1:N),'statistic'] <- mcres[[j]][,'statistic']
        }
    } else {
        for(j in 1:M){
            start.index <- (j - 1) * N
            for(i in 1:N){
                # add this to results if it's smaller than any previous
                pval <- mcres[[j]][i,'pvalue']
                if(!is.na(pval) & !is.nan(pval)){
                    if(max(res[,'pvalue']) > 1 || pval < max(res[,'pvalue'])){
                        ix <- which.max(res[,'pvalue'])
                        res[ix,'pvalue'] <- pval
                        res[ix,'direction'] <- mcres[[j]][i,'direction']
                        res[ix,'statistic'] <- mcres[[j]][i,'statistic']
                        res[ix,'gwasID'] <- colnames(gwas)[i]
                        res[ix,'mwasID'] <- colnames(mwas)[j]
                    }
                }
            }
        }
    }
    
    total.tests <- N * M
    if(test.type == 'zi-negbin') total.tests <- total.zi.tests

    if(verbose) cat('\n')
	ix <- !is.na(res[,'pvalue'])
	if(use.qvalue){
        require('qvalue')
        qvals <- try(qvalue(c(res[ix,'pvalue'],rep(1-10^-5,total.tests - sum(ix))))$qvalue[1:sum(ix)],silent=TRUE)
    }
    if(!use.qvalue || class(qvals)=='try-error'){
		qvals <- p.adjust(drop.na(res[,'pvalue']),method='fdr',n=total.tests)
    }
    res[ix,'qvalue'] <- qvals

    res <- res[order(res[,'pvalue']),]
    res <- res[order(res[,'qvalue']),]
    

    if(print.results){
         max.print <- 10
         min.print <- 5
         sig <- which(res[,'qvalue'] < .1)
        print(res[1:min(10,nrow(res)),,drop=F])         
#         if(length(sig) > 0){
#             print(res[sig[1:min(length(sig),max.print)],,drop=F])
#             if(length(sig) > max.print){
#                 cat(sprintf('+ %d more...\n',length(sig) - max.print))
#             }
#         } else {
#             print(res[1:5,,drop=F])
#         }
    }
    if(is.null(collapse.at)){
        return(list(assoc=res))
    } else {
        return(list(assoc=res,mb.groups=mwas.groups,gx.groups=gwas.groups,total.failures=total.failures))
    }
}


# Run stratified tests, combined them with fisher's method
# *.subsets are lists of binary subset indicators for each stratum
# if not a list, then just a single indicator vector for all strata
# subject.subsets must be a list, or NULL if no subsets
#
# covariates can either be a single matrix/data.frame or a list of matrices 
# if stratum.names is NULL, names taken from subject.subsets
"stratified.association.tests" <- function(gx, qtr, 
        subject.subsets=NULL, snp.subsets=NULL, qtr.subsets=NULL, covariates=NULL,
        stratum.names=NULL, test.type=c('np','linear')[1],
        test.against=0, drop.outliers.range=NULL,
        print.results=TRUE, warn=TRUE, verbose=FALSE,
        top.n.results=NULL, do.power.transform=FALSE, nperm=999, use.qvalue=FALSE,
        max.pct.zero=NULL,
        ncores=1){
        
    results <- list()
    if(is.null(dim(gx))) gx <- matrix(gx,ncol=1)
    if(is.null(dim(qtr))) qtr <- matrix(qtr,ncol=1)
    if(is.null(colnames(qtr))) colnames(qtr) <- sprintf('col%d',1:ncol(qtr))
    if(is.null(colnames(gx))) colnames(gx) <- sprintf('col%d',1:ncol(gx))

    if(is.null(subject.subsets)) subject.subsets <- list('No_subsets'=rep(TRUE,nrow(gx)))
    if(is.null(stratum.names)) stratum.names <- names(subject.subsets)
    if(is.null(snp.subsets)) snp.subsets <- rep(TRUE,ncol(gx))
    if(is.null(qtr.subsets)) qtr.subsets <- rep(TRUE,ncol(qtr))
    for(i in 1:length(subject.subsets)){

        snp.subset <- snp.subsets
        if(class(snp.subset) == 'list') snp.subset <- snp.subsets[[i]]
        qtr.subset <- qtr.subsets
        if(class(qtr.subset) == 'list') qtr.subset <- qtr.subsets[[i]]
        covariates.i <- NULL
        if(!is.null(covariates)) {
            covariates.i <- covariates[subject.subsets[[i]],,drop=F]
        }
        if(length(subject.subsets) > 1 && verbose) cat('\n',stratum.names[i],'\n',sep='')
        results[[stratum.names[i]]] <- 
                association.tests(gx[subject.subsets[[i]], snp.subset, drop=F],
                                  qtr[subject.subsets[[i]], qtr.subset, drop=F],
                                  verbose=verbose,
                                  covariates=covariates.i,
                                  test.type=test.type,
                                  test.against=test.against, drop.outliers.range=drop.outliers.range,
                                  print.results=length(subject.subsets)==1 & print.results,
                                  warn=warn, top.n.results=top.n.results,
                                  do.power.transform=do.power.transform, nperm=nperm,
                                  use.qvalue=use.qvalue,
                                  max.pct.zero=max.pct.zero,
                                  ncores=ncores)$assoc
    }
    # combine results and do fisher's
    if(length(results) > 1){
        cat('Combining results with fisher\'s method...\n')
        combined <- combine.results(results, colnames(gx), colnames(qtr))
        if(print.results) print(combined[1:3,])
    } else {
        combined <- results[[1]]
    }

    return(combined)
}

# assoc has gx id in first column, quant trait id in second column
# ... sent to plot
"plot.association.results" <- function(gx, qtr, assoc, qtr.names=NULL, covariates=NULL,
        filebase=NULL,drop.outliers.range=NULL,add.to.plot=FALSE, ymin=NULL,
        shorten.qtr.names=TRUE, shorten.qtr.names.delim=';',
        regression.type='linear',
        color.by=NULL,legend.names=NULL,
        color.by.color=NULL,
        xlab=NULL,
        plot.type=c('vioplot','boxplot','density','barplot')[1],
        nonzero.only=FALSE,
        do.power.transform=FALSE){
    if(is.null(dim(gx))) gx <- matrix(gx,ncol=1)
    if(is.null(colnames(qtr))) colnames(qtr) <- sprintf('col%d',1:ncol(qtr))
    if(is.null(colnames(gx))) colnames(gx) <- sprintf('col%d',1:ncol(gx))
    
     if(shorten.qtr.names){
         shortnames <- strsplit(colnames(qtr),shorten.qtr.names.delim)
         shortnames <- sapply(shortnames,function(xx) paste(rev(rev(xx)[1:min(2,length(xx))]),collapse=' | '))
     } else {
        shortnames <- colnames(qtr)
     }

    if(nrow(gx) != nrow(qtr)){
        stop(sprintf('Error: gx has %d rows, mb has %d rows, covariates has %d rows.\n',
            nrow(gx), nrow(qtr), nrow(covariates)))
    }
    non.na <- rep(TRUE, nrow(qtr))
    if(!is.null(covariates)){
        non.na <- rowSums(is.na(covariates)) == 0
        if(nrow(gx) != nrow(covariates)){
            stop(sprintf('Error: gx has %d rows, mb has %d rows, covariates has %d rows.\n',
                nrow(gx), nrow(qtr), nrow(covariates)))
        }
        mm <- covariate.model.matrix(covariates[non.na,,drop=FALSE])
    }

    if(is.null(dim(assoc))) assoc <- matrix(assoc, nrow=1)
    
    if(!is.null(filebase)) pdf(paste(filebase,'.pdf',sep=''), width=8,height=8)        
    n <- nrow(assoc)
    nc <- ceiling(sqrt(n)) # n rows in plot
    nr <- ceiling(sqrt(n)) # n cols in plot
#     if(n == 2) nr <- 1
#     nr <- ceiling(n / nc) # n cols in plot
    if(!add.to.plot) par(mfrow=c(nr, nc))
    outlier <- rep(FALSE, nrow(qtr))
    if(is.null(color.by)){
        cols <- rep('#99999944',nrow(qtr))
     } else {
        if(is.null(color.by.color)){
            cols <- sprintf('%s99',c("999999","#81ADD3","#649D2C","#7D661C","#7E2712","#380024",brewer.pal(9,'Set1')[-6]))            
        } else {
            cols <- color.by.color
        }
        cols <- cols[as.numeric(color.by)]
     }
    
    for(i in 1:nrow(assoc)){
        if(nonzero.only){
            # if nonzero only, include zeros as "outliers"
            outlier <- qtr[,assoc[i,2]] == 0
            nz <- !outlier
            if(!is.null(drop.outliers.range)){
                outlier[nz] <- is.outlier(qtr[nz,assoc[i,2]],range=drop.outliers.range)
            }
        } else if(!is.null(drop.outliers.range)){
            outlier <- is.outlier(qtr[,assoc[i,2]],range=drop.outliers.range)
        }

        if(!is.null(covariates)){
            if(regression.type == 'linear'){
				if(do.power.transform){
					require('car')
					y <- qtr[!outlier & non.na,assoc[i,2]]
					if(any(y==0)){
						y <- y + min(y[y>0])/2 
					}
					p1 <- powerTransform(y ~ mm[!outlier & non.na,,drop=F], method='BFGS')
					if(p1$convergence > 0){
						# if power transform failed, use pseudolog
						mlm <- lm(pseudolog(y) ~ mm[!outlier & non.na,,drop=F])
					} else {
						# fit linear model with transformed response:
						mlm <- lm(bcPower(y, p1$roundlam) ~ + mm[!outlier & non.na,,drop=F])
					}
				} else {
	                y <- lm(y ~ mm[!outlier[non.na],,drop=F])$residuals
				}
            } else if(regression.type == 'negbin') {
                require('MASS')
                y <- glm.nb(qtr[!outlier & non.na,assoc[i,2]] ~ mm[!outlier[non.na],,drop=F])$residuals
            } else {
                stop('Unknown regression type\n')
            }
        } else {
            y <- qtr[!outlier & non.na,assoc[i,2]]
        }
        

        x <- gx[!outlier & non.na,assoc[i,1]]
        na.i.ix <- is.na(y) | is.na(x)
        y <- y[!na.i.ix]
        x <- x[!na.i.ix]
        
        # if more than 10 unique values, treat as real
        is.discrete <- length(unique(x)) < 5
        ylim <- range(y)
        if(!is.null(ymin)) ylim[1] <- 0
        if(is.null(xlab)) xlab.i <- assoc[i,1]
        legend.space <- 0
        if(!is.null(legend.names)) legend.space <- .85
        if(is.discrete){
            
            x <- as.factor(x)
            if(plot.type == 'barplot') {
                bh <- barplot(sapply(split(y,x),mean), add=FALSE)
                my.error.bars(bh,sapply(split(y,x),mean),sapply(split(y,x),function(xx) sd(xx)/sqrt(length(xx))))
                
            } else if(plot.type == 'density'){
                    densities <- list()
                    all.y <- NULL
                    all.x <- NULL
                    for(j in 1:length(unique(x))){
                        densities[[j]] <- density(y[x==levels(x)[j]],adjust=.5)
                        all.y <- c(all.y,densities[[j]]$y)
                        all.x <- c(all.x,densities[[j]]$x)
                    }
                    
                    plot(density(y),type='n',xlim=range(all.x),ylim=range(all.y))
                    for(j in 1:length(unique(x))){
                        lines(densities[[j]],col=bong.cols[j],lwd=3)
                    }
            } else if(plot.type == 'boxplot'){
                stripchart(y ~ x, vert=T, method='jitter',
                        pch=21,bg='#FFFFFF00', col='#FFFFFF00',
                        xlab=xlab.i, ylab='Residual log relative abundance',
                        main=sprintf('%s\n(p=%f, q=%f)',shortnames[match(assoc[i,2],colnames(qtr))], assoc[i,'pvalue'], assoc[i,'qvalue']),
                        cex.main=.7, 
                        cex.lab=.75,ylim=ylim, at=1:length(unique(x)), add=FALSE,
                        xlim=c(0.5,length(unique(x))+.5+legend.space));
                 boxplot(y ~ x,pch=NA,add=T)
                
            } else {
                beeswarm <- FALSE
                if(beeswarm){
                    require('beeswarm')
                    beeswarm(y ~ x,
                            pch=21,bg='#FFFFFF00', col='#FFFFFF00',
                            xlab=xlab.i, ylab='Residual log relative abundance',
                            main=sprintf('%s\n(p=%f, q=%f)',shortnames[match(assoc[i,2],colnames(qtr))], assoc[i,'pvalue'], assoc[i,'qvalue']),
                            cex.main=.7, 
                            cex.lab=.75,ylim=ylim, at=1:length(unique(x)), add=FALSE,
                            xlim=c(0.5,length(unique(x))+.5+legend.space));
                } else {
                    stripchart(y ~ x, vert=T, method='jitter',
                            pch=21,bg='#FFFFFF00', col='#FFFFFF00',
                            xlab=xlab.i, ylab='Residual log relative abundance',
                            main=sprintf('%s\n(p=%f, q=%f)',shortnames[match(assoc[i,2],colnames(qtr))], assoc[i,'pvalue'], assoc[i,'qvalue']),
                            cex.main=.7, 
                            cex.lab=.75,ylim=ylim, at=1:length(unique(x)), add=FALSE,
                            xlim=c(0.5,length(unique(x))+.5+legend.space));
                }
                if(plot.type=='vioplot'){
                    require('vioplot')
                    for(j in 1:length(unique(x))){
                        vioplot(y[x==levels(x)[j]],at=j,add=TRUE,col='#00000000',drawRect=FALSE)
                    }
                }
                if(beeswarm){
                    beeswarm(y ~ x,
                            pch=21,bg=cols[!outlier & non.na][!na.i.ix], col='#00000044',
                            cex.lab=.75,ylim=ylim, at=1:length(unique(x)), add=TRUE,
                            );
                } else {
                    points(jitter(as.numeric(as.factor(x))),y,
#                             pch=21,bg=cols[!outlier & non.na][!na.i.ix], col='#00000044',
                            pch=16,col=cols[!outlier & non.na][!na.i.ix],
                            cex.lab=.75,ylim=ylim,
                            );

                }
                
                # add median
                for(j in 1:length(unique(x))){
                    lines(c(j-.25,j+.25),rep(mean(y[x==levels(x)[j]]),2),col='#000000FF',lwd=2)
                }
            }                
        } else {
            plot(x, y, pch=21, bg='#00000044', col='#00000099',
                    xlab=xlab.i, ylab='Residual log relative abundance',
                    main=sprintf('%s\n(p=%f, q=%f)',shortnames[match(assoc[i,2],colnames(qtr))], assoc[i,'pvalue'], assoc[i,'qvalue']),
                    cex.main=.7, 
                    cex.lab=.75,ylim=ylim);
            abline(lm(y ~ x))
#             print(cor.test(x,y))
        }
        
        # draw legend for colorby values if legend.names were given
        if(!is.null(legend.names)){
            legend('bottomright',legend.names,pch=16,cex=.5,col=color.by.color,pt.cex=1)
        }
    }
    
    if(!is.null(filebase)) dev.off()
}


# allows pre-regression on covariates
"plot.association" <- function(gx, qtr, covariates=NULL,
        filebase=NULL,drop.outliers.range=NULL,add.to.plot=FALSE, ymin=NULL){
    require('vioplot')
    if(is.null(dim(gx))) gx <- matrix(gx,ncol=1)
    if(is.null(dim(qtr))) qtr <- matrix(qtr,ncol=1)
    if(is.null(colnames(qtr))) colnames(qtr) <- sprintf('col%d',1:ncol(qtr))
    if(is.null(colnames(gx))) colnames(gx) <- sprintf('col%d',1:ncol(gx))

    if(nrow(gx) != nrow(qtr)){
        stop(sprintf('Error: gx has %d rows, mb has %d rows, covariates has %d rows.\n',
            nrow(gx), nrow(qtr), nrow(covariates)))
    }
    non.na <- rep(TRUE, nrow(qtr))
    if(!is.null(covariates)){
        non.na <- rowSums(is.na(covariates)) == 0
        if(nrow(gx) != nrow(covariates)){
            stop(sprintf('Error: gx has %d rows, mb has %d rows, covariates has %d rows.\n',
                nrow(gx), nrow(qtr), nrow(covariates)))
        }
        mm <- covariate.model.matrix(covariates[non.na,,drop=FALSE])
    }

    if(!is.null(filebase)) pdf(paste(filebase,'.pdf',sep=''), width=8,height=8)        
    
    outlier <- rep(FALSE, nrow(qtr))
    if(!is.null(drop.outliers.range)) outlier <- is.outlier(qtr[,1],range=drop.outliers.range)

    if(!is.null(covariates)){
        y <- lm(qtr[!outlier & non.na,1] ~ mm[!outlier[non.na],,drop=F])$residuals
    } else {
        y <- qtr[!outlier & non.na,1]
    }
    x <- gx[!outlier & non.na,1]
    # if more than 10 unique values, treat as real
    is.discrete <- length(unique(x)) < 10
    ylim <- range(y)
    if(!is.null(ymin)) ylim[1] <- 0
    if(is.discrete){
        x <- as.factor(x)
        stripchart(y ~ x, vert=T, method='jitter',
                    pch=21,bg='#FFFFFF00', col='#FFFFFF00',
                    xlab='', ylab='Residual log relative abundance',
                    cex.main=.75, 
                    cex.lab=.75,ylim=ylim, at=1:length(unique(x)), add=FALSE,
                    xlim=c(0.5,length(unique(x))+.5));
        for(i in 1:length(unique(x))){
            vioplot(y[x==levels(x)[i]],at=i,add=TRUE,col='#00000000',drawRect=FALSE)
        }
        stripchart(y ~ x, vert=T, method='jitter',
                    pch=21,bg='#FFFFFF99', col='#00000099',
                    ylab=sprintf('%s (residual)',1),
                    cex.lab=.75,ylim=ylim, at=1:length(unique(x)), add=TRUE);
    } else {
        plot(x, y, pch=21, bg='#00000044', col='#00000099',
                ylab=sprintf('%s (residual)',1),
                cex.lab=.75,ylim=ylim);
        abline(lm(y ~ x))
        print(cor.test(x,y))
    }
    if(!is.null(filebase)) dev.off()
}

# returns residuals from powertransformed regression
# outliers get NA in return vector
# pass outlier.range=NULL to keep all points
"get.residuals" <- function(x,y,outlier.range=3,do.power.transform=FALSE){
    ix <- !is.outlier(y,outlier.range)
    if(!any(ix)) stop('All points were outliers')
    res <- rep(NA,length(y))
    x <- covariate.model.matrix(x[ix,])
    y <- y[ix]
    if(do.power.transform){
        if(any(y < 0)) stop('y contains negative values, cannot perform power transform')
        if(any(y == 0)) y <- y + min(y[y>0])/2

        require('car')
        p1 <- powerTransform(y ~ x)
        if(p1$convergence > 0){
            warning('Warning, power transform failed; using log')
            # if power transform failed, use pseudolog
            mlm <- lm(pseudolog(y) ~ x)
        } else {
            # fit linear model with transformed response:
            mlm <- lm(bcPower(y, p1$roundlam) ~ x)
        }
    } else {
        mlm <- lm(y ~ x)
    }
    res[ix] <- mlm$residuals
    return(res)
}