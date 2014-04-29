"qq.pvals" <- function(pvals,pvals.null=NULL, filename=NULL,add=FALSE,...){
    observed <- sort(pvals)
    lobs <- -(log10(observed))

    expected <- c(1:length(observed)) 
    if(is.null(pvals.null)){
        lexp <- -(log10(expected / (length(expected)+1)))
    } else {
        lexp <- -log10(pvals.null)
    }
    lims <- c(0,1.05*max(c(lobs,lexp)))
    if(!is.null(filename)) pdf(filename,width=6, height=6)
    if(!add){
        plot(c(0,7), c(0,7), col="red", lwd=3, type="l",
    #         xlab="Expected (-logP)", ylab="Observed (-logP)",
            xlab=expression(Expected~~-log[10](italic(p))),
            ylab=expression(Observed~~-log[10](italic(p))),
            xlim=lims, ylim=lims, las=1, xaxs="i", yaxs="i", bty="l")
    }
    points(lexp, lobs, pch=18, cex=.75, ...) 
    if(!is.null(filename)) dev.off()
}

# returns the probability of seeing as many significant pvalues based on poisson and alpha
"poisson.enrichment" <- function(pvals,nexp=NULL,nobs=NULL, alpha=0.05){
    if(is.null(pvals) && (is.null(nexp) | is.null(nobs))) stop('Please pass a vector of numeric p-values\n')
    if(is.null(nexp)) nexp <- alpha * length(pvals)
    if(is.null(nobs)) nobs <- sum(pvals <= alpha)
    return(ppois(nobs, nexp, lower.tail=FALSE))
}

# linkage diseq., inputs are vectors of (0,1,2)
# g1 and g2 can be matrices, returns g1 x g2 matrix of LDs
"LD.numeric" <- function(g1,g2,verbose=FALSE,ncores=1){
    library('genetics')
    if(is.null(dim(g1))) g1 <- matrix(g1,ncol=1)
    if(is.null(dim(g2))) g2 <- matrix(g2,ncol=1)
    
    # convert each column to string
    s1 <- matrix('A/A',nrow=nrow(g1),ncol=ncol(g1))
    s2 <- matrix('A/A',nrow=nrow(g2),ncol=ncol(g2))
    for(i in 1:ncol(g1)){
        s1[g1[,i] == 1,i] <- 'A/B'
        s1[g1[,i] == 2,i] <- 'B/B'
    }
    for(i in 1:ncol(g2)){
        s2[g2[,i] == 1,i] <- 'A/B'
        s2[g2[,i] == 2,i] <- 'B/B'
    }
    
    res <- matrix(0,nrow=ncol(g1), ncol=ncol(g2))
    rownames(res) <- colnames(g1)
    colnames(res) <- colnames(g2)
    if(ncores==1){
        for(i in 1:ncol(g1)){
            if(verbose) cat(i,'')
            for(j in 1:ncol(g2)){
                if(verbose && j %% 100 == 0) cat('.')
                res[i,j] <- LD(genotype(s1[,i]),genotype(s2[,j]))$`D'`
            }
        }
    } else {
        # parallel
        library('multicore')
        resp <- mclapply(seq(ncol(g1)), function(ixx){
                if(verbose) cat(ixx,'')
                res.i <- numeric(ncol(g2))
                for(j in 1:ncol(g2)){
                    if(verbose && j %% 10 == 0) cat('.')
                    res.i[j] <- LD(genotype(s1[,ixx]),genotype(s2[,j]))$`D'`
                }
                cat('i =',ixx,'max ld =',max(res.i[j]), '| ')
                return(res.i)
            },
            mc.cores=ncores
        )
        for(i in 1:ncol(g1)){
            res[i,] <- resp[[i]]
        }
    }

    if(verbose) cat('\n')
#     g1c <- rep('A/A',length(g1)); g1c[g1==1] <- 'A/B'; g1c[g1==2] <- 'B/B'
#     g2c <- rep('A/A',length(g2)); g2c[g2==1] <- 'A/B'; g2c[g2==2] <- 'B/B'
#     ld <- try(LD(genotype(g1c),genotype(g2c)))
#     if(class(ld) == 'try-error') browser()
#     return(ld$`D'`)
    return(res[,])
}

# return association p-values for all pairs of columns
# can be data frame/matrix
"getPmat" <- function(X, adjust=TRUE){
	pvals <- matrix(0, ncol(X), ncol(X))
	colnames(pvals) <- colnames(X)
	rownames(pvals) <- colnames(X)
	for(i in 1:(ncol(X) - 1)){
		for(j in (i+1):ncol(X)){
			pvals[i,j] <- getP(X[,i], X[,j])
			pvals[j,i] <- pvals[i,j]
		}
	}
	if(adjust){
		qvals <- p.adjust(pvals[upper.tri(pvals)],'fdr')
		pvals[upper.tri(pvals)] <- qvals
		pvals[lower.tri(pvals)] <- t(pvals)[lower.tri(pvals)]
	}
	return(pvals)
}

# return association p-value for any pair of vectors
# each can be numeric or factor/character
"getP" <- function(x,y,parametric=FALSE){
 	drop.ix <- is.na(x) | is.na(y)
 	x <- x[!drop.ix]
 	y <- y[!drop.ix]
	if(!is.numeric(x)) x <- as.factor(x)
	if(!is.numeric(y)){
		# ensure that if only one is numeric, it's x
		z <- as.factor(y)
		y <- x
		x <- z
	}
		
 	if(parametric){
 	
 	} else {
 		if(is.numeric(x)){
 			if(is.numeric(y)){
 				# both numeric
 				return(cor.test(x,y,method='spear',exact=FALSE)$p.value)
 			} else {
 				# x numeric, y not
 				if(length(levels(y)) == 2){
 					ix <- y==levels(y)[1]
 					return(wilcox.tests(x[ix], x[!ix])$p.value)
				} else {
	 				return(kruskal.test(x,y)$p.value)
	 			}
 			}
 		} else {
 			# both non-numeric, do chi-sq
 			chisq.test(table(x,y))$p.value
 		}
 	}
 }



#
# return list of lists of the missense/frameshift/stop-gained variants
# for each gene
# requires an ichip map with columns "#ichipID", "geneID", "fxn_class"
"pathogenic.ichipIDs.by.gene" <- function(ichipmap,
         fxn.classes = c("missense","frameshift-variant", "stop-gained","intron-variant","nc-transcript-variant")
        ){
        
    if(is.null(fxn.classes)) fxn.classes <- c("missense","frameshift-variant", "stop-gained","intron-variant","nc-transcript-variant")
    # drop NA's
    non.na.ix <- !is.na(ichipmap[,'geneID']) & ichipmap[,'geneID'] != '?'
    ichipmap <- ichipmap[non.na.ix,]
    
    fxn.pattern <- paste(fxn.classes,collapse='|',sep='')
    patho.ix <- grep(fxn.pattern, ichipmap[,'fxn_class'])
    if(fxn.classes[1] == 'all') patho.ix <- 1:nrow(ichipmap)
    ichipmap <- ichipmap[patho.ix,]
    if(class(ichipmap) == 'data.frame') ichipmap <- droplevels(ichipmap)
    # split ichipIDs by geneID
    res <- split(ichipmap[,'#ichipID'],ichipmap[,'geneID'])
    return(res)
}


# convenience wrapper on matrix eqtl
# snps is a matrix of sample x snp counts
# variates is a matrix of sample x variate (e.g. taxa)
# covariates is a matrix of sample x covariate (real or factors OK)
# if full.result, includes the Matrix_eQTL_engine result object (for plotting, e.g.)
"eqtl" <- function(snps, variates, covariates=NULL, max.p=1e-10,
            full.result=FALSE, verbose=TRUE, filename=tempfile()){
    library('MatrixEQTL')
    snps1 = SlicedData$new( t( snps ) );
    gene1 = SlicedData$new( t( variates ) );
    if(is.null(covariates)){
        cvrt1 <- SlicedData$new()
    } else {
        # get dummy variables, drop intercept
        mm <- model.matrix(~ ., covariates)[,-1]
        cvrt1 = SlicedData$new( t( mm ) );
    }
    snps1$ResliceCombined(500);
    gene1$ResliceCombined(500);
    
    meh = Matrix_eQTL_engine(snps = snps1, gene = gene1, cvrt = cvrt1, 
        output_file_name = filename,
        pvOutputThreshold = max.p, useModel = modelLINEAR,
        errorCovariance = numeric(), verbose = verbose, pvalue.hist = 'qqplot');
    res <- data.frame(meh$all$eqtls[c('snps','gene','FDR')])
    res[,1] <- as.character(res[,1])
    res[,2] <- as.character(res[,2])
    if(full.result){
        return(list(res=res,me=meh))
    } else {
        return(res)
    }
}


# returns snp counts collapsed by gene
"collapse.ichip.by.gene" <- function(counts, ichipmap, verbose=TRUE,
        fxn.classes = NULL){
    # Collapse by function
    bad.snps <- pathogenic.ichipIDs.by.gene(ichipmap,fxn.classes=fxn.classes)
    
    # Collapse by gene
    cat('aggregating genes\n')
    counts.by.gene <- NULL
    for(i in 1:length(bad.snps)){
        if(verbose) if(i %% 100 == 0) cat(sprintf('%.2f ',i/length(bad.snps)))
        gene.name <- names(bad.snps)[i]
        common <- intersect(bad.snps[[i]],colnames(counts))
        if(length(common) > 0){
            non.na.ix <- apply(counts[,common,drop=FALSE],2,function(xx) sum(is.na(xx)) == 0)
            if(sum(non.na.ix) > 0){
                newx <- rowSums(counts[,common[non.na.ix],drop=FALSE])
                counts.by.gene <- cbind(counts.by.gene,newx)
                colnames(counts.by.gene)[ncol(counts.by.gene)] <- gene.name
            }
        }
    }
    if(verbose) cat('\n')
    return(counts.by.gene)
}

# returns vector of cluster ids for clusters with internal
# complete-linkage correlation of min.cor
"cluster.by.correlation" <- function(x, min.cor=.75){
#     library('fastcluster')
    cc <- cor(x,use='pairwise.complete.obs',method='pear')
    if(ncol(x) == 379) browser()
    cc <- as.dist(1-cc)
    hc <- hclust(cc)
    res <- cutree(hc,h=1-min.cor)
    names(res) <- colnames(x)
    return(res)
}

# returns vector of cluster ids for clusters with internal
# complete-linkage correlation of min.cor
#
# by default, chooses cluster reps as highest-variance member
# if select.rep.fcn=mean
"collapse.by.correlation" <- function(x, min.cor=.5, select.rep.fcn=c('var','mean','lowest.mean',
			'longest.name', 'shortest.name')[4],
        verbose=FALSE){
    if(verbose) cat('Clustering',ncol(x),'features...')
    gr <- cluster.by.correlation(x, min.cor=min.cor)
    if(verbose) cat('getting means...')
    if(select.rep.fcn == 'mean'){
        v <- apply(x,2,function(xx) mean(xx,na.rm=TRUE))
    } else if(select.rep.fcn == 'lowest.mean'){
        v <- apply(x,2,function(xx) -mean(xx,na.rm=TRUE))
    } else if(select.rep.fcn == 'longest.name'){
        v <- nchar(colnames(x))
    } else if(select.rep.fcn == 'shortest.name'){
        v <- -nchar(colnames(x))
    } else {
        v <- apply(x,2,function(xx) var(xx,use='complete.obs'))
    }
    if(verbose) cat('choosing reps...')
    reps <- sapply(split(1:ncol(x),gr),function(xx) xx[which.max(v[xx])])
    if(verbose)
        cat(sprintf('collapsed from %d to %d.\n',ncol(x), length(reps)))
    return(list(reps=reps, groups=gr))
}

# returns vector of cluster ids for clusters with internal
# complete-linkage correlation of min.cor
#
# by default, chooses cluster reps as highest-variance member
# if select.rep.fcn=mean
"collapse.by.correlation.recursive" <- function(x, min.cor=.5, select.rep.fcn=c('var','mean')[2],
        verbose=FALSE,max.n.per.step=1000){
    N <- ncol(x)
    if(N > max.n.per.step){
        if(verbose) cat('Clustering',ncol(x),'features...')
        # divide into bins of max.n.per.step
        nsteps <- ceiling(N / max.n.per.step)
        nsteps <- ceiling(N / max.n.per.step)
    } else {
        if(verbose) cat('Clustering',ncol(x),'features...')    
    }
    gr <- cluster.by.correlation(x, min.cor=min.cor)
    if(verbose) cat('getting means...')
    if(select.rep.fcn == 'mean'){
        v <- apply(x,2,function(xx) mean(xx,na.rm=TRUE))
    } else {
        v <- apply(x,2,function(xx) var(xx,use='complete.obs'))
    }
    if(verbose) cat('choosing reps...')
    reps <- sapply(split(1:ncol(x),gr),function(xx) xx[which.max(v[xx])])
    if(verbose)
        cat(sprintf('collapsed from %d to %d.\n',ncol(x), length(reps)))
    return(list(reps=reps, groups=gr))
}

# collapse.at = .6 means collapse any features sharing > 60% of eachother's bits
"covariate.model.matrix" <- function(covariates,collapse.at=NULL){
  if(is.null(covariates)) return(NULL)
  # drop covariates that are constant
  covariates <- covariates[,apply(covariates,2,function(xx) length(unique(xx))) > 1,drop=F]
  
  if(is.null(dim(covariates))) stop('No non-constant covariates\n')
  
  covariates <- as.data.frame(covariates)
  covariates <- droplevels(covariates)
  
  covariate.names <- colnames(covariates)
  non.numeric <- sapply(covariates, class) != 'numeric'
  mm <- NULL
  
  if(any(non.numeric)){
    for(i in which(non.numeric)){
      covariates[,i] <- as.factor(covariates[,i])
    }
    mm <- model.matrix(~ ., covariates[,non.numeric])[,-1]
  }
  if(any(!non.numeric)){
    mm <- cbind(mm, as.matrix(covariates[,!non.numeric]))
  }
  if(is.null(dim(mm))) mm <- matrix(mm,nc=1)
  if(ncol(mm) == length(covariate.names)) colnames(mm) <- covariate.names
  mm <- mm[,apply(mm,2,var) > 0,drop=F]
  
  if(!is.null(collapse.at)){
    mm <- collapse.by.mutual.information(mm,threshold=collapse.at)
  }
  return(mm)
}

# collapses features if they explain more than 50% of eachothers bits
# if subset not NULL, uses subset for determination of mutual information, but still
# returns full model matrix
# nbins is the number of bins for collapsing continuous variables (default sqrt(nrow(x)))
"collapse.by.mutual.information" <- function(x,threshold=.5,
        plot.filename=NULL, do.plot=FALSE,
        table.filename=NULL,
        subset.ix=NULL, nbins=NULL,
        return.muc=FALSE){
    library('infotheo')
    if(is.null(subset.ix)) subset.ix <- 1:nrow(x)
    if(is.null(nbins)) nbins <- sqrt(length(subset.ix))
    mm <- covariate.model.matrix(x[subset.ix,])
    
    mmd <- mm;
    mmd[,'Years_since_diagnosis'] <- discretize(mmd[,'Years_since_diagnosis'],nbins=nbins)[,1]
    mmd[,'Age'] <- discretize(mmd[,'Age'],nbins=nbins)[,1]
    mmdf <- as.data.frame(mmd)
    mi <- mutinformation(mmdf)

    # uncertainty coefficient
    uc <- sweep(mi,1,apply(mmdf,2,entropy),'/')

    # maximum uncertainty coefficient
    # 1 - muc is maximum fraction of bits of one not explained by other
    muc <- uc
    muc[,] <- ifelse(uc < t(uc),uc,t(uc))

    hc <- hclust(as.dist(1-muc))
    groups <- cutree(hc,h=threshold)
    group.list <- split(names(groups), groups)
    keep.names <- NULL
    for(i in seq_along(group.list)){
        if(length(group.list[[i]]) > 2){
            err.msg <- sprintf('Collapsing groups of size > 2 not implemented. Cannot collapse group containing [%s].\n',
                    paste(group.list[[i]],collapse=', '))
            stop(err.msg)
        } else if(length(group.list[[i]]) == 2){
            group.names <- 
            if(uc[group.list[[i]][1],group.list[[i]][2]] > uc[group.list[[i]][2],group.list[[i]][1]]){
                keep.names <- c(keep.names, group.list[[i]][1])
            } else {
                keep.names <- c(keep.names, group.list[[i]][2])
            }
        } else {
            keep.names <- c(keep.names,group.list[[i]])
        }
    }
    
    # do plot
    if(do.plot || !is.null(plot.filename)){
        if(!is.null(plot.filename)) pdf(plot.filename,width=5,height=7)
        par(cex=.6)
        plot(hc,xlab='',ylab='Symmetric uncertainty coefficient',lty=1,lwd=1,sub='',axes=FALSE)
        axis(side=2, at=seq(.4,1,.1), col="#F38630",labels=FALSE, lwd=2)
        # add text in margin
        mtext(seq(.4,1,.1), side=2, at=seq(.4,1,.1),
              line=1, col="#A38630", las=2,cex=.7)
        if(!is.null(plot.filename)) dev.off()
    }
        
    # write association table
    if(!is.null(table.filename)){
        sink(table.filename)
        cat('\t')
        write.table(1-uc,sep='\t',quote=F)
        sink(NULL)
    }
    
    # get full model matrix in case a subset was used
    mm <- covariate.model.matrix(droplevels(as.data.frame(x)))
    if(is.data.frame(x)) mm <- as.data.frame(mm)

    if(return.muc) return(muc)
    return(mm[,keep.names])
}



"fishersMethod" <- function(pvals) pchisq(-2 * sum(log(pvals)),df=2*length(pvals),lower=FALSE)


# returns table of gx.feature name, qt.feature name, pval1 pval2 pval3 ... fisher.pval fisher.qval
# takes only top 1000 test in each result table
# note: just returns 0 for all
"combine.results" <- function(results, gx.features, qt.features){
    if(any(sapply(results, function(xx) any(drop.na(xx) == 0)))) stop('Zero values found in results list')
   
    N <- length(gx.features)
    M <- length(qt.features)
    combined.results <- as.data.frame(matrix(NA,nrow=N * M,ncol=2+length(results)))
    colnames(combined.results) <- c('gx','qtrait',names(results))
    combined.results[,1] <- rep(gx.features,each=M)
    combined.results[,2] <- rep(qt.features,N)
#     cat('Merging results into one data frame...\n')
#     for(i in 1:length(results)){
#         for(j in 1:nrow(results[[i]])){
#             gx.id <- results[[i]][j,1]
#             qt.id <- results[[i]][j,2]
#             row.ix <- combined.results[,1] == gx.id & combined.results[,2] == qt.id
#             if(sum(row.ix) > 1) stop(paste(gx.id, qt.id, 'found in multiple rows.\n'))
#             combined.results[row.ix,i+2] <- as.numeric(results[[i]][j,'pvalue'])
#         }
#     }
    
    cat('Merging results...\n')
    rownames(combined.results) <- paste(combined.results[,1], combined.results[,2], sep='__')
    for(i in 1:length(results)){
        rownames.i <- paste(results[[i]][,1], results[[i]][,2],sep='__')
        combined.results[rownames.i,i+2] <- results[[i]][,'pvalue']
    }

    # extract top n in each class
    results.topn <- order(combined.results[,3])
    for(i in 4:ncol(combined.results)){
        results.topn <- intersect(results.topn,order(combined.results[,i]))
    }
    ntests <- nrow(combined.results)
    cat(length(results.topn),'out of',ntests,'tests\n')
    combined.results <- combined.results[results.topn,]
    combined.results <- combined.results[rowSums(is.na(combined.results)) == 0,]
    combined.pvals <- apply(combined.results[,-(1:2)], 1, fishersMethod)
    combined.qvals <- p.adjust(combined.pvals, 'fdr', n=ntests)
    combined.results <- cbind(combined.results, combined.pvals, combined.qvals)
    colnames(combined.results)[1:2 + (ncol(combined.results) - 2)] <- c('pvalue','qvalue')
    combined.results <- combined.results[order(combined.results[,'pvalue']),]
    combined.results <- combined.results[order(combined.results[,'qvalue']),]
    return(combined.results)
}


"my.error.bars" <- function(x,centers,spread,barw.pct=NULL,...){
    
	xlim <- range(x)
    if(is.null(barw.pct)){
        width = min(.010,.25/length(x))
    	barw <- diff(xlim) * width
    } else {
        barw <- diff(xlim) * barw.pct 
    }

    upper <- centers + spread
    lower <- centers - spread

	segments(x, upper, x, lower, ...)
	segments(x - barw, upper, x + barw, upper, ...)
	segments(x - barw, lower, x + barw, lower, ...)
}



"scatterhist" <- function(x, y, xlab="", ylab=""){
  zones=matrix(c(2,0,1,3), ncol=2, byrow=TRUE)
  layout(zones, widths=c(4/5,1/5), heights=c(1/5,4/5))
  xhist = hist(x, plot=FALSE)
  yhist = hist(y, plot=FALSE)
  top = max(c(xhist$counts, yhist$counts))
  par(mar=c(3,3,1,1))
  plot(x,y)
  par(mar=c(0,3,1,1))
  barplot(xhist$counts, axes=FALSE, ylim=c(0, top), space=0)
  par(mar=c(3,0,1,1))
  barplot(yhist$counts, axes=FALSE, xlim=c(0, top), space=0, horiz=TRUE)
  par(oma=c(3,3,0,0))
  mtext(xlab, side=1, line=1, outer=TRUE, adj=0, 
    at=.8 * (mean(x) - min(x))/(max(x)-min(x)))
  mtext(ylab, side=2, line=1, outer=TRUE, adj=0, 
    at=(.8 * (mean(y) - min(y))/(max(y) - min(y))))
}

# greedily drops rows and columns to remove NA's 
# from a matrix with as few drops as possible
# returns a list of boolean row.drop.ix, col.drop.ix
# set prefer.columns to 2 means columns must have twice the NA's as rows to be
# removed first
# stop.at = stop when this rate na's or less on average
"drop.na.2d" <- function(x,prefer.columns=1,max.row.na.rate=.05,max.col.na.rate=.05){
    drop.row <- rep(FALSE,nrow(x))
    drop.col <- rep(FALSE,ncol(x))
    
    na.x <- is.na(x[!drop.row, !drop.col,drop=F])
    while(ncol(na.x) > 1 && nrow(na.x) > 1 
            && (max(colMeans(na.x)) >= max.col.na.rate || max(rowMeans(na.x)) >= max.row.na.rate)){
        row.na.sums <- rowSums(na.x > 0)
        col.na.sums <- colSums(na.x > 0)
        row.max <- max(row.na.sums)
        col.max <- max(col.na.sums)
        if(row.max >= col.max/prefer.columns && ncol(na.x) > 2){
            drop.row[which(!drop.row)[which.max(row.na.sums)]] <- TRUE
        } else {
            drop.col[which(!drop.col)[which.max(col.na.sums)]] <- TRUE
        }
        na.x <- is.na(x[!drop.row, !drop.col,drop=F])
    }
    return(list(drop.row=drop.row, drop.col=drop.col))
}


# if zero.replace is NULL, uses half smallest nonzero value
"power.transform" <- function(qtr, zero.replace=NULL){
    require('MASS')
    require('car')
    if(is.null(dim(qtr))) qtr <- matrix(qtr,ncol=1)
    colnames.save <- colnames(qtr)
    if(any(qtr == 0)){
        if(is.null(zero.replace)){
            qtr <- qtr + min(qtr[qtr>0])/2
        } else {
            qtr[qtr==0] <- zero.replace
        }
    }
    res <- apply(qtr,2,function(xx) boxcox(xx ~ 1,plotit=FALSE))
    lambda <- sapply(res,function(xx) xx$x[which.max(xx$y)])
    qtr <- bcPower(qtr, lambda)
    colnames(qtr) <- colnames.save
    return(qtr)
}

# coef of determination
"r.squared" <- function(y,yhat){
    sse <- sum( (y -    yhat)**2 )
    sst <- sum( (y - mean(y))**2 )
    return(max(0,1 - sse/sst))
}

# convenience wrapper for p.distribution.test that takes association tables
# instead of p-value matrices
"p.distribution.test.from.association.tables" <- function(obs,null=NULL, snp.names=NULL, pvalue.column='pvalue',
        rs2gene=NULL, ...){
    if(!is.null(snp.names)) obs[,1] <- snp.names
    p.sets <- t(sapply(split(obs[,pvalue.column], obs[,1]),function(xx) xx))
    if(!is.null(rs2gene)) rownames(p.sets) <- rs2gene[rownames(p.sets)]
    null.sets <- NULL
    if(!is.null(null)) null.sets <- t(sapply(split(null[,pvalue.column], null[,1]),function(xx) xx))
    return(p.distribution.test(p.sets, null.sets, ...))
}

# takes a matrix of p-values
# each row is a different "set" of p values
# we want to find sets that are outliers
# optionally takes a separate set of sets to form
# the null distribution of sets
# null.sets and p.sets must have same # columns
# note: alpha is for one-tailed test
# returns indices of hits
# 
# hits.to.show is an optional list of rownames of p.sets that
# are permitted to be displayed as hits
# legend.order determines order of genes for legend text
"p.distribution.test" <- function(p.sets, null.sets=NULL, 
        alpha=.05,do.plot=TRUE, plot.alpha=0.05,
        plot.null=FALSE, plot.nonsig=TRUE,
        n.plot=NULL,
        hits.to.show=NULL,
        test.type=c('average','slope','ratio','diff','geom','geom50','harm','nsignif','rankproduct','qvalue-p0','kstest-parametric','kstest')[1],
        filename=NULL,
        use.qvalue=FALSE,
        important.hits.only=FALSE
    ){
    if(is.null(null.sets)){
        null.sets <- p.sets
    }
    if(any(is.na(p.sets))) {
        cat(sprintf('Warning: NAs in observed p.sets, setting NAs to 1 in %d rows\n',sum(rowSums(is.na(p.sets)) > 0)))
        p.sets[is.na(p.sets)] <- 1
    }
    if(any(is.na(null.sets))) {
        cat(sprintf('Warning: NAs in observed null.sets, setting NAs to 1 in %d rows\n',sum(rowSums(is.na(null.sets)) > 0)))
        null.sets[is.na(null.sets)] <- 1
    }
    
    if(important.hits.only){
        hits.to.show <- c(
                            'MTPAP',
                            'HINT1',
                            'HCK',
                            'PCNXL3',
                            'IL18RAP',
                            'NOD2',
                            'RIMBP3C',
                            'RIPK2',
                            'TMEM189-UBE2V1',
                            'STAT3',
                            'JAK2',
                            'CARD9',
                            'C21orf54',
                            'GLS',
                            'TNFSF15'
                        )                          
    }
    if(ncol(null.sets) != ncol(p.sets)) stop('Null sets and observed sets must have same number of p-values.\n')
    n <- ncol(p.sets)
    # x is the expected list of pvals
    x <- sort(-log10((1:n)/(1+n)))
    x.pvals <- sort((1:n)/(1+n))
    p.sets.pvals <- p.sets
    null.sets.pvals <- null.sets
    for(i in 1:nrow(p.sets)) p.sets[i,] <- sort(-log10(p.sets[i,]))
    for(i in 1:nrow(null.sets)) null.sets[i,] <- sort(-log10(null.sets[i,]))
    qlow <- apply(null.sets,2,quantile,alpha)
    qhigh <- apply(null.sets,2,quantile,1-alpha)

    if(test.type == 'average'){
        pvals <- apply(p.sets,1,function(xx) mean(colMeans(sweep(null.sets,2,xx,'>='))))
    } else if (test.type == 'ratio'){
        ratios <- sweep(p.sets.pvals, 2, x.pvals, '/')
        null.ratios <- sweep(null.sets.pvals, 2, x.pvals, '/')
        pvals <- apply(ratios,1,function(xx) mean(colMeans(sweep(null.ratios,2,xx,'<='))))
    } else if (test.type == 'geom'){
        geoms <- apply(p.sets.pvals,1,function(xx) exp(mean(sum(log(xx)))))
        null.geoms <- apply(null.sets.pvals,1,function(xx) exp(mean(sum(log(xx)))))
        pvals <- sapply(geoms,function(xx) mean(c(null.geoms,xx) <= xx))
    } else if (test.type == 'harm'){
        harms <- apply(p.sets.pvals,1,function(xx) 1/mean(1/xx))
        null.harms <- apply(null.sets.pvals,1,function(xx) 1/mean(1/xx))
        pvals <- sapply(harms,function(xx) mean(c(null.harms,xx) <= xx))
    } else if (test.type == 'geom50'){
        geoms <- apply(p.sets.pvals,1,function(xx) exp(mean(sum(log(xx[xx<median(xx)])))))
        null.geoms <- apply(null.sets.pvals,1,function(xx) exp(mean(sum(log(xx[xx<median(xx)])))))
        pvals <- sapply(geoms,function(xx) mean(c(null.geoms,xx) <= xx))
    } else if (test.type == 'slope-intercept' || test.type == 'slope'){
        if(test.type == 'slope-intercept') {
            slopes <- apply(p.sets,1,function(xx) coef(lm(xx ~ x+1))['x'])
        } else {
            slopes <- apply(p.sets,1,function(xx) coef(lm(xx ~ x+0))['x'])
        }
        #  get slopes for every increment of 1/(number of null sets)
        # or 1/999, whichever is smaller
        nsteps <- min(999, nrow(null.sets))
        null.slopes <- numeric(nsteps)
        null.pvals <- numeric(nsteps)
        for(i in 1:nsteps){
            pval <- i / (nsteps + 1)
            qhigh <- apply(null.sets,2,quantile,1-pval)
            if(test.type == 'slope-intercept') {
                null.slope.i <- coef(lm(qhigh ~ x+1))['x']
            } else {
                null.slope.i <- coef(lm(qhigh ~ x+0))['x']
            }
            null.slopes <- c(null.slopes, null.slope.i)
            null.pvals <- c(null.pvals, pval)
        }
        pvals <- 1-sapply(slopes, function(xx) null.pvals[which.max(which(null.slopes <= xx))])
    } else if (test.type == 'nsignif') {
        counts <- rowSums(p.sets.pvals <= alpha)
        null.counts <- rowSums(null.sets.pvals <= alpha)
        pvals <- sapply(counts,function(xx) mean(c(null.counts,xx) >= xx))
    } else if (test.type == 'kstest-parametric') {
        pvals <- apply(p.sets.pvals,1,function(xx) ks.test(xx,y='punif',alternative='greater')$p.value)
    } else if (test.type == 'kstest') {
        pvals <- apply(p.sets.pvals,1,function(xx) ks.test(xx,y='punif',alternative='greater')$statistic)
        null.pvals <- apply(null.sets.pvals,1,function(xx) ks.test(xx,y='punif',alternative='greater')$statistic)
        pvals <- sapply(pvals,function(xx) mean(c(null.pvals,xx) >= xx))
    } else if (test.type == 'qvalue-p0'){
    	require('qvalue')
    	p0 <- apply(p.sets.pvals,1,function(xx) qvalue(xx)$pi0)
    	null.p0 <- apply(null.sets.pvals,1,function(xx) qvalue(xx)$pi0)
        pvals <- sapply(p0,function(xx) mean(c(null.p0,xx) >= xx))
    } else if (test.type == 'rankproduct'){
		ranks <- rank(c(as.numeric(t(p.sets.pvals)), as.numeric(t(null.sets.pvals))))
		memberships <- rep(1:nrow(p.sets.pvals),each=ncol(p.sets.pvals))
		null.memberships <- rep(1:nrow(null.sets.pvals),each=ncol(p.sets.pvals))
		rank.prods <- sapply(split(ranks[1:length(memberships)], memberships), 
				function(xx) sum(log(xx)))
		null.rank.prods <- sapply(split(ranks[length(memberships) + seq(length(null.memberships))], null.memberships), 
				function(xx) sum(log(xx)))
		pvals <- sapply(rank.prods,function(xx) mean(c(null.rank.prods,xx) <= xx))
		names(pvals) <- rownames(p.sets)
    } else {
        stop('Unknown test type\n')
    }

    signif <- pvals <= plot.alpha
    if(!is.null(hits.to.show)) signif <- signif & names(pvals) %in% hits.to.show
    sigix <- which(signif)[order(pvals[signif])]
    if(!is.null(n.plot)) if(length(sigix) > n.plot) sigix <- sigix[1:n.plot]
    # hack to put nod2 on top
#     if('NOD2' %in% rownames(p.sets)[sigix]){
#         nod2ix <- which(rownames(p.sets)[sigix] == 'NOD2')
#         sigix <- c(sigix[nod2ix], sigix[-nod2ix])
#     }
#     cat('95% upper confidence interval:',boundary.slope,'\n')
#     browser()
    if(do.plot){
        if(!is.null(filename)) {
            pdf(filename,width=4,height=4)
            par(mar=c(5,5,1.5,1.5));
        }
        sig.col <- sprintf('%sbb',c(bong.cols,brewer.pal(9,'Set1'),brewer.pal(12,'Set3')))
        sig.col <- sprintf('%sbb',c(brewer.pal(9,'Set1')[-6],brewer.pal(12,'Set3')))
        sig.col <- c(sprintf('%sbb',bong.cols[]), sig.col)
#         sig.col <- c(sprintf('%sff',bong.cols[-1]), sig.col)
        xlim <- range(x)
        if(length(sigix) > 0){
            ylim <- range(p.sets[sigix,])
        } else {
            ylim <- range(p.sets)
        }
#         xlim <- c(0,2.3)
#         ylim <- c(0,4.5)
        plot(x,null.sets[1,],xlim=xlim,ylim=ylim,type='n',
            xlab='Expected -log10(p)',ylab='Observed -log10(p)',
            cex.axis=.75, cex.lab=.85)
        nonsig.col <- '#00000004'
        if(plot.null){
            # plot at most 200 of the null set, otherwise it's too noisy
            for(i in sample(nrow(null.sets),min(400,nrow(null.sets)))){
                # null
                points(x,null.sets[i,],col=nonsig.col,pch=16,cex=.9)
            }
            devs <- apply(null.sets,2,sd)
            y <- colMeans(null.sets)
#             devs <- apply(null.sets,2,mad) * 1.4826
#             y <- apply(null.sets,2,median)
#             polycols <- c('#999999','#cccccc','#eeeeee')
#             for(i in length(polycols):1){
#                 yi <- c(y + i * devs, rev(y - i * devs))
#                 polygon(c(x,rev(x)), ifelse(yi > 0,yi,0), col=polycols[i],border=NA) 
#             }
            
        }
        # plot nonsignificant observed sets if requested
        if(any(!signif)){
            if(plot.nonsig){
                
                for(i in 1:sum(!signif)){
                    ix <- which(!signif)[i]
                    points(x,p.sets[ix,],
                        col=nonsig.col,pch=16,cex=.9)
                }
            }
        }
        
        # intervals
        if(test.type=='slopes' || test.type == 'average'){
            int.col <- '#ff2244'
            lty=2
            lines(x,qhigh,col=int.col,lwd=2,lty=lty)
            lines(x,qlow,col=int.col,lwd=2,lty=lty)
            abline(0,1,col='black',lwd=2,lty=lty)
        }
        # hits
        if(length(sigix) > 0){
            for(i in length(sigix):1){
                y <- p.sets[sigix[i],]

                # print slope of best fit for this hit
#                 cat(rownames(p.sets)[sigix[i]],':',pvals[sigix[i]],'\n')
                lines(smooth.spline(x,y, df=5),
                    bg=sig.col[(i - 1) %% length(sig.col) + 1],
                    col=sig.col[(i - 1) %% length(sig.col) + 1],
                    pch=16,
                    lwd=4,type='l',cex=2)
#                 points(x,y,
#                     bg=sig.col[(i - 1) %% length(sig.col) + 1],
#                     col='black',
#                      pch=21, 
# #                     col=sig.col[(i - 1) %% length(sig.col) + 1],
# #                      pch=16, 
#                     lwd=1,type='p',cex=1)
# #                 lines(x,predict(loess(y ~ x)),
# #                     bg=sig.col[(i - 1) %% length(sig.col) + 1],
# #                     col=sig.col[(i - 1) %% length(sig.col) + 1],
# #                      pch=16,
# #                     lwd=2,type='l',cex=1)
            }
            if(!is.null(rownames(p.sets))){
#                 legend('topleft',c(rownames(p.sets)[sigix],'y = x'),col=c(rep('black',length(sigix)),'red'),
#                     pt.bg=c(sig.col[(1:length(sigix) - 1) %% length(sig.col) + 1],NA),
#                     lty=c(rep(NA,length(sigix)),2),lwd=c(rep(NA,length(sigix)),1.5),cex=.7,pt.cex=1,pch=c(rep(21,length(sigix)),NA))
                legend('topleft',c(rownames(p.sets)[sigix],'y = x'),col=c(sig.col[(1:length(sigix) - 1) %% length(sig.col) + 1],'red'),
                    lty=c(rep(1,length(sigix)),2),lwd=c(rep(3,length(sigix)),2),cex=.7,pt.cex=1,pch=c(rep(NA,length(sigix)),NA))
            }
        }
#         lines(x,x,col='red',lwd=3,lty=2)
#         abline(0,1,col='red',lwd=3,lty=2)
        yy <- as.numeric(null.sets)
        xx <- rep(x,each=nrow(null.sets))
        abline(lm(yy ~ xx),col='red',lwd=3,lty=2)
        if(!is.null(filename)) dev.off()
    }
    if(use.qvalue){
        require('qvalue')
        qvals <- try(qvalue(pvals)$qvalue,silent=TRUE)
    }
    if(!use.qvalue || class(qvals)=='try-error'){
		qvals <- p.adjust(pvals,method='fdr')
    }

    return(list(
        hits=names(sigix),
        pvals=sort(pvals),
        qvals=qvals[names(sort(pvals))],
        alpha=alpha))
}



# note: this method takes a data list containing map, covariates, gx, mb, mbf, dx, etc.
# inflamed=na means include all
# cohort = list('prism'=ix1, 'msh'=ix2) means do stratified + fisher's method
# cohort = na means include all
# cohort = logical for inclusion, must match dim(map)
# min.years = 2.0001 means only include patients diagnosed > 2 years ago
# note: test.type='chisq' implies mbf.transform and mb.transform = 'presence-absence'
# if filebase = NULL, uses predictor_test.type
# include.gx.pc=3 means include the first 3 PC's of genetics from d$gxx (alt: NULL)
"test.wrapper" <- function(d, predictor=c('NOD2','163+NOD2','IChip','null','null-500','Inflamed')[1],
                           targets=c('module', 'mb', 'alpha', 'beta','pathway3',),
                           test.type=c('linear','linear-nonzero','np','chisq','perm.x','perm.y','perm.xy','negbin','poisson')[1],
                           mbf.transform=c('asin-sqrt','pseudolog','power','none','presence-absence','split-median')[3],
						   mb.transform=c('asin-sqrt','pseudolog','power','none','presence-absence','split-median')[1],
						   min.mb.prevalence=NULL,
						   max.mb.prevalence=NULL, # can be useful for logistic model
						   min.mb.average=NULL,
                           cohorts=NA,
                           min.age=18, max.age=75, min.years=0,
                           inflamed=NA,
                           np.test.against=0,
                           nperm=999,
                           target.names=NA,
                           alpha=.2,
                           drop.outliers=3,
                           n.bdiv.axes=3,
                           min.maf=.2, # e.g. min.maf=0.2 means only keep snps with > 20% minor allele freq
                           verbose=TRUE,
#                            print.results=TRUE,
                            print.results=FALSE,
                           plot.results=FALSE,
                           export.results=FALSE,
                           linear.power.transform=FALSE,
                           filebase=NULL, 
                           covariates.override=NULL,
                           prefilter.housekeeping.pathways=FALSE,
                           prefilter.outlier.pathways=FALSE,
                           use.qvalue=TRUE,
                           max.pct.zero=NULL,
                           include.gx.pc=TRUE,
                           top.n.results=NULL,
                           ncores=1
                    ) with(d,{
                     
    results <- list()
    attr(results,'call') <- sys.call(1)
	ignore.pathways <- c("Transcription machinery","Translation proteins","DNA replication","Function unknown","General function prediction only","Homologous recombination","Membrane and intracellular structural molecules","Transcription machinery","Translation proteins","Base excision repair","Cell cycle - Caulobacter","Chaperones and folding catalysts","Chromosome","DNA repair and recombination proteins","DNA replication proteins","Mismatch repair","Others","Plant-pathogen interaction","Replication, recombination and repair proteins","Ribosome Biogenesis","RNA polymerase","Transcription factors","Translation factors",'Ribosome',"Protein folding and associated processing","RNA degradation","Drug metabolism - other enzymes","Drug metabolism - other enzymes","Nucleotide excision repair","Photosynthesis proteins", "Carbon fixation in photosynthetic organisms","Photosynthesis proteins")
	outlier.pathways <- c('Methane metabolism','Glycolysis / Gluconeogenesis','Phosphotransferase system (PTS)',    'Porphyrin and chlorophyll metabolism',    'Selenocompound metabolism',    'Vitamin B6 metabolism',    'Fatty acid biosynthesis',    'Cysteine and methionine metabolism',    'Riboflavin metabolism',    'ABC transporters',    'Ubiquinone and other terpenoid-quinone biosynthesis',    'Thiamine metabolism',    'Aminoacyl-tRNA biosynthesis',    'Ribosome',    'Amino sugar and nucleotide sugar metabolism',    'Streptomycin biosynthesis',    'Glutathione metabolism',    'Two-component system',    'beta-Alanine metabolism',    'Carbon fixation in photosynthetic organisms',    'Biotin metabolism',    'Purine metabolism',    'Arginine and proline metabolism',    'D-Glutamine and D-glutamate metabolism',    'D-Alanine metabolism',    'Pyrimidine metabolism',    'Pantothenate and CoA biosynthesis',    'One carbon pool by folate',    'Lipoic acid metabolism',    'C5-Branched dibasic acid metabolism')
	
	if(prefilter.housekeeping.pathways){
		mbf[['pathway3']] <- mbf[['pathway3']][,!(colnames(mbf[['pathway3']]) %in% ignore.pathways)]
	}
	if(prefilter.outlier.pathways){
		mbf[['pathway3']] <- mbf[['pathway3']][,!(colnames(mbf[['pathway3']]) %in% outlier.pathways)]
	}
	
    params <- list(
        d=d,
        predictor=predictor,
        targets=targets,
        test.type=test.type,
        mbf.transform=mbf.transform,
        mb.transform=mb.transform,
        cohorts=cohorts,
        min.age=min.age,
        max.age=max.age,
        min.years=min.years,
        min.mb.prevalence=min.mb.prevalence,
 	    max.mb.prevalence=max.mb.prevalence,
        min.mb.average=min.mb.average,
        inflamed=inflamed,
        np.test.against=np.test.against,
        nperm=nperm,
        target.names=target.names,
        alpha=alpha,
        drop.outliers=drop.outliers,
        n.bdiv.axes=n.bdiv.axes,
        min.maf=min.maf,
        verbose=verbose,
        print.results=print.results,
        plot.results=plot.results,
        export.results=export.results,
        linear.power.transform=linear.power.transform,
        filebase=filebase, 
        covariates.override=covariates.override,
        prefilter.housekeeping.pathways=prefilter.housekeeping.pathways,
        prefilter.outlier.pathways=prefilter.outlier.pathways,
        use.qvalue=use.qvalue,
        max.pct.zero=max.pct.zero,
        include.gx.pc=include.gx.pc,
        top.n.results=top.n.results,
        ncores=ncores
    )

    attr(results,'params') <- params


    subject.filter <- rep(TRUE,nrow(map))
    if(!is.na(inflamed)) subject.filter <- subject.filter & is.inflamed == inflamed
    
    subject.filter <- subject.filter & map$Age <= max.age
    subject.filter <- subject.filter & map$Age >= min.age
    subject.filter <- subject.filter & map$Years_since_diagnosis >= min.years
    cat(sum(subject.filter),' subjects after filtering.\n')
    if(is.na(cohorts[1])) cohorts <- rep(TRUE, nrow(map))
    if(!is.list(cohorts)) cohorts <- list('cohorts'=cohorts)
    for(i in 1:length(cohorts)) cohorts[[i]] <- cohorts[[i]] & subject.filter
    nsubj <- sapply(cohorts,sum)
    cat('Total', sum(nsubj),' subjects in cohorts after filtering',
    	paste(nsubj,collapse=', '),'\n')

    # collapse covariates, add gx PC'S
    covariates.i <- covariates
    if(!is.null(covariates.override)) covariates.i <- covariates.override
    mm <- collapse.by.mutual.information(covariates.i, subset.ix=subject.filter)
    # include genetic principal components
    if(include.gx.pc){
        mm <- cbind(d$gxx.pc,mm)
        colnames(mm)[1:3] <- sprintf('PC%d',1:3)
    }
    # mb filtering
    if(!is.null(min.mb.prevalence)) mb <- mb[,colMeans(mb > 0) >= min.mb.prevalence]
    if(!is.null(max.mb.prevalence)) mb <- mb[,colMeans(mb > 0) <= max.mb.prevalence]
    if(!is.null(min.mb.average)) mb <- mb[,colMeans(mb) >= min.mb.average]
    
    
    # set up predictor
    x.list <- list()
    if(class(predictor) != 'character'){
        x.list[['Custom']] <- predictor
    } else if(predictor == 'NOD2'){
        x.list[['NOD2']] <- as.matrix(data.frame('NOD2'=rowSums(gx[,daly.snps])))
        x.list[['NOD2']][x.list[['NOD2']] > 2] <- 2
    } else if (predictor == '163+NOD2'){
        x.list[['163+NOD2']] <- cbind(gx[,jostins.snps],as.matrix(data.frame('NOD2'=rowSums(gx[,daly.snps]))))
        maf <- (.5 - abs(apply(x.list[['163+NOD2']],2,function(xx) mean(drop.na(xx))/2)-.5))
        pct.na <- colMeans(is.na(x.list[['163+NOD2']]))
        x.list[['163+NOD2']][,'NOD2'][x.list[['163+NOD2']][,'NOD2'] > 2] <- 2
        x.list[['163+NOD2']] <- x.list[['163+NOD2']][,(maf >= min.maf & pct.na < .1) | colnames(x.list[['163+NOD2']]) == 'NOD2']
    } else if (predictor == 'IChip'){
        x.list[['IChip']] <- gxx
        maf <- (.5 - abs(apply(x.list[['IChip']],2,function(xx) mean(drop.na(xx))/2)-.5))
        x.list[['IChip']] <- x.list[['IChip']][,maf >= min.maf]
    } else if (predictor == 'null'){
        null.uids <- scan('~/drive/prism/data/genetic/snp-subsets/null-ichip-uid.txt',w='c')
        uids <- sapply(strsplit(colnames(gxx.null),'_'),function(xx) paste(xx[1:3],collapse='_'))
        keep.ix <- drop.na(match(null.uids,uids))        
        x.list[['IChip-null']] <- gxx.null[,keep.ix]
        maf <- (.5 - abs(apply(x.list[['IChip-null']],2,function(xx) mean(drop.na(xx))/2)-.5))
        x.list[['IChip-null']] <- x.list[['IChip-null']][,maf >= min.maf]
    } else if (predictor == 'null-500'){
        null.uids <- scan('~/drive/prism/data/genetic/snp-subsets/null-ichip-uid.txt',w='c')
        uids <- sapply(strsplit(colnames(gxx.null),'_'),function(xx) paste(xx[1:3],collapse='_'))
        keep.ix <- drop.na(match(null.uids,uids))        
        x.list[['IChip-null']] <- gxx.null[,sample(keep.ix,500)]
        maf <- (.5 - abs(apply(x.list[['IChip-null']],2,function(xx) mean(drop.na(xx))/2)-.5))
        x.list[['IChip-null']] <- x.list[['IChip-null']][,maf >= min.maf]
    } else if (predictor == 'null-1000'){
        null.uids <- scan('~/drive/prism/data/genetic/snp-subsets/null-ichip-uid.txt',w='c')
        uids <- sapply(strsplit(colnames(gxx.null),'_'),function(xx) paste(xx[1:3],collapse='_'))
        keep.ix <- drop.na(match(null.uids,uids))        
        x.list[['IChip-null']] <- gxx.null[,sample(keep.ix,1000)]
        maf <- (.5 - abs(apply(x.list[['IChip-null']],2,function(xx) mean(drop.na(xx))/2)-.5))
        x.list[['IChip-null']] <- x.list[['IChip-null']][,maf >= min.maf]
    } else if (predictor == 'RX Pathways'){
        disease.scores.combined <- disease.scores$CD[,c('FULL','RX2','RX3','RX4')]
        disease.scores.combined[has.uc] <- disease.scores$UC[has.uc,c('FULL','RX2','RX3','RX4')]
        x.list[['Host.Pathways']] <- disease.scores.combined
    } else if(predictor == 'DK Pathways'){
        disease.scores.combined <- disease.scores$CD[,c('FULL','FUT','AUTOPHAGY','LIGAND_SENSING')]
        disease.scores.combined[has.uc] <- disease.scores$UC[has.uc,c('FULL','FUT','AUTOPHAGY','LIGAND_SENSING')]
        x.list[['Host.Pathways']] <- disease.scores.combined
    } else {
        stop(paste('Unknown predictor',predictor))
    }   

    if(is.null(filebase)) filebase <- paste(predictor,test.type,sep='_')
    
    t.type <- test.type
    f.types <- targets
    f.type.names <- target.names
    if(is.na(f.type.names)) f.type.names <- f.types

    for(x.type in names(x.list)){
        for(f.type.i in 1:length(f.types)){
            f.type <- f.types[f.type.i]
            if(f.type == 'mb') {
                if(mb.transform == 'asin-sqrt'){
                    y <- asin(sqrt(mb))
                } else if(mb.transform == 'pseudolog'){
                    y <- pseudolog(mb)
                } else if(mb.transform == 'power'){
                    y <- power.transform(mb)
                    y <- (y - mean(y)) / sd(y)
                } else if(mb.transform == 'presence-absence'){
                    y <- (mb > 0) + 0
                } else if(mb.transform == 'split-median'){
                    y <- apply(mb,2,function(xx) xx > median(xx))
                } else if(mb.transform == 'none'){
                    y <- mb
                } else {
                    stop(paste('Unknown mb.transform:',mb.transform))
                }
                
                if(t.type == 'chisq' && mb.transform != 'presence-absence' && mb.transform != 'split-median') y <- mb > 0
                
            } else if(f.type == 'alpha'){
                y <- as.matrix(dx[,1,drop=F])
            } else if(f.type == 'beta' || f.type == 'beta-pathway3' || f.type == 'beta-uuf' || f.type == 'beta-wuf' || f.type == 'beta-uuf-wuf'){
            	if(f.type == 'beta'){
	                bdiv.types <- c('pathway3-bc','wuf','uuf')
	            } else if(f.type == 'beta-pathway3') {
	                bdiv.types <- c('pathway3-bc')
	            } else if(f.type == 'beta-uuf'){
	            	bdiv.types <- c('uuf')
   	            } else if(f.type == 'beta-wuf'){
	            	bdiv.types <- c('wuf')
	            } else if(f.type == 'beta-uuf-wuf'){
	            	bdiv.types <- c('uuf','wuf')
	            }

                n.axes <- n.bdiv.axes
                y <- matrix(0,nrow(map),n.axes*length(bdiv.types))
                colnames(y) <- character(n.axes * length(bdiv.types))
                for(ii in seq_along(bdiv.types)){
                    # get principal coords for only this subset of samples
                    if(is.list(cohorts)){
                        for(jj in seq_along(cohorts)){
                            ix.jj <- cohorts[[jj]]
                            y[ix.jj,1:n.axes + n.axes*(ii-1)] <- cmdscale(bdiv[[bdiv.types[ii]]][ix.jj,ix.jj],n.axes)
                        }
                    } else {
                        y[cohorts,1:n.axes + n.axes*(ii-1)] <- cmdscale(bdiv[[bdiv.types[ii]]][cohorts,cohorts],n.axes)
                    }
                    colnames(y)[1:n.axes + n.axes*(ii-1)] <- sprintf('%s-PC%d',bdiv.types[ii],1:n.axes)
                }
            } else {
                if(mbf.transform == 'asin-sqrt'){
                    y <- asin(sqrt(mbf[[f.type]]))
                } else if(mbf.transform == 'pseudolog'){
                    y <- pseudolog(mbf[[f.type]])
                } else if(mbf.transform == 'power'){
                    y <- power.transform(mbf[[f.type]])
                    y <- (y - mean(y)) / sd(y)
                } else if(mbf.transform == 'presence-absence'){
                    y <- mbf[[f.type]] > 0
                } else if(mbf.transform == 'split-median'){
                    y <- apply(mbf[[f.type]],2,function(xx) xx > median(xx))
                } else if(mbf.transform == 'none'){
                    y <- mbf[[f.type]]
                } else {
                    stop(paste('Unknown mbf.transform:',mbf.transform))
                }
                
                if(t.type == 'chisq' && mbf.transform != 'presence-absence' && mbf.transform != 'split-median') y <- mbf > 0
                
            }

            if(print.results) cat(f.type.names[f.type.i], '\n')			

            # run.test            
            res <- stratified.association.tests(
                x.list[[x.type]],
                y,
                print.results=FALSE,
                warn=TRUE,
                test.type=t.type,
                snp.subsets=NULL, qtr.subsets=NULL,
                covariates=mm, 
                subject.subsets=cohorts,
                test.against=np.test.against,
                drop.outliers=drop.outliers,
                do.power.transform=linear.power.transform,
                nperm=nperm,
                verbose=verbose,
                top.n.results=top.n.results,
                use.qvalue=use.qvalue,
                max.pct.zero=max.pct.zero,
                ncores=ncores
            )
                
            # print results if any
            sig <- res[,'qvalue'] < alpha
            if(f.type == 'alpha' | f.type == 'beta'){
                sig <- res[,'pvalue'] < .05
            }

            if(sum(sig,na.rm=TRUE) > 0){
                if(print.results) print(res[which(sig)[1:min(10,sum(sig,na.rm=T))],])

                if(export.results){
                    if(!is.null(filebase)) sink(sprintf('%s-%s.xls',filebase,f.type.names[f.type.i]))
                    write.table(res[which(sig),],col=T,row=F,quo=F,sep='\t')
                    if(!is.null(filebase)) sink(NULL)
                }

                # plot
                if(plot.results){
                    x <- x.list[[x.type]]
                    if(t.type == 'np' | t.type == 'np-nonzero') x <- x > 0
                    if(is.null(dim(x))) x <- matrix(x,ncol=1)
                    if(!is.null(filebase)) pdf(sprintf('%s-%s.pdf',filebase,f.type.names[f.type.i]),width=8,height=8)
                    n.per.plot <- 4
                    k <- 0
                    while(k * n.per.plot < sum(sig,na.rm=TRUE) & sum(sig,na.rm=TRUE) > 0){
                        new.ix <- k * n.per.plot + 1:n.per.plot
                        new.ix <- new.ix[which(new.ix <= sum(sig,na.rm=TRUE))]
                        ix <- cohorts[[1]]
                        plot.association.results(
                            x[ix,,drop=F],
                            y[ix,,drop=F], covariates=mm[ix,,drop=F],
                            assoc=res[new.ix,,drop=F],add.to.plot=FALSE,
#                             drop.outliers.range=ifelse(is.null(drop.outliers),2,drop.outliers),
                            drop.outliers.range=drop.outliers,
                            regression.type == ifelse(t.type=='linear' | t.type == 'linear-nonzero','linear','negbin'),
                            plot.type=c('vioplot','boxplot')[1]
                                , nonzero.only= t.type == 'linear-nonzero' | t.type == 'np-nonzero'
#                               , nonzero.only= t.type == 'linear-nonzero' | (f.type == 'mb' && t.type == 'linear') | t.type == 'np-nonzero'
 #                                      , color.by=factor(apply(gx[ix,daly.snps],1,function(xx) ifelse(xx[6]>0,6,ifelse(any(xx > 0),which.max(xx),0))))
                                   , color.by=factor(apply(gx[ix,daly.snps], 1, function(xx) ifelse(any(xx > 0),which.max(xx),0)))
                                 , color.by.color=c('#99999944',sprintf('%s99',c("#CFB8DB","#81ADD3","#649D2C","#7D661C","#7E2712","#380024",brewer.pal(9,'Set1')[-6]))),
                                  , legend.names=c('Unaffected',colnames(gx)[daly.snps])
#                                      , color.by=factor(gx[ix,'rs5743293'] > 0)
#                                      , color.by.color=c('#99999944',sprintf('%s99',c("#380024",brewer.pal(9,'Set1')[-6]))),
#                                      , legend.names=c('Unaffected','rs5743293')
#                                 , color.by=factor(rep('0',sum(ix)))
#                                 , color.by.color='#99999944'
#                                 , legend.names='Unaffected'
                            )
                        k <- k + 1
                    }
                    if(!is.null(filebase)) dev.off()
                }
            }
            res$Gene <- topgenes[res[,1]]
            results[[f.type]] <- res
        }
    }
    return(results)
})

"rarefy" <- function(x, rare.depth=5000, keep.low.samples=FALSE){
  # ensure x integer valued
  if(any((x - floor(x)) > 0))stop('input data x must be integer-valued')
  if(is.null(rare.depth)) return(x)
  
  if(!is.element(class(x), c('matrix', 'data.frame','array')))
    x <- matrix(x,nrow=1)
  nr <- nrow(x)
  nc <- ncol(x)
  
  for(i in 1:nrow(x)){
    if(sum(x[i,]) > rare.depth){
      # get taxon index for each sequence
      D <- sum(x[i,])
      tax.cumsum <- cumsum(x[i,])
      tax.ix <- sapply(1:D,function(x) min(which(x[i,]<=tax.cumsum)))
      s <- tax.ix[sample(D,size=max.depth)]
      x[i,] <- hist(s,breaks=seq(.5,nc+.5,1), plot=FALSE)$counts
    }
  }
  return(x)
}
