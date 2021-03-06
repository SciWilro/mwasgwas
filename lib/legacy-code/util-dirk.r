# rv.id is the ID of the random variable (subject ID column header is the intention)
"association.tests.mixed" <- function(X,Y,random.id=NULL,mixed=TRUE,
									 remove.vowels.from.taxon.names=TRUE,
									 include.residuals=FALSE,
									 drop.outlier.range=3,
									 use.qvalue=TRUE
                            ){
	rY <- Y # residuals
	rY[,] <- NA
	.X.72828 <<- as.data.frame(X)
	na.ix <- rowSums(is.na(.X.72828)) > 0
	.X.72828 <- droplevels(.X.72828[!na.ix,])
	Y <- Y[!na.ix,]
    
	fixed.ids <- setdiff(colnames(X),random.id)
	# assumes dependent var will be "y"
	fixed.formula <- as.formula(paste("y.72828", paste(fixed.ids, collapse=" + "), sep=" ~ "))
	if(mixed){
		random.formula <- as.formula(sprintf('~ 1 | %s', random.id))
	}

	index <- 1

	# make one model matrix to determine number of covariates
	y.72828 <<- Y[,1]
	n.covar <- ncol(model.matrix(fixed.formula, data=.X.72828)) - 1
	res <- matrix(0,nrow=n.covar * ncol(Y), ncol=3); # will be numeric
	desc <- matrix("", nrow=nrow(res), ncol=2) # will be character
	colnames(res) <- c('pvalue','qvalue','coefficient')
	colnames(desc) <- c('Covariate','Taxon')
	rx <- Y # residuals
	rx[,] <- NA
    count <- 1
	for(i in 1:ncol(Y)){
		y.72828 <<- Y[,i]
		outlier.72828 <<- is.outlier(y.72828,range=drop.outlier.range)
		if(!(var(y.72828[!outlier.72828]) == 0 | sum(outlier.72828) > .1 * length(y.72828))){

			if(mixed){
				m <- try(lme(fixed=fixed.formula,
						random=random.formula,
						data=.X.72828,
						subset=!outlier.72828))
				if(class(m)=='try-error'){
					res[ix,c('pvalue','coefficient')] <- c(NA,NA)
					rY[rownames(Y),i][!outlier.72828] <- NA
					cat('Error fitting',colnames(Y)[i],'\n')
				} else {
					res.i <- summary(m)[[20]]
					ix <- count:(count + nrow(res.i) - 2)
					res[ix,c('pvalue','coefficient')] <- res.i[-1,c('p-value','Value')]
					rY[rownames(Y),i][!outlier.72828] <- resid(m)
				}
			} else {
				m <- lm(fixed.formula, data=.X.72828, subset=!outlier.72828)
				res.i <- summary(m)[[4]]
				ix <- count:(count + nrow(res.i) - 2)
				res[ix,c('pvalue','coefficient')] <- res.i[-1,c('Pr(>|t|)','Estimate')]
				rY[rownames(Y),i][!outlier.72828] <- resid(m)
			}
		}
		desc[ix,'Taxon'] <- colnames(Y)[i]
		desc[ix,'Covariate'] <- rownames(res.i)[-1]
    	count <- max(ix) + 1
	}
	# if res is not full, drop empty rows
	if(count <= nrow(res)){
	    res <- res[-(count:nrow(res)),]
	    desc <- desc[-(count:nrow(desc)),]
	}
	desc <- desc[!is.na(res[,'pvalue']),]
	res <- res[!is.na(res[,'pvalue']),]
    if(use.qvalue){
        require('qvalue')
        qvals <- try(qvalue(res[,'pvalue'])$qvalue,silent=TRUE)
    }
    if(!use.qvalue || class(qvals)=='try-error'){
		qvals <- p.adjust(res[,'pvalue'],'fdr')
    }
    res[,'qvalue'] <- qvals

	desc <- desc[order(res[,'qvalue']),]
	res <- res[order(res[,'qvalue']),]
	res.df <- as.data.frame(desc)
	res.df <- cbind(res.df,res)
	
	# make taxon names shorter
	taxon.names <- res.df$Taxon
	for(to.replace in c('k__',' p__',' c__',' o__',' f__',' g__',' s__')){
		taxon.names <- gsub(to.replace,'',taxon.names)
	}
	if(remove.vowels.from.taxon.names){
		taxon.names <- gsub('[aeiou]','',taxon.names)
	}
	res.df$Taxon_abbrev <- taxon.names

	if(include.residuals){
		return(list(res=res.df, residuals=rY))
	} else {
		return(res.df)
	}
}

"association.tests.mixed.with.partial.residuals" <-function(X,Y,random.id=NULL,mixed=TRUE,
									 remove.vowels.from.taxon.names=TRUE, drop.outlier.range=3,
									 use.qvalue=TRUE){
	full.res <- association.tests.mixed(X,Y,random.id,mixed,
					remove.vowels.from.taxon.names,include.residuals=TRUE, drop.outlier.range=drop.outlier.range, use.qvalue=use.qvalue)
	partial.residuals <- list()
	random.var.ix <- which(colnames(X) == random.id)
	for(i in (1:ncol(X))[-random.var.ix]){
		cat('Getting partial residuals for covariate ',colnames(X)[i],'...\n',sep='')
		partial.residuals[[colnames(X)[i]]] <- 
			association.tests.mixed(X[,-i],Y,random.id,mixed,
				remove.vowels.from.taxon.names,include.residuals=TRUE)$residuals
	}
	return(list(res=full.res$res,
			residuals=full.res$residuals,
			partial.residuals=partial.residuals
	))
}



# runs cross-validation 
# if predict.fun is NULL, uses S3 predict method
# if nfolds > length(y) or nfolds==-1, uses leave-one-out cross-validation
# ...: additional parameters for train.fun
#
# value:
# y: true values
# predicted: cv predicted values
# probabilities: cv predicted class probabilities (or NULL if unavailable)
# confusion.matrix: confusion matrix (true x predicted)
# nfolds: nfolds
# params: list of additional parameters
# importances: importances of features as predictors
"rf.cross.validation.regression" <- function(x, y, nfolds=10, verbose=verbose, ...){
    if(nfolds==-1) nfolds <- length(y)
    folds <- sample(rep(1:nfolds,ceiling(length(y)/nfolds)))
    
    result <- list()
    result$y <- y
    result$predicted <- result$y
    result$importances <- matrix(0,nrow=ncol(x),ncol=nfolds)
    result$errs <- numeric(length(unique(folds)))

    # K-fold cross-validation
    for(fold in sort(unique(folds))){
        if(verbose) cat(sprintf('Fold %d...\n',fold))
        foldix <- which(folds==fold)
        model <- randomForest(x[-foldix,], factor(result$y[-foldix]), importance=TRUE, do.trace=verbose, ...)
        newx <- x[foldix,]
        if(length(foldix)==1) newx <- matrix(newx,nrow=1)
        result$predicted[foldix] <- predict(model, newx)
        probs <- predict(model, newx, type='prob')
        result$probabilities[foldix,colnames(probs)] <- probs
        result$errs[fold] <- mean(result$predicted[foldix] != result$y[foldix])
        result$importances[,fold] <- model$importance[,'MeanDecreaseAccuracy']
    }

	result$nfolds <- nfolds
    result$params <- list(...)
    result$confusion.matrix <- t(sapply(levels(y), function(level) table(result$predicted[y==level])))
    return(result)    
}


# d is the data list
# d should contain at least:
# map: with 'PCDAI_06' column real values, to be predicted
# cvr: covariates table
# taxa: sample x taxon relative abundance matrix
# 
# if x is null, uses taxa object
"mc.test.taxa.prediction" <- function(d, y, x=NULL, nreps.obs=10, nreps.mc=10, nfolds=10, ntree=2000,
		covariate.names=NULL,
		subset.ix=NULL){ with(d, {
	# y <- rowMeans(map[,c('PCDAI_12','PCDAI_06')]) - map$PCDAI
	# y <- rowMeans(map[,c('PCDAI_12','PCDAI_06')])

	if(is.null(covariate.names)){
		cvr.i <- cvr
	} else {
		cvr.i <- map[,covariate.names]
	}
	
	if(is.null(subset.ix)){
		subset.ix <- 1:length(y)
	}
	
	if(!is.null(x)) taxa <- x
	obs.errs <- NULL
	probabilities <- matrix(NA,nrow=length(subset.ix),ncol=nreps.obs)
	rownames(probabilities) <- rownames(taxa)[subset.ix]
	
	for(i in 1:nreps.obs){
		cat('Obs iteration',i,'\n')
		y.i <- y
		x <- cbind(taxa, cvr.i[,colnames(cvr.i) != 'AnSubID'])
		y.i <- y.i[subset.ix]
		x <- x[subset.ix,]
		
		na.ix <- is.na(y.i) | rowSums(is.na(x)) > 0
		y.i <- y.i[!na.ix]
		x <- x[!na.ix,]

		ix <- x$Race != "Other"
		y.i <- factor(y.i[ix])
		x <- droplevels(x[ix,])

		rf.cv <- rf.cross.validation(x,y.i,
				nfolds=nfolds, verbose=FALSE,ntree=ntree,mtry=ncol(x)**.75)
		cat('OBS error:',mean(rf.cv$errs),'\n')
		probabilities[rownames(rf.cv$probabilities),i] <- rf.cv$probabilities[,2]
		if(i == 1) {
			y.to.return <- y.i
			names(y.to.return) <- rownames(x)
			predicted <- rf.cv$predicted
			names(predicted) <- rownames(x)
			importances <- rf.cv$importances
			rownames(importances) <- colnames(x)
			print(table(y.i,predicted))
			
		}
		
		obs.errs <- c(obs.errs, mean(rf.cv$errs))
	}
	cat('Error:',mean(obs.errs),'\n')
	probabilities <- rowMeans(probabilities)
	probabilities <- probabilities[!is.na(probabilities)]
	
	baseline.errs <- NULL
	mc.probabilities <- matrix(NA,nrow=length(subset.ix),ncol=nreps.mc)
	rownames(mc.probabilities) <- rownames(taxa)[subset.ix]

	if(nreps.mc > 0){
		for(i in 1:nreps.mc){
			cat('MC iteration',i,'\n')
			y.i <- y
			x <- cbind(taxa[sample(1:nrow(taxa)),], cvr.i[,colnames(cvr.i) != 'AnSubID'])
			y.i <- y.i[subset.ix]
			x <- x[subset.ix,]
		
			na.ix <- is.na(y.i) | rowSums(is.na(x)) > 0
			y.i <- y.i[!na.ix]
			x <- x[!na.ix,]

			ix <- x$Race != "Other"
	 		y.i <- factor(y.i[ix])
			x <- droplevels(x[ix,])
	
			rf.cv.mc <- rf.cross.validation(x,y.i, nfolds=nfolds, verbose=FALSE,ntree=ntree,mtry=ncol(x) ** .75)
			cat('MC error:',mean(rf.cv.mc$errs),'\n')
			baseline.errs <- c(baseline.errs, mean(rf.cv.mc$errs))
			mc.probabilities[rownames(cvr.i)[subset.ix][!na.ix][ix],i] <- rf.cv.mc$probabilities[,2]
		}
	}
	mc.probabilities <- rowMeans(mc.probabilities)
	mc.probabilities <- mc.probabilities[!is.na(mc.probabilities)]
	
	cat('\n')
	cat(sprintf('Observed error: %.3f +/- %.3f (%d reps)\n',
				mean(obs.errs),
				sd(obs.errs)/sqrt(nreps.obs),
				nreps.obs))
			
	cat(sprintf('MC error: %.3f +/- %.3f (%d reps)\n',
				mean(baseline.errs),
				sd(baseline.errs)/sqrt(nreps.mc),
				nreps.mc))
			
	pvalue <- (sum(baseline.errs <= mean(obs.errs)) / (1 + nreps.mc))
	cat(sprintf('MC p-value = %.5f\n', pvalue))

	return(list(obs.errs=obs.errs, mc.errs=baseline.errs,  pvalue=pvalue, predicted=predicted, y=y.to.return,
				importances=importances, probabilities=probabilities,
				mc.probabilities=mc.probabilities))
	
})}
