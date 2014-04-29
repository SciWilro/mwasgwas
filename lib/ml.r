# Tests whether QTs can be predicted from genetics
# returns yhat, rsq/acc, pval, qval
# pvals are permutation-based (permuting genetics)
# allows prediction of residuals after regressing on covariates (for regression)
# note: qtrs must be a numeric to perform regression
"predict.complex.trait" <- function(gx,qtr,covariates=NULL,nperm=999,warn=TRUE,
        regression.type=c('linear','lars')[1]){
    require('randomForest')
    if(is.null(dim(qtr))) qtr <- matrix(qtr,ncol=1)
    is.regression <- is.numeric(qtr[,1])
    
    if(!is.null(covariates)){
        na.ix <- colMeans(is.na(covariates)) > .5
        # drop any metadata columns with > .5 NA
        if(any(na.ix)){
            if(warn) cat(sprintf('Warning: dropping metadata column %s due to >50 %% NA.\n',colnames(covariates)[na.ix]))
            covariates <- covariates[,!na.ix,drop=F]
        }
        na.ix <- rowSums(is.na(covariates)) > 0
        if(any(na.ix)){
         if(warn) cat(sprintf('Warning: dropping %d samples due to NA values in covariates.\n',sum(na.ix)))
            gx <- gx[!na.ix,,drop=F]
            qtr <- qtr[!na.ix,,drop=F]
        }
        covariates <- droplevels(covariates[!na.ix,,drop=F])
        if(nrow(gx) != nrow(covariates)){
            stop(sprintf('Error: gx has %d rows, mb has %d rows, covariates has %d rows.\n',
                nrow(gx), nrow(qtr), nrow(covariates)))
        }
        mm <- covariate.model.matrix(covariates)
    }


    yhat <- as.data.frame(qtr)
    perfs <- rep(NA,ncol(qtr)) # performances (rsq or accuracy)
    imps <- matrix(NA,ncol(gx),ncol(qtr)) # per-QT gx importances
    
    for(i in 1:ncol(qtr)) {
        y <- qtr[,i]
        if(is.regression) {
            # use linear regression
            if(!is.null(covariates)){
                mlm <-lm(y ~ mm)
                y <- mlm$residuals;
            }
            if(regression.type == 'linear'){
                mlm <- lm(y ~ gx)
                imps[,i] <- abs(coef(mlm))[-1]
                perfs[i] <- 1 - sum(mlm$residuals**2)/sum((y-mean(y))**2)
                yhat[,i] <- mlm$fitted.values
            } else {
                
            }
        } else {
            y <- factor(y)        
            mrf <- randomForest(gx,y, importance=T)
            yhat[,i] <- mrf$predicted
            imps[,i] <- mrf$importance[,'MeanDecreaseAccuracy']
            baseline <- 1-(max(table(y))/length(y))
            perfs[i] <- baseline - mrf$err.rate[mrf$ntree,1]
            if(!is.finite(perfs[i])) print(baseline)
        }
        if(!is.null(nperm)) cat(i,perfs[i],'')
    }
 
    if(!is.null(nperm)){
        cat('\nPERM:\n')
        perm.perfs <- matrix(NA,ncol(qtr),nperm+1)
        for(i in 1:nperm){
            cat(i,'')
            perm.perfs[,i] <- predict.complex.trait(gx[sample(nrow(gx)),],qtr,covariates,NULL)$perf
        }
        cat('\n')
        perm.perfs[,nperm+1] <- perfs
        
        pvals <- rowMeans(sweep(perm.perfs,1,perfs,'>='))
        qvals <- p.adjust(pvals,method='fdr')
    } else {
        pvals <- NA
        qvals <- NA
    }
    return(list(yhat=yhat,perf=perfs,importance=imps,pvals=pvals,qvals=qvals))
}


# Hold-out predictions using cv-tuned lasso
# nfolds=-1 means leave-one-out; nfolds=10 means 10-fold
# nfolds = NA means don't do hold-out, just tune and predict on all data
"loo.lasso" <- function(x,y,nfolds=10){
    require('lars')
    yhat <- y
    n <- length(y)

    if(is.na(nfolds)){
        # cross-validation tuning of lasso fraction
        cvl <- cv.lars(x,y,normalize=FALSE,plot=F)
        # build single lasso model
        l <- lars(x,y,normalize=FALSE)
        # predict y_i
        yhat <- predict(l,x,s=cvl$index[which.min(cvl$cv)],mode='fraction')$fit
    } else {
        if(nfolds == -1) nfolds <- n
        
        folds <- sample(rep(1:nfolds,length=n))
        for(i in 1:nfolds){
            cat(i,'')
            foldix <- which(folds==i)
            # cross-validation tuning of lasso fraction
            cvl <- cv.lars(x[-foldix,,drop=F],y[-foldix],normalize=FALSE,plot=F)
            # build single lasso model
            l <- lars(x[-foldix,,drop=F],y[-foldix],normalize=FALSE)
            # predict y_i
            yhat[foldix] <- predict(l,x[foldix,,drop=F],s=cvl$index[which.min(cvl$cv)],mode='fraction')$fit
        }
        cat('\n')
    }
    return(yhat)
}

