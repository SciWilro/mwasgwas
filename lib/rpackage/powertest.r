# n samples
# true.cor true correlation (can be a vector)
# conf pct confidence or fraction of true correlations recovered
# FWER desired family-wise error rate
"powertest" <- function(n, true.cor=seq(0,1,.01), fwer=.05,do.plot=FALSE, return.log=TRUE){
    require(mvtnorm)
    
    # generates a p-value for correlation test of two normals
    # with given true correlation
    pFun <- function(n, true.cor){
                        x <- rmvnorm(n, sigma=matrix(c(1, true.cor, true.cor, 1), ncol=2))
                        cor.test(x[,1], x[,2], method="spearman")$p.value
                    }
    
    cor.range <- true.cor
    
    # For a given true correlation,
    # find the 90% quantile of p-values
    # i.e. 90% of p-values are better than this value
    # i.e. we have a 90% chance of being significant if this is the cutoff
    p.val <- sapply(cor.range, function(rho) quantile(replicate(100, pFun(150, rho)), probs=0.90) )
    names(p.val) <- as.character(cor.range)
    
    if(do.plot){
        # for a given number of tests, n, with true correlation, say 0.5,
        # we will be testing at alpha = 0.05 / n
        # so n = 0.05 / alpha
        # alpha is chosen so that we'll recover 90% of the associations at that level
        pdf("pval_vs_cor.pdf", useDingbats=FALSE)
        plot(-log(p.val, 10) ~ cor.range,
             xlab="Spearman correlation",
             main=sprintf("-log10(nominal p-value) vs true Spearman correlation for %d samples",n) )
        abline(v=seq(0.1, 0.9, by=0.1), lty=3)
        abline(h=1:12, lty=3)
        plot(log(0.05/p.val, 10) ~ cor.range,
             xlab="Spearman correlation",
             ylab="log10(# pairwise comparisons possible)",
             main=sprintf("-log10(# comparisons possible for FWER < 0.05) vs true Spearman correlation\n for %d samples",n) )
        abline(v=seq(0.1, 0.9, by=0.1), lty=3)
        abline(h=1:12, lty=3)
        dev.off()
    }

    if(return.log){
        return(log10(0.05) - log10(p.val))
    } else {
        return(0.05/p.val)
    }

}



"powertest.t" <- function(n, delta=seq(0.5,2,.01), fwer=.05){
	pvals <- sapply(delta,function(xx) power.t.test(n,delta=xx,sd=1,sig.level=NULL,power=.9)$sig.level)
	ncomps <- 0.05 / pvals
	return(ncomps)
}


# bootstrap estimate of wilcoxon power
# x: independent variable
# y: grouping
"powertest.mc" <- function(x,y,n=NULL,alpha=0.05,pow=.9,nperm=10000, filepath='powertest.mc.pdf'){
    if(is.null(n)) n <- max(table(y))
    # generate fake microbiomes
    classx <- split(x,y)
    pvals <- numeric(nperm)
    for(i in 1:nperm){
        simx <- sapply(classx,function(xx) sample(xx,size=n,replace=TRUE))
        pvals[i] <- wilcox.test(simx[,1],simx[,2],exact=FALSE)$p.val
    }
    
    ncomps <- seq(100,10000,1000)
    pows <- sapply(ncomps, function(xx) mean(pvals < alpha/xx))

    fn <- function(ncomp) abs(pow - mean(pvals < alpha/ncomp))
    res <- optimize(fn,lower=1,upper=1e7)$minimum
    
    # plot
    pdf(filepath,width=4,height=4)
#     browser()
    plot(sort(c(ncomps,res)),c(pows,pow)[order(c(ncomps,res))],xlab='Number of comparisons',ylab='Statistical power',
            cex.axis=.8, cex.lab=.9, type='b')
    lines(c(0,res,res),c(pow,pow,0),lty=2)
    dev.off()
    
    
    return(floor(res))
    
}