source('~/drive/research/prism/src/lib/rpackage/subsets.r')


# compare.obs <- run.subsets(x,gx,subsetx=subset9 & m$Collection_Center != 'RISK',alpha=.1,min.prevalence=.5)
# compare.null <- run.subsets(x,gx,subsetx=subset9 & m$Collection_Center != 'RISK',alpha=.1,min.prevalence=.5, permute=TRUE)

# mean.cors.obs <-  sapply(compare.obs$cors.cohort.cor,mean)
# mean.cors.null <-  sapply(compare.null$cors.cohort.cor,mean)


alpha <- .05
hit.ix.2way <- rep(FALSE,nrow(compare.obs$res.all))
hit.ix.3way <- rep(FALSE,nrow(compare.obs$res.all))
for(i in which(rowSums(compare.obs$res.all[,c('pvalue.mgh','pvalue.msh','pvalue.neth')] < alpha) >= 2)){
# 	if(compare.obs$res.all$Covariate[i] == 'NOD2') browser()
	which.hits <- which(compare.obs$res.all[i,c('pvalue.mgh','pvalue.msh','pvalue.neth')] < alpha)
	nhits <- length(which.hits)
	signs <- sign(as.numeric(compare.obs$res.all[i,c('coefficient.mgh','coefficient.msh','coefficient.neth')]))
	if(nhits >= 2){
# 		if(i == 1) browser()
		hit.ix.2way[i] <- max(table(signs[which.hits])) == 2
		if(nhits == 3) hit.ix.3way[i] <- max(table(signs[which.hits])) == 3
	}
}