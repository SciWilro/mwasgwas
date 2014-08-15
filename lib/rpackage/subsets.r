# compares results from MGH, MSH, NETH
"run.subsets" <- function(x, gx, subsetx=NULL, permute=FALSE,
						   cor.type=c('value','sign')[2],
						   cor.method=c('pear','spear')[1],
						   min.prevalence=.5,
						   max.missing.gx=.05,
						   alpha=.05, do.plot=FALSE){

	MIN.PREV <- min.prevalence
	OUTLIER.RANGE <- 3
	ALPHA <- alpha
	pc.k <- 0

	snps <- read.table('~/drive/research/prism/data/genetic/2013-01-09-generate-snps/all-snps-metadata-with-pos-with-categories-with-NOD2.txt',sep='\t',head=T,row=1,check=F)

	ix <- subsetx
	pvals.crosscohort <- list()
	cors.crosscohort <- list()
	xyz.all <- NULL

	snps.to.try <- colnames(gx)[colSums(is.na(gx[ix,,drop=F])) < ceiling(max.missing.gx * sum(ix))]
	for(snp in snps.to.try){
		# for(snp in 'NOD2'){
		cat('Comparing cohorts for SNP',snp,'\n')

		gxx <- gx[,snp,drop=F]
		if(permute){
			gxx <- gxx[sample(1:nrow(gxx),size=nrow(gxx)),,drop=F]
		}
		xx <- x

		# if(pc.k > 0){
		# 	pcoa.gx1 <- cmdscale(dgx[ix,ix],k=pc.k); colnames(pcoa.gx1) <- sprintf('PC%d',1:pc.k)
		# 	Z <- cbind(m[ix,clin],pcoa.gx1)
		# } else {
		# 	Z <- m[ix,clin]
		# }

		ix <- subsetx & m$is.mgh
		if(pc.k > 0){
			pcoa.gx1 <- cmdscale(dgx[ix,ix],k=pc.k); colnames(pcoa.gx1) <- sprintf('PC%d',1:pc.k)
			Z <- cbind(m[ix,clin],pcoa.gx1)
		} else {
			Z <- m[ix,clin]
		}



		cat('MGH:',sum(ix),'\n')
		gx.stats <- mwas.xwas.association.tests(gxx[ix,,drop=F], asin(sqrt(xx[ix,colMeans(xx > 0) > MIN.PREV])), Z,drop.outlier.range=OUTLIER.RANGE,verbose=FALSE)
		gx.stats$gene <- snps$topgene[match(gx.stats$Covariate,rownames(snps))]
		gx.stats.mgh <- gx.stats;

		# print(gx.stats.mgh[1:10,])

		ix <- subsetx & m$is.msh
		if(pc.k > 0){
			pcoa.gx1 <- cmdscale(dgx[ix,ix],k=pc.k); colnames(pcoa.gx1) <- sprintf('PC%d',1:pc.k)
			Z <- cbind(m[ix,clin],pcoa.gx1)
		} else {
			Z <- m[ix,clin]
		}
		cat('MSH:',sum(ix),'\n')
		gx.stats <- mwas.xwas.association.tests(gxx[ix,,drop=F], asin(sqrt(xx[ix,colMeans(xx > 0) > MIN.PREV])), Z,drop.outlier.range=OUTLIER.RANGE,verbose=FALSE)
		gx.stats$gene <- snps$topgene[match(gx.stats$Covariate,rownames(snps))]
		gx.stats.msh <- gx.stats;

		ix <- subsetx & m$is.neth
		if(pc.k > 0){
			pcoa.gx1 <- cmdscale(dgx[ix,ix],k=pc.k); colnames(pcoa.gx1) <- sprintf('PC%d',1:pc.k)
			Z <- cbind(m[ix,clin],pcoa.gx1)
		} else {
			Z <- m[ix,clin]
		}
		cat('NETH:',sum(ix),'\n')
		gx.stats <- mwas.xwas.association.tests(gxx[ix,,drop=F], asin(sqrt(xx[ix,colMeans(xx > 0) > MIN.PREV])), Z,drop.outlier.range=OUTLIER.RANGE,verbose=FALSE)
		gx.stats$gene <- snps$topgene[match(gx.stats$Covariate,rownames(snps))]
		gx.stats.neth <- gx.stats;

		
		pdf(sprintf('compare-%s-across-cohorts.pdf',snp),width=9,height=3)
		par(mfrow=c(1,3))

		# correlations and plotting
		xy <- merge(gx.stats.mgh, gx.stats.msh, by=c('Covariate','Taxon'))
		sigix <- xy$pvalue.x < ALPHA | xy$pvalue.y < ALPHA
		if(sum(sigix) > 2){
			plot(xy$coefficient.x[sigix], xy$coefficient.y[sigix],main='MGH v MSH')
			if(length(unique(xy$coefficient.x[sigix])) == 1){
				abline(v = xy$coefficient.x[sigix][1])
			} else if(length(unique(xy$coefficient.y[sigix])) == 1){
				abline(h = xy$coefficient.y[sigix][1])
			} else {
				abline(lm(xy$coefficient.y[sigix] ~ xy$coefficient.x[sigix]))
			}
			if(cor.type == 'sign'){
				p.mgh.msh <- cor.test(sign(xy$coefficient.x[sigix]), sign(xy$coefficient.y[sigix]),method=cor.method,exact=F)$p.value
				cor.mgh.msh <- cor.test(sign(xy$coefficient.x[sigix]), sign(xy$coefficient.y[sigix]),method=cor.method,exact=F)$estimate
			}
# 			if(cor.type == 'fisher' || (cor.type == 'sign' && is.na(p.mgh.msh))){
			if(cor.type == 'fisher'){
				sign.table <- table(factor(sign(xy$coefficient.x[sigix]),levels=c(-1,1)), factor(sign(xy$coefficient.y[sigix]),levels=c(-1,1)))
				if(any(rowSums(sign.table) == 0) || any(colSums(sign.table) == 0)) sign.table <- sign.table + 1
				p.mgh.msh <- fisher.test(sign.table,alternative='greater')$p.value
				cor.mgh.msh <- fisher.test(sign.table,alternative='greater')$estimate
			} else {
				p.mgh.msh <- cor.test(xy$coefficient.x[sigix], xy$coefficient.y[sigix],method=cor.method,exact=F)$p.value
				cor.mgh.msh <- cor.test(xy$coefficient.x[sigix], xy$coefficient.y[sigix],method=cor.method,exact=F)$estimate
			}
		} else {
			p.mgh.msh <- NA
			cor.mgh.msh <- NA
		}


		xy <- merge(gx.stats.mgh, gx.stats.neth, by=c('Covariate','Taxon'))
		sigix <- xy$pvalue.x < ALPHA | xy$pvalue.y < ALPHA
		if(sum(sigix) > 2){
			plot(xy$coefficient.x[sigix], xy$coefficient.y[sigix],main='MGH v NETH')
			if(length(unique(xy$coefficient.x[sigix])) == 1){
				abline(v = xy$coefficient.x[sigix][1])
			} else if(length(unique(xy$coefficient.y[sigix])) == 1){
				abline(h = xy$coefficient.y[sigix][1])
			} else {
				abline(lm(xy$coefficient.y[sigix] ~ xy$coefficient.x[sigix]))
			}
			# print(cor.test(xy$coefficient.x[sigix], xy$coefficient.y[sigix],method=cor.method,exact=F))
			if(cor.type == 'sign'){
				p.mgh.neth <- cor.test(sign(xy$coefficient.x[sigix]), sign(xy$coefficient.y[sigix]),method=cor.method,exact=F)$p.value
				cor.mgh.neth <- cor.test(sign(xy$coefficient.x[sigix]), sign(xy$coefficient.y[sigix]),method=cor.method,exact=F)$estimate
			}
# 			if(cor.type == 'fisher' || (cor.type == 'sign' && is.na(p.mgh.neth))){
			if(cor.type == 'fisher'){
				sign.table <- table(factor(sign(xy$coefficient.x[sigix]),levels=c(-1,1)), factor(sign(xy$coefficient.y[sigix]),levels=c(-1,1)))
				if(any(rowSums(sign.table) == 0) || any(colSums(sign.table) == 0)) sign.table <- sign.table + 1
				p.mgh.neth <- fisher.test(sign.table,alternative='greater')$p.value
				cor.mgh.neth <- fisher.test(sign.table,alternative='greater')$estimate
			} else {
				p.mgh.neth <- cor.test(xy$coefficient.x[sigix], xy$coefficient.y[sigix],method=cor.method,exact=F)$p.value
				cor.mgh.neth <- cor.test(xy$coefficient.x[sigix], xy$coefficient.y[sigix],method=cor.method,exact=F)$estimate
			}
		} else {
			p.mgh.neth <- NA
			cor.mgh.neth <- NA
		}
		
		xy <- merge(gx.stats.msh, gx.stats.neth, by=c('Covariate','Taxon'))
		sigix <- xy$pvalue.x < ALPHA | xy$pvalue.y < ALPHA
		if(sum(sigix) > 2){
			plot(xy$coefficient.x[sigix], xy$coefficient.y[sigix],main='MSH v NETH')
			if(length(unique(xy$coefficient.x[sigix])) == 1){
				abline(v = xy$coefficient.x[sigix][1])
			} else if(length(unique(xy$coefficient.y[sigix])) == 1){
				abline(h = xy$coefficient.y[sigix][1])
			} else {
				abline(lm(xy$coefficient.y[sigix] ~ xy$coefficient.x[sigix]))
			}
			# print(cor.test(xy$coefficient.x[sigix], xy$coefficient.y[sigix],method=cor.method,exact=F))
			if(cor.type == 'sign'){
				p.msh.neth <- cor.test(sign(xy$coefficient.x[sigix]), sign(xy$coefficient.y[sigix]),method=cor.method,exact=F)$p.value
				cor.msh.neth <- cor.test(sign(xy$coefficient.x[sigix]), sign(xy$coefficient.y[sigix]),method=cor.method,exact=F)$estimate
			}
# 			if(cor.type == 'fisher' || (cor.type == 'sign' && is.na(p.msh.neth))){
			if(cor.type == 'fisher'){
				sign.table <- table(factor(sign(xy$coefficient.x[sigix]),levels=c(-1,1)), factor(sign(xy$coefficient.y[sigix]),levels=c(-1,1)))
				if(any(rowSums(sign.table) == 0) || any(colSums(sign.table) == 0)) sign.table <- sign.table + 1
				p.msh.neth <- fisher.test(sign.table,alternative='greater')$p.value
				cor.msh.neth <- fisher.test(sign.table,alternative='greater')$estimate
			} else {
				p.msh.neth <- cor.test(xy$coefficient.x[sigix], xy$coefficient.y[sigix],method=cor.method,exact=F)$p.value
				cor.msh.neth <- cor.test(xy$coefficient.x[sigix], xy$coefficient.y[sigix],method=cor.method,exact=F)$estimate
			}
		} else {
			p.msh.neth <- NA
			cor.msh.neth <- NA
		}
		
		colnames(gx.stats.mgh)[3:5] <- sprintf('%s.mgh',colnames(gx.stats.mgh)[3:5])
		colnames(gx.stats.msh)[3:5] <- sprintf('%s.msh',colnames(gx.stats.msh)[3:5])
		colnames(gx.stats.neth)[3:5] <- sprintf('%s.neth',colnames(gx.stats.neth)[3:5])

		
		xy <- merge(gx.stats.mgh, gx.stats.msh[,1:5], by=c('Covariate','Taxon'))
		xyz <- merge(xy, gx.stats.neth[,1:5], by=c('Covariate','Taxon'))
		xyz$pvalue.fisher <- apply(xyz[,c('pvalue.mgh','pvalue.msh','pvalue.neth')],1, fishersMethod)
		xyz$qvalue.fisher <- p.adjust(xyz$pvalue.fisher, method='fdr')
		xyz <- xyz[order(xyz$pvalue.fisher),]


		dev.off()

		xyz.all <- rbind(xyz.all,xyz)

		pvals.crosscohort[[snp]] <- c(p.mgh.msh,p.mgh.neth,p.msh.neth)
		names(pvals.crosscohort[[snp]]) <- c('p.mgh.msh','p.mgh.neth','p.msh.neth')
		cors.crosscohort[[snp]] <- c(cor.mgh.msh,cor.mgh.neth,cor.msh.neth)
		names(cors.crosscohort[[snp]]) <- c('cor.mgh.msh','cor.mgh.neth','cor.msh.neth')
		print(cbind(round(pvals.crosscohort[[snp]],5), round(cors.crosscohort[[snp]],3))[1:2,])
	}

	# cancel out cases where all three cohorts don't agree?
	for(i in 1:nrow(xyz.all)) {
		signs <- c(sign(xyz.all$coefficient.mgh[i]),
				   sign(xyz.all$coefficient.msh[i]),
				   sign(xyz.all$coefficient.neth[i]))
		if(signs[1] != signs[2] && signs[1] != signs[3]){
			xyz.all$pvalue.fisher[i] <- fishersMethod(xyz.all[i,c('pvalue.msh','pvalue.neth')])
		} else if(signs[2] != signs[1] && signs[2] != signs[3]){
			xyz.all$pvalue.fisher[i] <- fishersMethod(xyz.all[i,c('pvalue.mgh','pvalue.neth')])
		} else if(signs[3] != signs[1] && signs[3] != signs[2]){
			xyz.all$pvalue.fisher[i] <- fishersMethod(xyz.all[i,c('pvalue.mgh','pvalue.msh')])
		}
	}
	xyz.all$qvalue.fisher <- p.adjust(xyz.all$pvalue.fisher,'fdr')
	xyz.all <- xyz.all[order(xyz.all$pvalue.fisher),]
	xyz.all <- xyz.all[,c("Covariate","gene","Taxon","Taxon_abbrev", "pvalue.mgh", "pvalue.msh", "pvalue.neth", "coefficient.mgh", "coefficient.msh", "coefficient.neth", "qvalue.mgh", "qvalue.msh", "qvalue.neth","pvalue.fisher", "qvalue.fisher")]
	xyz.all$qvalue.mgh <- p.adjust(xyz.all$pvalue.mgh,'fdr')
	xyz.all$qvalue.msh <- p.adjust(xyz.all$pvalue.msh,'fdr')
	xyz.all$qvalue.neth <- p.adjust(xyz.all$pvalue.neth,'fdr')
	
	xyz.all$gene <- as.character(xyz.all$gene)
	xyz.all$gene[is.na(xyz.all$gene)] <- 'Unknown'
	
	for(i in 1:length(pvals.crosscohort)) {pvals.crosscohort[[i]][cors.crosscohort[[i]] <= 0] <- 1}
	
	# qvals.combined takes worst qvalue within-cohort for each SNP
# 	qvals.combined <- numeric(length(pvals.crosscohort))
# 	for(i in 1:length(pvals.cross.cohort)){
# 		cors.i <- cors.crosscohort[[i]]
# 		
# 		pvals.crosscohort[[i]][cors.crosscohort[[i]] <= 0] <- 1
# 	}
	qvals.combined <- apply(cbind(p.adjust(sapply(pvals.crosscohort,'[',1),'fdr'), p.adjust(sapply(pvals.crosscohort,'[',2),'fdr')),1,function(xx) if(sum(is.na(xx))==2) NA else min(xx[!is.na(xx)]))
	names(qvals.combined) <- sapply(strsplit(names(qvals.combined),'\\.'),'[',1)

	return(list(res.all=xyz.all,
				pvals.cohort.cor=pvals.crosscohort,
				qvals.combined=qvals.combined,
				cors.cohort.cor=cors.crosscohort))
}