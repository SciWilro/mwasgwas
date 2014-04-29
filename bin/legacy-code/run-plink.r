source('~/drive/research/prism/code/load.r')

# run analysis
if(args$verbose) cat('Getting plink pvals...\n') 
res <- get.plink.pvals(mwas, gwas, normtype='asin', verbose=args$verbose)
hitix <- which(res[,'qvalue'] < 0.01)
if(args$verbose) cat(sprintf('%d hits...\n',length(hitix)))

# save result to file
if(args$verbose) cat('Saving results...\n')
write.table(res,file=sprintf('%s/results.txt',args$outdir),sep='\t',quote=F)
# save.mwas.gwas.results(res)

if(length(hitix) > 0){
        
    # print result to screen
    print(res[hitix,])
    # print.mwas.gwas.results(res)
    # make plots of result
    # plot.mwas.gwas.results(res)
    for(i in 1:length(hitix)){
        pdf(sprintf('%s/hit%02d.pdf',args$outdir,i),width=5,height=5)
        par(mar=c(5.1,5.6,4.1,2.1))
        
        stripchart(mwas[,res[hitix[i],2]] ~ gwas$counts[,res[hitix[i],1]],
                method='jitter',
                main=sprintf('%s\nq= %.2g',res[hitix[i],2], res[hitix[i],4]),
                cex.main=.5, vertical=TRUE,
                xlab=sprintf('Risk allele count, %s',res[hitix[i],1]),
                ylab='Relative abundance',
                cex.lab=.7,
                cex.axis=.6,
                col='#000000aa')
        dev.off()
    }

    for(i in 1:length(hitix)){
        pdf(sprintf('%s/boxplot-hit%02d-.pdf',args$outdir,i),width=5,height=5)
        par(mar=c(5.1,5.6,4.1,2.1))
        boxplot(mwas[,res[hitix[i],2]] ~ gwas$counts[,res[hitix[i],1]],
                main=sprintf('%s\nq= %.2g',res[hitix[i],2], res[hitix[i],4]),
                cex.main=.5, horizontal=FALSE,
                xlab=sprintf('Risk allele count, %s',res[hitix[i],1]),
                ylab='Relative abundance',
                cex.lab=.7,
                cex.axis=.6,
                col='#000000aa')
        dev.off()
    }

} else {
    cat('No hits with FDR < 0.1\n')
}


if(!is.null(args$groupfile)){
    group.snps <- scan(args$groupfile,quiet=TRUE,what='character')
    mwas.collapse <- mwas
    print(group.snps)
    print(group.snps %in% colnames(gwas$counts))
    snpix <- match(group.snps[group.snps %in% colnames(gwas$counts)], colnames(gwas$counts))
    print(snpix)
    # collapse taxa
    # genera <- sapply(strsplit(colnames(mwas),';'), '[',6)
    # genera
    # mwas.collapse <- matrix(0,nrow=nrow(mwas), ncol = length(unique(genera)))
    # colnames(mwas.collapse) <- unique(genera)
    # rownames(mwas.collapse) <- rownames(mwas)
    #for(i in 1:length(unique(genera))){mwas.collapse[,i] <- rowSums(mwas[,genera == unique(genera)[i],drop=F])}
    
    group.counts <- rowSums(gwas$counts[,snpix])
    pvals <- apply(mwas.collapse,2,function(xx){
        mlm <- lm(asin(sqrt(xx)) ~ group.counts)
        summary(mlm)[[4]][1,4]
    })
    qvals <- p.adjust(pvals, m='fdr')
    ix <- order(qvals)
    pdf(sprintf('%s/%s_group.pdf',args$outdir,args$groupfile),width=10,height=10)
    par(mfrow=c(3,3))
    for(i in 1:min(c(9,ncol(mwas.collapse)))){
        plot(jitter(group.counts,amount=.1), 
                asin(sqrt(mwas.collapse[,ix[i]])),
                main=sprintf('%s: %.2g',colnames(mwas.collapse)[ix[i]],qvals[ix[i]]),
                cex.main=.4)
    }
    dev.off()
}