# returns summary statistics of microbe interaction
# for each icsnp of that snp, the best snp (or NA if same),
# and the group of fine mapping snps
"test.fine.mapping" <- function(d){
    gxhd <- fast.read.table('~/drive/research/prism/data/genetic/merged/merge-prism-msh-neth-risk-hd-region-snps.raw',header=T,includes.rownames=T,type='character',verbose=T)
    rownames(gxhd) <- gsub('-','__',gxhd[,1]); gxhd <- gxhd[,-(1:5)]; class(gxhd) <- 'numeric'
    colnames(gxhd) <- sapply(strsplit(colnames(gxhd),'_'),function(xx) paste(xx[1:3],collapse='_'))
    gxhd <- gxhd[rownames(d$map),]
    
    finemap <- read.table('~/drive/research/prism/data/genetic/snp-subsets/fine_mapping_results.txt',sep='\t',head=T)
    rs2uid <- read.table('~/drive/research/prism/data/genetic/snp-subsets/hd-region-snps-rs2uid.txt',sep='\t',row=1)
    
    results <- matrix(NA,nrow(finemap),9)
    colnames(results) <- outer(c('ICSNP','BESTSNP','SNP.GROUP'),c('module','mb','beta'),paste,sep='-')
    null.sets <- list()
    null.sets[['module']] <- load('res.joint.module.null.rdt')
    null.sets[['mb']] <- load('res.joint.mb.null.rdt')
    null.sets[['beta']] <- load('res.joint.beta.null.1000.rdt')
    
    for(i in 1:nrow(finemap)){
        cat(i,'of',nrow(finemap),'\n')
        icsnp <- finemap[i,'ICSNP']
        bestsnp <- finemap[i,'Best.SNP']
        if(bestsnp == '') bestsnp <- NA
        snp.group <- ifelse(is.na(bestsnp) || !(bestsnp %in% colnames(gxhd)),icsnp, bestsnp)
        add.snp.columns <- grep('Additional.snp',colnames(finemap))
        for(j in 1:length(add.snp.columns)){
            if(!is.na(finemap[i,add.snp.columns[j]])){
                if(finemap[i,add.snp.columns[j]] != ''){
                    snp.group <- c(snp.group,finemap[i,add.snp.columns[j]])
                }
            }
        }

        if(length(snp.group) == 1) snp.group <- NA
        
        # translate to uid
        if(!(icsnp %in% rownames(rs2uid))) stop(sprintf('ICSNP %s not in rs2uid\n',icsnp))
        icsnp <- rs2uid[icsnp,1]
        if(!is.na(bestsnp)) bestsnp <- rs2uid[bestsnp,1]
        if(!is.na(snp.group)[1]) snp.group <- rs2uid[snp.group,1]
        
        # print update
        cat(sprintf('icsnp=%s\tbest.snp=%s\tsnp.group=%s\n',
            icsnp,
            bestsnp,
            ifelse(is.na(snp.group)[1],'NA',paste(snp.group,collapse='-'))
        ))

        # run test on all applicable
        if(!(icsnp %in% colnames(gxhd))) stop(sprintf('ICSNP %s not in gxhd\n',icsnp))
        x <- matrix(gxhd[,icsnp],ncol=1)
        colnames(x)[1] <- 'ICSNP'
        if(!is.na(bestsnp)){
            if(!(bestsnp %in% colnames(gxhd))){ 
                cat(sprintf('BESTSNP %s not in gxhd\n',bestsnp))
            } else {
                x <- cbind(x,gxhd[,bestsnp])
                colnames(x)[2] <- 'BESTSNP'
            }
        }
        if(!is.na(snp.group)[1]){
#             if(i == 3) browser()
            if(all(!(snp.group %in% colnames(gxhd)))){
                cat(sprintf('SNPGROUP %s not in gxhd\n',snp.group))
            } else if (any(!(snp.group %in% colnames(gxhd)))) {
                snp.group <- snp.group[snp.group %in% colnames(gxhd)]
            } else {
                x <- cbind(x,rowSums(gxhd[,snp.group,drop=F]))
                colnames(x)[ncol(x)] <- 'SNPGROUP'
                x[x[,'SNPGROUP'] > 2,'SNPGROUP'] <- 2
                print(table(x[,'SNPGROUP']))
            }
        }

        
        target.types <- c('module','mb','beta')
        test.types <- c('linear','linear-nonzero','linear')
        for(k in 1:length(target.types)){
            null.set <- 
            baseix <- (k-1) * 3
            res <- test.wrapper(d,predictor=x,targets=target.types[k], test.type=test.types[k])
#             res <- p.distribution.test.from.association.tables(res, null.sets[[target.tyeps]],
#                         alpha=0.05,test.type='geom');
            res.pvals <- sapply(split(res[[1]]$pvalue, res[[1]][,1]),function(xx) exp(mean(log(xx))))
            print(res.pvals)
            results[i, baseix + 1] <- res.pvals['ICSNP']
            if('BESTSNP' %in% names(res.pvals)){
                results[i, baseix + 2] <- res.pvals['BESTSNP']
            }
            if('SNPGROUP' %in% names(res.pvals)){
                results[i, baseix + 3] <- res.pvals['SNPGROUP']
            }
        }

        # test ic snp
        # test best snp (if different)

        # test snp group (best + additional)

    }
    return(results)
}