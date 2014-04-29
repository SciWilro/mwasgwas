# runs plink externally on each mwas column
# qtr is a quant trait matrix
# plink.file is a file path to the basename of plink map/ped
# or if plink.file is null, use plink.ped and plink.map
# or plink.tped and plink.tfam
# if FID is NULL, assumes fid is same as iid
"get.plink.pvals" <- function(qtr, plink.file=NULL,
        plink.ped=NULL, plink.map=NULL,
        plink.tped=NULL, plink.tfam=NULL,
        covariates=NULL, snp.subset=NULL,
        normtype = c("none", "log", "asin", "boxcox")[2],
        fid=NULL,     print.results=TRUE,
        verbose=FALSE, drop.outliers.range=NULL,
        hwe=.001,maf=0.05,
        geno=.05,mind=0.02, # minimum genotyping rate for snps and individuals
        pfilter=1e-03,
        prefix='r-plink-out',
        remove.intermediate.files=TRUE
    ){
    # process link input file names
    plink.pedflag <- '--ped'
    plink.mapflag <- '--map'
    if(!is.null(plink.file)){
        plink.ped <- sprintf('%s.ped',plink.file)
        plink.map <- sprintf('%s.map',plink.file)
    } else {        
        if(!is.null(plink.tped)){
            plink.ped <- plink.tped
            plink.map <- plink.tfam
            plink.pedflag <- '--tped'
            plink.mapflag <- '--tfam'
        }
    }
    
    if(is.null(dim(qtr))) qtr <- matrix(qtr, ncol=1)
    mm <- NULL
    if(!is.null(covariates)){
        if(nrow(qtr) != nrow(covariates)){
            stop(sprintf('Error: qtr has %d rows, covariates has %d rows.\n',
                nrow(qtr), nrow(covariates)))
        }
        mm <- covariate.model.matrix(covariates)
        mm <- mm[,apply(mm,2,var) > 0,drop=F]
    }

    if(normtype == "log"){
        qtr <- qtr + min(qtr[qtr>0])/2
        qtr <- log(qtr)
    } else if(normtype == "asin"){
        qtr <- asin(sqrt(qtr))
    } else if (normtype == 'boxcox'){
        qtr <- power.transform(qtr)

#         require('geoR')
#         require('car')
#         res <- apply(qtr, 2, function(xx) boxcoxfit(xx,lambda2=min(xx[xx>0])/2))
#         lambda1 <- sapply(res,function(xx) xx$lambda[1])
#         lambda2 <- sapply(res,function(xx) xx$lambda[2])
#         qtr <- sweep(qtr, 2, lambda2,'+')
#         qtr <- bcPower(qtr, lambda1)
    }

    if(is.null(colnames(qtr))) colnames(qtr) <- sprintf('col%d',1:ncol(qtr))

    M <- ncol(qtr)
    grow.by <- 1000
    res <- data.frame(gwasID=rep('',grow.by),
                      qtrID=rep('',grow.by),
                      pvalue=numeric(grow.by),
                      qvalue=numeric(grow.by),
                      stringsAsFactors=FALSE)
    res[,c('pvalue','qvalue')] <- NA
    if(verbose) cat('Of ',M,': ',sep='')
    
    # get plink pvals
    # plink.pvals will hold list of pvals for each taxon
    # run plink externally on each taxon
    plink.pvals <- list()
    outlier <- rep(FALSE,nrow(qtr))
    for(i in 1:ncol(qtr)){
        if(verbose) cat(sprintf('%d of %d\t',i, ncol(qtr)))
        if(verbose && i %% 100 == 0) cat(sprintf('%d of %d\t',i, ncol(qtr)))

        if(!is.null(drop.outliers.range)) outlier <- is.outlier(qtr[,i],range=drop.outliers.range)

        phenovals <- qtr[!outlier,i]
        # make non-negative and non-zero
        phenovals <- phenovals - min(phenovals) + .1
            
        plink.filenames <- write.plink.files(phenovals=phenovals, covar=mm[!outlier,,drop=F], fid=fid, prefix=prefix)
        cmd <- sprintf("plink --noweb %s %s %s %s --pheno %s --linear --out %s --pfilter %f", 
                plink.pedflag, plink.ped,
                plink.mapflag, plink.map,
                plink.filenames$phenofile, prefix, pfilter)
        if(!is.null(mm)) cmd <- sprintf('%s --covar %s', cmd, plink.filenames$covarfile)
        if(!is.null(hwe)) cmd <- paste(cmd,sprintf('--hwe %f',hwe),sep=' ')
        if(!is.null(maf)) cmd <- paste(cmd,sprintf('--maf %f',maf),sep=' ')
        if(!is.null(geno)) cmd <- paste(cmd,sprintf('--geno %f',geno),sep=' ')
        if(!is.null(mind)) cmd <- paste(cmd,sprintf('--mind %f',mind),sep=' ')
#          print(cmd)
#         if(i==1) browser()
        system(cmd,ignore.stdout=TRUE)
        res <- read.table(sprintf('%s.assoc.linear',prefix),head=T,strip.white=T,check=F)
        keepix <- res[,'TEST'] == 'ADD'
        tmp <- res[keepix,'P']
        names(tmp) <- res[keepix,'SNP']
        if(!is.null(snp.subset)){
            tmp <- tmp[names(tmp) %in% snp.subset]
        }
        plink.pvals[[i]] <- tmp
        names(plink.pvals)[i] <- colnames(qtr)[i]
        if(verbose) cat(sprintf('min p=%g; min overall p=%g\n',min(tmp), min(unlist(plink.pvals))))
        n.tests <- scan(paste(prefix,'log',sep='.'),what='character',quiet=T)
        n.tests <- as.numeric(n.tests[rev(grep('SNPs',n.tests))[1]-1])
    }
    n.tests <- n.tests * ncol(qtr)
    
    if(verbose) cat('\n')
    to.remove <- c(sprintf('%s.%s',prefix,c('log','pheno','assoc.linear')))
    if(!is.null(mm)) to.remove <- c(to.remove, sprintf('%s.%s',prefix,'covar'))
    if(remove.intermediate.files)  sapply(to.remove, file.remove)
    return(list(pvals=plink.pvals,qvals=p.adjust(unlist(plink.pvals),'fdr',n=n.tests),n.tests=n.tests))

    plink.pvals <- sapply(plink.pvals,function(xx) xx)    
    colnames(plink.pvals) <- colnames(qtr)
    
    # fdr adjustment if requested
    if(verbose) cat('Adjusting p-values...\n')
    plink.qvals <- plink.pvals
    plink.qvals[,] <- p.adjust(plink.qvals,'BH')
    
    # convert to result dataframe with each row = variant, taxon, p-value
    if(verbose) cat('Converting to data frame...\n')
    N <- prod(dim(plink.pvals)) # number of comparisons
    res.mat <- data.frame(variant=character(N), taxon=character(N),
                          pvalue=numeric(N), qvalue=numeric(N), stringsAsFactors=F)
    grid.ix <- expand.grid(1:nrow(plink.pvals), 1:ncol(plink.pvals))
    
    res.mat[,1] <- rownames(plink.pvals)[grid.ix[,1]]
    res.mat[,2] <- colnames(plink.pvals)[grid.ix[,2]]
    res.mat[,3] <- as.numeric(plink.pvals)
    res.mat[,4] <- as.numeric(plink.qvals)
    
    res.mat <- res.mat[order(res.mat[,'qvalue']),]
    
    if(print.results){
        max.print <- 10
        min.print <- 5
        sig <- which(res.mat[,'qvalue'] < .1)
        if(length(sig) > 0){
            print(res.mat[sig[1:max(min.print,min(length(sig),max.print))],,drop=F])
            if(length(sig) > max.print){
                cat(sprintf('+ %d more...\n',length(sig) - max.print))
            }
        } else {
            print(res.mat[1:5,,drop=F])
        }
    }

    
    return(res.mat)
}


# runs plink externally on each qtr column
# mwas is a taxon table
# gwas is a gwas object with at least a $ped and $map object
"get.plink.pvals.old" <- function(mwas, gwas,
        set.file=NULL, normtype = c("none", "log", "asin")[1],
        verbose=FALSE){

    quant <- TRUE
    if(normtype == "log"){
        mwas[mwas==0] <- min(mwas[mwas>0])/2
        mwas <- -log(mwas)
    } else if(normtype == "asin"){
        mwas <- asin(sqrt(mwas))
    }

    # if mwas is a vector, convert to a 1-column matrix
    if(is.null(dim(mwas))){
        tmp <- names(mwas)
        mwas <- matrix(mwas,ncol=1)
        rownames(mwas) <- tmp
        colnames(mwas) <- "MWAS"
    }

    if(identical(sort(rownames(mwas)), sort(rownames(gwas$ped)))){
        mwas <- mwas[rownames(gwas$ped),,drop=FALSE]
    } else {
        stop('MWAS and GWAS tables do not contain the same samples')
    }

    # get plink pvals
    # plink.pvals will hold list of pvals for each taxon
    # run plink externally on each taxon
    plink.pvals <- list()
    for(i in 1:ncol(mwas)){
        if(verbose) cat('.')
        if(verbose && i %% 100 == 0) cat('\n')
        if(quant){
            phenovals <- mwas[,i]
        } else {
            phenovals <- as.numeric(mwas[,i] > 0)
        }
            
        plink.filenames <- write.plink.files(gwas$ped, gwas$map, phenovals)
        
        outfile <- 'r-plink-out'
        cmd <- sprintf("plink --ped %s --map %s --pheno %s --assoc --out %s", 
                plink.filenames$pedfile, plink.filenames$mapfile,
                plink.filenames$phenofile, outfile)
        if(!is.null(set.file)){
            cmd <- paste(sprintf('%s --set %s ',cmd, set.file))
        }
        system(cmd,ignore.stdout=TRUE)
        if(quant){
            res <- read.table(sprintf('%s.qassoc',outfile),head=T,strip.white=T,check=F)
        } else {
            res <- read.table(sprintf('%s.assoc',outfile),head=T,strip.white=T,check=F)
        }
        tmp <- res[,'P']
        names(tmp) <- res[,'SNP']
        plink.pvals[[i]] <- tmp
    }
    if(verbose) cat('\n')
    
    to.remove <- c('r-plink-out.log','r-plink.map','r-plink.ped','r-plink.pheno')
    if(quant){
        to.remove <- c(to.remove,'r-plink-out.qassoc')
    } else {
        to.remove <- c(to.remove,'r-plink-out.assoc')
    }
    sapply(to.remove, file.remove)

    plink.pvals <- sapply(plink.pvals,function(xx) xx)    
    colnames(plink.pvals) <- colnames(mwas)
    
    # fdr adjustment if requested
    if(verbose) cat('Adjusting p-values...\n')
    plink.qvals <- plink.pvals
    plink.qvals[,] <- p.adjust(plink.qvals,'BH')
    
    # convert to result dataframe with each row = variant, taxon, p-value
    if(verbose) cat('Converting to data frame...\n')
    N <- prod(dim(plink.pvals)) # number of comparisons
    res.mat <- data.frame(variant=character(N), taxon=character(N),
                          pvalue=numeric(N), qvalue=numeric(N), stringsAsFactors=F)
    grid.ix <- expand.grid(1:nrow(plink.pvals), 1:ncol(plink.pvals))

    
    res.mat[,1] <- rownames(plink.pvals)[grid.ix[,1]]
    res.mat[,2] <- colnames(plink.pvals)[grid.ix[,2]]
    res.mat[,3] <- as.numeric(plink.pvals)
    res.mat[,4] <- as.numeric(plink.qvals)
    
    res.mat <- res.mat[order(res.mat[,'qvalue']),]
    return(res.mat)
}

    
# load ped file and convert to risk allele counts
# also requires a "map" file corresponding to the ped file
# optionally, variant.subset and sample.subset are files containing
# whitespace-separated lists of variant IDs or sampleIDs to use
#
# assumes that each map entry is a column of ped file
# variant.subset/sample.subset can either be NULL, a character vector, or a filename
#
# returns a list of
# counts: a table of sample x risk allele counts
# ped: the (subset of) ped data as a table
# map: the (subset of) map entries as a table
# major.alleles: a list of the major alleles
# minor.alleles: a list of the minor alleles
"load.ped.file" <- function(ped.filepath, map.filepath,
            sample.subset=NULL, variant.subset=NULL, verbose=FALSE){
    # load gwas data
    system.time(ped <- fast.read.table(ped.filepath, type='character',
            includes.rownames=FALSE, header=FALSE))
    rownames(ped) <- ped[,2]

    # load map file corresponding to ped file
    cat('Loading map...\n')
    map <- read.table(map.filepath,sep='\t',head=F,row=NULL,check=F)
    rownames(map) <- map[,2]
    
    colnames(ped) <- c('Family ID','Individual ID',
                        'Paternal ID','Maternal ID',
                        'Sex','Phenotype', rownames(map))

#     common <- intersect(colnames(ped), rownames(map))
#     map <- map[common,]
#     ped <- ped[,c(colnames(ped)[1:6],common)]

    cat('Extracting sample subset...\n')
    if(!is.null(sample.subset)){
        if(file.exists(sample.subset[1])) sample.subset <- scan(sample.subset,'',quiet=T)
        ix <- intersect(sample.subset, rownames(ped))
#         if(verbose && length(ix) < length(sample.subset)){
#             warning(sprintf(
#                 'Only %d of %d samples in sample subset are present in ped file.',
#                 length(ix),
#                 length(sample.subset)))
#         }
        ped <- ped[ix,]
    }
    
    cat('Extracting variant subset...\n')
    if(!is.null(variant.subset)){
        if(file.exists(variant.subset[1])) variant.subset <- scan(variant.subset,'',quiet=T)
        ix <- intersect(variant.subset, rownames(map))
#         if(verbose && length(ix) < length(variant.subset)){
#             warning(sprintf(
#                 'Only %d of %d variants in variants subset are present in ped file.',
#                 length(ix),
#                 length(variant.subset)))
#         }
         ped <- cbind(ped[,1:6], ped[,ix])
         map <- map[ix,]
    }

    cat('Building allele freqs...\n')
    # build list of letter frequencies in each column
    # we will assume most frequent letter is most commont variant
    # alleles will list most common variant first
    denom <- 2*nrow(ped)
    alleles <- list()
    counts <- matrix(0,nrow=nrow(ped), ncol=ncol(ped)-6)
    if(ncol(ped) > 10000) cat(sprintf('Of %d, \n', ncol(ped)-6))
    for(i in 7:ncol(ped)){
#         if((i-6) %% round(((ncol(ped)-6)/10)) == 0) cat('.')
        if((i-6) %% 10000 == 0) cat(i-6, ' ', sep='')
        allele.pair.types.i <- table(ped[,i])
        alleles.splits <- strsplit(names(allele.pair.types.i),' ')
        allele.types.i <- unique(unlist(alleles.splits))        
        allele.freqs.i <- numeric(length(allele.types.i))
        for(j in 1:length(allele.types.i)){
            allele.freqs.i[j] <- sum(sapply(alleles.splits,function(xx) allele.types.i[j] %in% xx) * allele.pair.types.i)
        }
        names(allele.freqs.i) <- allele.types.i
        allele.freqs.i <- sort(allele.freqs.i,dec=T)/denom
        major.allele <- names(allele.freqs.i)[1]
        allele.pair.values <- sapply(alleles.splits,function(xx) sum(xx != major.allele))
        counts[,i-6] <- allele.pair.values[match(ped[,i], names(allele.pair.types.i))]
        counts[ped[,i] == '0 0',i-6] <- NA
        alleles[[i-6]] <- allele.freqs.i
    }; cat('\n')

#     alleles <- sapply(apply(ped[,-(1:6)],2, function(xx){
#                     table(unlist(strsplit(xx,' ')))/denom
#                 }),
#                 sort,dec=T)
# 
#      cat('Building counts...\n')
#     # build matrix of risk allele counts per subject
#      counts <- sapply(1:(ncol(ped) - 6),
#          function(xx) sapply(strsplit(ped[,-(1:6)][,xx],' '),
#              function(xxx) sum(xxx != names(alleles[[xx]])[1])
#          )
#      )

    colnames(counts) <- rownames(map)
    rownames(counts) <- ped[,2]
    
    return(list(counts=counts, ped=ped, map=map, alleles=alleles))
}



# load ped file and convert to risk allele counts
# also requires a "map" file corresponding to the ped file
# optionally, variant.subset and sample.subset are files containing
# whitespace-separated lists of variant IDs or sampleIDs to use
#
# assumes that each map entry is a column of ped file
# variant.subset/sample.subset can either be NULL, a character vector, or a filename
#
# returns a list of
# counts: a table of sample x risk allele counts
# ped: the (subset of) ped data as a table
# map: the (subset of) map entries as a table
# major.alleles: a list of the major alleles
# minor.alleles: a list of the minor alleles
"load.opticall.file" <- function(ped.filepath,
            sample.subset=NULL, variant.subset=NULL, verbose=FALSE,
            delim=' '){
    # load gwas data
    system.time(ped <- fast.read.table(ped.filepath,type='character',
                    includes.rownames=TRUE, header=TRUE, delim=' '))
    map <- ped[,1:4]
    ped <- t(ped[,-(1:4)])

    cat('Extracting sample subset...\n')
    if(!is.null(sample.subset)){
        if(file.exists(sample.subset[1])) sample.subset <- scan(sample.subset,'',quiet=T)
        ix <- intersect(sample.subset, rownames(ped))
        ped <- ped[ix,]
    }
    
    cat('Extracting variant subset...\n')
    if(!is.null(variant.subset)){
        if(file.exists(variant.subset[1])) variant.subset <- scan(variant.subset,'',quiet=T)
        ix <- intersect(variant.subset, rownames(map))
        ped <- ped[,ix]
        map <- map[ix,]
    }

    cat('Building allele freqs...\n')
    counts <- ped
    counts[counts == '4'] <- NA
    counts <- apply(counts, 2, as.numeric) - 1
    rownames(counts) <- rownames(ped)
    
    allele.names <- sapply(map[,3],function(xx) strsplit(xx,''))
    alleles <- list()
    if(verbose && ncol(ped) > 10000) cat(sprintf('Of %d, \n', ncol(ped)))
    for(i in 1:length(allele.names)){
        if(verbose && i %% 10000 == 0) cat(i, ' ', sep='')
        alleles[[i]] <- numeric(length(allele.names[[i]]))
        names(alleles[[i]]) <- allele.names[[i]]
        alleles[[i]][1] <- sum(2-drop.na(counts[,i]))
        alleles[[i]][2] <- sum(drop.na(counts[,i]))
        alleles[[i]] <- alleles[[i]] / sum(alleles[[i]])
    }
    if(verbose && ncol(ped) > 10000) cat('\n')

    colnames(counts) <- rownames(map)
    rownames(counts) <- rownames(ped)
    
    return(list(counts=counts, ped=ped, map=map, alleles=alleles))
}


# writes ped file, map file, and pheno file
# with <suffix>.ped, <suffix>.map, <suffix>.pheno
# 
# covar should be a model matrix (factors as dummy vars)
# phenovals is a vector of phenotype values with -9 for "missing"
"write.plink.files" <- function(ped=NULL, map=NULL, phenovals=NULL, covar=NULL,
        prefix='r-plink',fid=NULL){
    
    # write ped
    pedfile <- mapfile <- phenofile <- covarfile <- NULL
    if(!is.null(ped)){
        pedfile <- sprintf('%s.ped',prefix)
        sink(pedfile)
        cat('# R-generated ped file\n')
        write.table(ped,sep='\t',quote=F,row.names=FALSE,col.names=TRUE)
        sink(NULL)
    }    

    # write map
    if(!is.null(map)){
        map <- map[colnames(ped)[-c(1:6)],]
        mapfile <- sprintf('%s.map',prefix)
        sink(mapfile)
        write.table(map,sep='\t',quote=F,row.names=FALSE,col.names=FALSE)
        sink(NULL)
    }
        
    # write pheno
    if(!is.null(phenovals)){
        if(is.null(fid)) fid <- names(phenovals)
        phenofile <- sprintf('%s.pheno',prefix)
        sink(phenofile)
        write.table(cbind(fid, names(phenovals),phenovals),
                sep='\t',quote=F,row.names=FALSE,col.names=FALSE)
        sink(NULL)
    }

    # write covariates
    if(!is.null(covar)){
        if(is.null(fid)) fid <- rownames(covar)
        covarfile <- sprintf('%s.covar',prefix)
        covar <- cbind(fid, rownames(covar), covar)
        colnames(covar)[1:2] <- c('FID','IID')
        sink(covarfile)
        write.table(covar,
                sep='\t',quote=F,row.names=FALSE,col.names=TRUE)
        sink(NULL)
    }    

    return(list(pedfile=pedfile, mapfile=mapfile, phenofile=phenofile, covarfile=covarfile))
}

# if snps is NULL, gets minor allele counts for all snps
"extract.snps.from.ped" <- function(plinkbase, snps=NULL){
    if(is.null(snps)){
        # extract all snps in LD with hits
        cmd <- sprintf("time plink --noweb --file %s --out tmp_plink_extract --recodeA --tab", plinkbase)    
    } else {
        # write hits to file for extraction
        sink('tmp_plink_extract.txt');
        cat(unique(snps),sep='\n')
        sink(NULL)
        
        # extract all snps in LD with hits
        cmd <- sprintf("time plink --noweb --file %s --extract tmp_plink_extract.txt --out tmp_plink_extract --recodeA --tab", plinkbase)
    }
    print(cmd)
    system(cmd,ignore.stdout=TRUE)
    
    res <- read.table('tmp_plink_extract.raw',sep='\t',head=T,check=F)
    rownames(res) <- res[,2]
    res <- as.matrix(res[,-(1:6)])
    
    unlink(c('tmp_plink_extract.raw','tmp_plink_extract.log'))
    return(res)
 }