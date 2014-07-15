# returns a mixwas feature set object
"load.microbiome.data" <- function(
  fp,
  type=c('taxa','otus','functions','betadiv','alphadiv')[1],
  annotations.fp=NULL, # optional annotations with rownames = feature names
  min.prevalence=0.1,       # mb features must be present in at least this % samples
  min.mean.fraction=0,      # mb features must be average at least this fraction
  collapse.at=0,            # collapse features at this level of spearman's correlation
  normalize=FALSE,          # normalize each row to sum to 1
  rarefy.at=NULL,           # rarefy at this depth (NULL/NA to skip rarefying)
  description=NULL          # description, default = type
){
  if(!is.numeric(min.prevalence)) stop('min.prevalence must be a numeric value between 0 and 1')
  if(!is.numeric(min.mean.fraction)) stop('min.mean.fraction must be a numeric value between 0 and 1')
  if(min.prevalence < 0 || min.prevalence > 1) stop('min.prevalence must be a numeric value between 0 and 1')  
  if(min.mean.fraction < 0 || min.mean.fraction > 1) stop('min.mean.fraction must be a numeric value between 0 and 1')
  
  res <- list()
  annot <- list()
  if(type == 'taxa'){
    res$x <- load.qiime.taxon.table(fp)
  } else if(type == 'otus'){
    tmp <- load.qiime.otu.table(fp,include.lineages=TRUE)
    res$x <- tmp$otus
    annot$lineage <- tmp$lineages
  } else if(type == 'functions'){
    res$x <- load.qiime.otu.table(fp)
    # drop "None" from functions
    res$x <- res$x[,!grepl('None|Unclassified',colnames(res$x))]
  } else {
    stop('Unknown data type specified.')
  }
  
  # rarefy
  if(!is.null(rarefy.at)){
    if(!is.na(rarefy.at)) res$x <- rarefy(res$x,rare.depth)
  }
  
  # get normalized version
  normx <- sweep(res$x,1,rowSums(res$x),'/')

  # filter by prevalence
  normx <- normx[,colMeans(res$x > 0) > min.prevalence]
  res$x <- res$x[,colMeans(res$x > 0) > min.prevalence]
  
  # filter by relative abundance
  res$x <- res$x[,colMeans(normx) > min.mean.fraction]

  # normalize
  if(normalize) res$x <- normx
  
  # collapse
  if(!is.na(collapse.at) & !is.null(collapse.at) & collapse.at > 0){
    collapse <- collapse.by.correlation(res$x, min.cor=collapse.at, "mean", verbose=TRUE)
    res$x <- res$x[,collapse$reps]
    res$collapse.groups <- collapse$groups
  } else {
    res$collapse.groups <- NULL
  }
  
  # load annotations if requested
  if(!is.null(annotations.fp)){
    d <- read.table(annotations.fp,sep='\t',head=T,row=1,comment='',stringsAsFactors=FALSE)
    if(!all(rownames(res$x) %in% rownames(d))) {
        stop('Not all input features are present in the annotation file.')
    }
    d <- d[rownames(res$x),]
    annot <- cbind(annot, d)
  }
  # add annotations to return object
  res$annot <- as.data.frame(annot)
  
  if(is.null(description)) description <- type
  res$description <- description
  return(res)
}


"load.genotype.data" <- function(
  fp,
  type=c('snps','pc')[1],
  annotations.fp=NULL, # optional annotations with rownames = feature names
  min.maf=0.1, # e.g. min.maf=0.1 means only keep snps with > 10% minor allele freq
  max.missing=0.1, # e.g. max.missing=0.1 means only keep snps with < 10% NA
  drop.na.method = c('none','by.row','by.column','2D')[1], # how (if) to remove NA's 
  description=NULL          # description, default = type
){
  res <- list()
  annot <- list()

  if(type == 'snps'){
    res$x <- read.table(fp,head=T,row=1,check=F,comment='',sep='\t')
  } else {
    stop('Unknown data type specified.')
  }

  # filter by MAF
  res$x <- res$x[, get.maf(res$x) > min.maf]
  
  # filter by NA
  res$x <- res$x[, colMeans(is.na(res$x)) < max.missing]
  
  # drop NA's
  if(drop.na.method == 'by.row'){
    res$x <- res$x[rowSums(is.na(res$x))==0,,drop=F]
  } else if(drop.na.method == 'by.column'){
    res$x <- res$x[,colSums(is.na(res$x))==0,drop=F]
  } else if(drop.na.method == '2D'){
    drop.ix <- drop.na.2d(res$x,prefer.columns=2)
    res$x <- res$x[!drop.ix$drop.row, !drop.ix$drop.col,drop=F]
  } else if(drop.na.method != 'none'){
    stop('Unknown NA filtering method specified')
  }
  if(nrow(res$x) == 0 || ncol(res$x) == 0){
    stop("All rows and columns filtered due to NA's")
  }
  
  # convert to principal components if requested
  if(type=='pc'){
    pc <- replicate(10,{
      tmp <- res$x
      tmp[is.na(tmp)] <- 0
      princomp(tmp[,sample(ncol(tmp),size=ceiling(.99*nrow(tmp)))])$scores[,1:3]
    })
    pc <- apply(pc,1:2,mean)
    res$x <- pc
    colnames(res$x) <- sprintf('PC%02d',1:ncol(pc))
  } 
  
  # load annotations if requested
  if(!is.null(annotations.fp)){
    d <- read.table(annotations.fp,sep='\t',head=T,row=1,comment='',stringsAsFactors=FALSE)
    if(!all(colnames(res$x) %in% rownames(d))) {
      stop('Not all input features are present in the annotation file.')
    }
    d <- d[colnames(res$x),]
    if(length(annot) > 0){
      annot <- cbind(annot, d)      
    } else {
      annot <- d
    }
  }

  # add annotations to return object
  res$annot <- as.data.frame(annot)
  
  if(is.null(description)) description <- type
  res$description <- description
  return(res)
}

# simple wrapper to load QIIME mapping file
"load.metadata" <- function(fp){
  
  return(list(x=load.qiime.mapping.file(fp)))
}

# returns MAF for each column
# input is a table of subjects x SNP counts (0,1,2)
"get.maf" <- function(gx){
  if(is.null(dim(gx))) gx <- as.matrix(gx)
  return (.5 - abs(apply(gx,2,function(xx) mean(drop.na(xx))/2)-.5))
}

# drop columns with zero variance?
"drop.empty.mwas.features" <- function(mwas){
    mwas <- mwas[,apply(mwas, 2, var) > 0]
    return(mwas)
}

# drop columns with zero variance
"drop.empty.gwas.features" <- function(gwas){
    variances <- apply(gwas$counts,2,var)
    zero.variance.ix <- which(variances==0 | is.na(variances))
    if(length(zero.variance.ix) > 0){
        gwas$counts <- gwas$counts[,-zero.variance.ix]
        gwas$alleles <- gwas$alleles[-zero.variance.ix]
        gwas$ped <- gwas$ped[,-(zero.variance.ix + 6)]
        gwas$map <- gwas$map[-zero.variance.ix,]
    }
    return(gwas)
}

# reads a QIIME otu/metadata/taxon/distance table.
# Support legacy formats, where
# the header may or may not start with '#', 
# and comment lines can be anywhere in the file.
# return value is a matrix unless as.data.frame is TRUE
"read.qiime.table" <- function(filepath, as.data.frame=FALSE){
    header.index <- get.header.index(filepath)
    # read the header
    f <- file(filepath,'r')
    header <- scan(filepath, what='character', sep='\t',comment='',skip=header.index-1,quote='"',
                    nlines=1,quiet=TRUE)
    close(f)
    # read the rest of the table
    datatable <- read.table(filepath,sep='\t',skip=header.index, comment='#',quote='"',
                        head=F,row.names=1,check=FALSE)
    
    # set column names using header
    colnames(datatable) <- header[-1]
    
    if(!as.data.frame) datatable <- as.matrix(datatable)
    return(datatable)
}


"load.qiime.mapping.file" <- function(filepath){
    return(read.qiime.table(filepath, as.data.frame=TRUE))
}

"load.qiime.otu.table" <- function(filepath,include.lineages=FALSE){
    otus <- read.qiime.table(filepath, as.data.frame=TRUE)

    # drop "Consensus Lineage" column if present
    if(otu.table.has.metadata(colnames(otus))){
        C <- ncol(otus)
        lineages <- as.character(otus[,C])
        otus <- otus[,-C]
    } else {
        lineages <- NULL
    }
    otus <- as.matrix(t(otus))
    
    if(include.lineages){
        return(list(otus=otus,lineages=lineages))
    } else {
        return(otus=otus)
    }
}

# TRUE if last column is "Consensus Lineage" or "OTU Metadata" or "Taxonomy"
# case-insensitive
"otu.table.has.metadata" <- function(headers){
    C <- length(headers)
    has.metadata <- grepl('consensus[ ]lineage|otu[ ]*metadata|taxonomy|KEGG_Pathways',
                          headers[C], ignore.case=TRUE)
    return(has.metadata)
}

# returns the index of the header line
# note: lines after the header may be comments with '#'
# read.table should therefore be called with (skip=header.index, comment='#')
"get.header.index" <- function(filepath){
    ncolumns.per.line <- NULL
    
    # read lines until the first line without a '#'
    # for each line, obtain the number of tab-delimited columns
    linecount <- 0
    start.character <- '#'
    while(start.character == '#'){
        linecount <- linecount + 1
        f <- file(filepath,'r') # open file in read mode
        line <- scan(f,what='character',skip=linecount-1,nlines=1, sep='\t', quiet=TRUE)
        close(f)
        # ncolumns is the number of entries in this line
        # not including trailing empties
        ncolumns <- max(which(sapply(line,nchar) > 0))
        ncolumns <- max(ncolumns, ncolumns.per.line)
        ncolumns.per.line <- c(ncolumns.per.line, ncolumns)
        start.character <- substring(line[1],1,1)
    }
    
    # first non-comment line gives the number of columns
    C <- ncolumns.per.line[linecount]
    if(linecount == 1){
        # if there are no comment lines, then the first line is the header
        header.index <- 1
    } else {
        if(any(ncolumns.per.line[-linecount] == C)){
            # if there is a comment line with the correct number of columns,
            # it is the header
            header.index <- max(which(ncolumns.per.line[-linecount] == C))
        } else {
            # if there is no comment line with the correct number of columns,
            # the first non-comment line is the header
            header.index <- linecount
        }
    }

    return(header.index)
}


# loads qiime-like taxon table
# extracts only given samples/taxa if requested.
# taxon.subset/sample.subset can either be NULL, a character vector, or a filename
"load.taxon.table" <- function(filepath, taxon.subset=NULL, sample.subset=NULL){
    taxa <- load.qiime.taxon.table(filepath)
    
    # extract only requested taxa
    if(!is.null(taxon.subset)){
        if(file.exists(taxon.subset[1])) taxon.subset <- scan(taxon.subset,'',quiet=T)
        taxon.subset <- taxon.subset[taxon.subset %in% colnames(taxa)]
        if(length(taxon.subset) == 0){
            stop('No taxa in subset are present in table.')
        }
        taxa <- taxa[,taxon.subset,drop=F]
    }
    # extract only requested taxa
    if(!is.null(sample.subset)){
        if(file.exists(sample.subset[1])) sample.subset <- scan(sample.subset,'',quiet=T)
        sample.subset <- sample.subset[sample.subset %in% rownames(taxa)]
        if(length(sample.subset) == 0){
            stop('No samples in subset are present in table.')
        }
        taxa <- taxa[sample.subset,,drop=F]
    }
    return(taxa)
}



"load.qiime.taxon.table" <- function(filepath){
    taxa <- as.matrix(t(read.table(filepath,sep='\t',head=T,row.names=1,check=FALSE,quote='"')))
    return(taxa)
}

"load.qiime.distance.matrix" <- function(filepath){
    d <- as.matrix(read.table(filepath,sep='\t',head=T,row.names=1,check=FALSE,quote='"'))
    return(d)
}

# ensure map, data table, etc., contain the same samples in the same order
# option sample subset takes a further intersection
"remove.nonoverlapping.samples" <- function(gwas, mwas, metadata, sample.subset){
    ix <- intersect(rownames(mwas), rownames(gwas$ped))
    ix <- intersect(ix, rownames(metadata))
    ix <- intersect(ix, sample.subset)
    mwas <- mwas[ix,]
    gwas$ped <- gwas$ped[ix,]
    gwas$counts <- gwas$counts[ix,]
    metadata <- metadata[ix,]
    return(list(mwas=mwas,
                gwas=gwas,
                metadata=metadata))
}



# gene.ids is a list of gene ids for each snp
# returns pathways.snp=pathway x snp, pathways.gene = pathway x gene
"pathway.membership" <- function(gene.ids,
        pathway.fp='~/drive/usr/ref/MSigDB/c2.cp.v3.1.symbols-reactome-biocarta-kegg.gmt',
        verbose=FALSE){
    # build snp by gene matrix
    uniq.genes <- sort(unique(unlist(gene.ids)))
    snp.by.gene <- matrix(FALSE,nrow=length(gene.ids), ncol=length(uniq.genes))
    rownames(snp.by.gene) <- names(gene.ids)
    colnames(snp.by.gene) <- uniq.genes
    for(j in 1:length(gene.ids)){
        snp.by.gene[j,uniq.genes %in% gene.ids[[j]]] <- TRUE
    }

    # load snp by pathway matrix
    con = file(pathway.fp,'rt')
    grow.by <- 1000
    pathways.snp <- matrix(FALSE,nrow=grow.by,ncol=length(gene.ids))
    pathways.gene <- matrix(FALSE,nrow=grow.by,ncol=length(uniq.genes))
    pathway.sizes <- NULL
    rownames(pathways.snp) <- rep('',nrow(pathways.snp))
    colnames(pathways.snp) <- names(gene.ids)
    rownames(pathways.gene) <- rep('',nrow(pathways.gene))
    colnames(pathways.gene) <- uniq.genes
    i <- 1
    if(verbose) cat('Loading pathways')
    while(length(line <- scan(con, what='character',nlines=1,sep='\t',quiet=TRUE)) > 0){
        if(i %% 100 == 0 & verbose) cat('.')
        # add more rows to pathway mat if needed
        if(i > nrow(pathways.snp)){
            pathways.snp <- rbind(pathways.snp, matrix(FALSE,nrow=grow.by, ncol=ncol(pathways.snp)))
            pathways.gene <- rbind(pathways.gene, matrix(FALSE,nrow=grow.by, ncol=ncol(pathways.gene)))
        }
        # set membership for each snp based on its genes
        for(j in 1:length(gene.ids)){
            pathways.snp[i,j] <- any(snp.by.gene[j,] & uniq.genes %in% line[-(1:2)])
        }
        
        # set membership for each gene
        pathways.gene[i,colnames(pathways.gene) %in% line[-(1:2)]] <- TRUE
        if(any(pathways.snp[i,])){
            rownames(pathways.snp)[i] <- line[1]
            rownames(pathways.gene)[i] <- line[1]
            pathway.sizes <- c(pathway.sizes,length(line) - 2)
            names(pathway.sizes)[i] <- line[1]
            i <- i + 1
        }
    }
    close(con)
    if(verbose) cat('\n')
    pathways.snp <- pathways.snp[rowSums(pathways.snp) > 0,]
    pathways.gene <- pathways.gene[rowSums(pathways.gene) > 0,]
    return(list(pathways.snp=pathways.snp,pathways.gene=pathways.gene,pathway.sizes=pathway.sizes))
}


# pick a single sample per subject. Note this function is very specific to 
# the particular data set for which it was developed.
# It assumes that some subjects have multiple samples from different biopsy sites
# and some inflamed/non-inflamed.
"pick.unique.samples" <- function(map,subject.column='Gx_Subject_ID',
		sample.site.column='Biopsy_Location_General',
		sample.site.order=c('Ileum','Colon','Rectum','Pouch','Pre_Pouch_Ileum'),
		inflamed.column='Inflamed',
		inflamed.order=c('Yes','No'),
		return.logical=TRUE){
	keep.ix <- NULL
	subj.ids <- split(rownames(m),m[,subject.column])
	for(i in seq_along(subj.ids)){
		subj.ix <- match(subj.ids[[i]],rownames(m))
		if(length(subj.ix) == 1) {
			keep.ix <- c(keep.ix,subj.ix)
		} else {
			infl <- m[subj.ix,inflamed.column]
			loc <- m[subj.ix,sample.site.column]
			keep.ix.i <- NULL
			for(loc.type in sample.site.order){
				for(infl.type in inflamed.order){
					ix.i <- infl == infl.type & loc == loc.type
					if(any(ix.i)) keep.ix.i <- c(keep.ix.i, subj.ix[ix.i][1])
				}
			}
			keep.ix.i <- c(keep.ix.i,subj.ix)
			keep.ix <- c(keep.ix,keep.ix.i[1])
		}
	}
	keep.ix.logical <- rep(FALSE,nrow(map))
	keep.ix.logical[keep.ix] <- TRUE
	if(return.logical) return(keep.ix.logical)
	return(keep.ix)
}