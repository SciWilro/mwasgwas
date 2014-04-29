"is.outlier" <- function(x,range=1.5){
	# normalize first
	x <- (x - mean(x)) / sd(x)
    sx <- summary(x)
    q1 <- sx[2]
    q3 <- sx[5]
    iqr <- q3 - q1
    ret <- x > q3  + range * iqr | x < q1 - range * iqr
    if(any(is.na(ret))) ret <- rep(FALSE,length(x))
    if(all(ret)) ret <- !ret # don't call all outliers
    return(ret)
}

"drop.na" <- function(x){
 return(x[!is.na(x)])
}

# drops file extension
"file.basename" <- function(fp){
    fp <- strsplit(fp,'\\.')[[1]]
    fp <- paste(fp[1:max(1,(length(fp)-1))],collapse='.')
    return(fp)
}

# returns indices of those points > range * IQR beyond the 1st/3rd quartile
"n" <- function(x,range=1.5) {
    if(is.null(range) || is.na(range)){
        return(rep(TRUE,length(x)))
    } else {
        quartiles <- boxplot(x,range=range,plot=FALSE)$stats[c(1,5)]
        return(x < quartiles[1] | x > quartiles[2])
    }
}

# pooled variance of x split by y
"pooled.variance" <- function(x, y){
    v <- sapply(split(x, y), var)
    len <- sapply(split(x, y),length) - 1
    return(sum(v * len) / sum(len))
}

# pseudo-log, using small eps to handle zeros
# if eps is null, uses 1/2 of smallest nonzero value
"pseudolog" <- function(x,eps=NULL,log.base=10){
    if(is.null(eps)){
        eps <- min(x[x>0])/2
    }
    x[x==0] <- eps
    return(log(x,base=log.base))
}

# note: passing a column.subset list of columns
# causes this to read line-by-line, extracting just the columns of interest
"fast.read.table" <- function(filepath,header=TRUE, includes.rownames=TRUE,
        column.subset=NULL, type=c('character','numeric')[2], verbose=FALSE,
        nlines=NULL, delim='\t'
    ){
    if(!header && !is.null(column.subset)) stop('Must have header if using column subset.')
    
    # read past header line
    if(header){
        lines <- list()
        count <- 0
        while(count == 0 || substring(lines[[count]],1,1) == '#'){
            count <- count + 1
            lines[[count]] <- scan(filepath,what='character',sep='\n',multi.line=FALSE, nlines=1,skip=count-1, quiet=TRUE);
        }
        if(count == 1){
            count <- 2
        }
        headers <- strsplit(lines[[count - 1]],delim)[[1]]
        skip.count <- count - 1
    } else {
        skip.count <- 0
        headers <- scan(filepath,what='character',sep=delim,multi.line=FALSE, nlines=1, quiet=TRUE);
    }
    # skip first column if includes.rownames
    if(includes.rownames) headers <- headers[-1]
    if(is.null(column.subset)){
        keepix <- 1:length(headers)
    } else {
        overlap <- intersect(column.subset, headers)
        keepix <- match(overlap,headers)
    }

    if(verbose) cat('Fast reading file', filepath,'\n')
    numcolumns <- length(keepix) + includes.rownames
    if(length(keepix) == length(headers)){
        if(verbose) cat('Reading data...\n')
        if(is.null(nlines)) nlines <- numlines <- nlines(filepath) - skip.count
        dd <- scan(filepath,what='character',sep=delim,multi.line=FALSE,
                    skip=skip.count, quiet=TRUE, nlines=nlines)
    } else {
        if(verbose) cat('Counting number of lines...')
        numlines <- nlines(filepath) - skip.count
        if(verbose) cat(sprintf('%d line(s) found.\n',numlines))
        if(verbose) cat('Reading line by line...\n')
        dd <- character(numcolumns * numlines) 
        
        con = file(filepath,'rt')
        linecount <- 0
        if(skip.count > 0){
            invisible(readLines(con,n=skip.count))
        }
        line <- scan(con, what='character',nlines=1,sep='\t',quiet=TRUE)

        while(length(line) > 0 && linecount < nlines){
            if(includes.rownames){
                line <- line[c(1,1+keepix)]
            } else {
                line <- line[keepix]
            }
            ddix <- 1:numcolumns + linecount * numcolumns
            dd[ddix] <- line
            line <- scan(con, what='character',nlines=1,sep='\t',quiet=TRUE)
            linecount <- linecount + 1
        }
        close(con)
    }

    if(verbose) cat('Casting to matrix...\n')
    ddf <- matrix(dd,ncol=numcolumns, byrow=TRUE)
    if(header){
        if(includes.rownames){
            colnames(ddf) <- c('rownames',headers[keepix])
        } else {
            colnames(ddf) <- headers[keepix]
        }
    } else {
        colnames(ddf) <- sprintf('Column_%d',1:length(headers))
    }
    if(includes.rownames){
        rownames(ddf) <- ddf[,1]
        ddf <- ddf[,-1]
    }
    if(type=='numeric'){
        if(verbose) cat('Casting to numeric...\n')
        rn <- rownames(ddf)
        ddf <- apply(ddf,2,as.numeric)
        rownames(ddf) <- rn
    }
    return(ddf)
}

# returns the number of lines in a file
"nlines" <- function(filepath){
    x <- strsplit(system(sprintf('wc -l %s',filepath),intern=TRUE),' ')[[1]]
    return(as.numeric(x[min(which(sapply(x,nchar) > 0))]))
}

"startswith" <- function(string, pattern){
    if(length(string) > 1){
        return(sapply(string, startswith, pattern))
    }
    if(nchar(pattern) > nchar(string)) return(FALSE)
    return(substr(string,1,nchar(pattern)) == pattern)
}

"endswith" <- function(string, pattern){
    if(length(string) > 1){
        return(sapply(string, endswith, pattern))
    }
    if(nchar(pattern) > nchar(string)) return(FALSE)
    return(substr(string,nchar(string)-nchar(pattern)+1,nchar(string)) == pattern)
}

# extracts only samples that match the given metadata filter
# filter is ;-separated, e.g.:
# 'GENDER:Male;COUNTRY:*,!USA'
#
# metadata is a table
"extract.samples.by.metadata" <- function(metadata, states){
    filter.strings <- strsplit(states,';')[[1]]
    headers <- as.character(sapply(sapply(filter.strings, strsplit, ':'),'[',1))
    values <- as.character(sapply(sapply(filter.strings, strsplit, ':'),'[',2))
    filters <- sapply(values, strsplit, ',')
    names(filters) <- headers
    keep.sample <- rep(TRUE, nrow(metadata))
    for(header in names(filters)){
        values <- filters[[header]]
        keep.sample.i <- rep(FALSE, nrow(metadata))
        if(is.element('*',values)){
            keep.sample.i <- rep(TRUE, nrow(metadata))
        }
        for(value in values){
            if(value != '*'){
                # get all matching indices                
                if(substring(value, 1, 1) == '!'){
                    ix <- metadata[[header]] == substring(value,2)
                    keep.sample.i[ix] <- FALSE
                } else {
                    ix <- metadata[[header]] == value
                    keep.sample.i[ix] <- TRUE
                }
            }
        }
        keep.sample <- keep.sample & keep.sample.i
    }
    return(rownames(metadata)[keep.sample])
}

resetPar <- function() {
    dev.new()
    op <- par(no.readonly = TRUE)
    dev.off()
    op
}

# plots an ROC curve from vector of TRUE probabilities
# (only works for TRUE/FALSE classification)
# 
# y has the true response values - must be TRUE/FALSE
#
# p must be vector of probabilities of TRUE, or
# matrix where each column is a different receiver
#
# if p is a matrix, then a separate ROC curve is plotted for each column
# of probabilities
"p.ROC" <- function(p,y,step.size=0.01,filepath='ROC_curve.pdf',cols=NULL){
    if(is.null(cols)) cols <- c('black',brewer.pal(9,'Set1')[-6])
    if(is.null(dim(p))) p <- as.matrix(p,ncol=1)
    if(is.null(colnames(p))) colnames(p) <- sprintf('Classifier %d',1:ncol(p))
    fprs <- list(ncol(p))
    tprs <- list(ncol(p))
    for(i in 1:ncol(p)){
        fprs[[i]] <- sapply(seq(1,0,-abs(step.size)), function(xx) mean(p[!y,i] >= xx))
        tprs[[i]] <- sapply(seq(1,0,-abs(step.size)), function(xx) mean(p[y,i] >= xx))
        tprs[[i]][fprs[[i]] > tprs[[i]]] <- fprs[[i]][fprs[[i]] > tprs[[i]]]
    }

	pdf(filepath,width=4,height=4)
	par(mar=c(5,4,3,1))
	plot(fprs[[1]],tprs[[1]],xlim=c(0,1),ylim=c(0,1),type='l',lwd=2,
			xlab='False positive rate',ylab='True positive rate',col=cols[1])
	abline(0,1,lty=2,col='#000000aa')
	grid()

    if(ncol(p) > 1){
        for(i in 2:ncol(p)){
            lines(fprs[[i]],tprs[[i]],xlim=c(0,1),ylim=c(0,1),type='l',lwd=2,col=cols[i-1])
        }
    }
	legend('topleft',colnames(p),lwd=2,col=cols,cex=.75)
	dev.off()
}

"kld" <- function(p,q,eps=0,normalize=FALSE,use.log2=TRUE){
    if(eps == 0){
        nonzero <- p>0 & q>0
        p <- p[nonzero]
        q <- q[nonzero]
    } 
    else {
        p[p==0] <- eps
        q[q==0] <- eps
        if(normalize){
            p <- p / sum(p)
            q <- q / sum(q)
        }
    }
    if(use.log2) return(sum(p * (log2(p) - log2(q))))
    return(sum(p * (log(p) - log(q))))
}


"jsd" <- function(p,q,eps=0,normalize=FALSE, use.log2=TRUE){
    if(normalize){
        p <- p / sum(p)
        q <- q / sum(q)
    }
    if(eps == 0){
        m <- (p + q)/2
        return((kld(p,m,use.log2=use.log2) + kld(q,m,use.log2=use.log2))/2)
    } else {
        p[p==0] <- eps
        q[q==0] <- eps
        if(normalize){
            p <- p / sum(p)
            q <- q / sum(q)
        }
        m <- (p + q)/2
        return((kld(p,m,use.log2=use.log2) + kld(q,m,use.log2=use.log2))/2)
    }
}

# returns a matrix of pairwise Jensen-Shannon distances between samples
"jsdist" <- function(x,eps=0.000001,normalize=TRUE,use.log2=TRUE){
    N <- nrow(x)
    d <- matrix(0, N, N)

    if(normalize){
        x <- sweep(x,1,rowSums(x),'/')
    }

    for(i in 1:(N-1)){
        for(j in (i+1):N){
            d[i,j] <- d[j,i] <- jsd(x[i,], x[j,],eps=eps,use.log2=use.log2)
        }
    }
    rownames(d) <- rownames(x)
    colnames(d) <- rownames(x)
    return(as.dist(d))
}

# spiral-based colors lifted from bong wang's http://www.nature.com/nmeth/journal/v7/n8/pdf/nmeth0810-573.pdf
require('grDevices',quiet=TRUE)
bong.cols <- rgb(c(207,129,100,125,126,56),c(184,173,157,102,39,0),c(219,211,44,28,18,36),maxColorValue=255)
