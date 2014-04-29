source('~/drive/research/code/r/general-utilities.r')
source('~/drive/research/prism/code/util.load.r')

# usage Rscript opticall.fp > out.ped
# converts opticall file to ped
# by default uses opticall base file name + '.ped'
"opticall2ped" <- function(opticall.fp, output.fp=NULL, verbose=TRUE){
    cat('reading table...')
    o <- fast.read.table(opticall.fp,head=T,includes.rownames=T, type='character', delim='\t')
    cat('\n')
    browser()
    ped <- matrix('0 0',nrow=ncol(o) - 3, ncol=nrow(o) + 6)
    
    for(i in 1:nrow(o)){
        if(verbose && i %% 1000 == 0) cat('.')
        if(o[i,'Alleles'] == 'DI' || o[i,'Alleles'] == 'ID'){
            # insertion/deletion, process separately
            ins <- strsplit(rownames(o)[i],'-')[[1]]
            ins <- rev(ins)[2]
            alleles <- c(ins,'-')
            if(o[i,'Alleles'] == 'DI') alleles <- rev(alleles)
        } else {
            alleles <- strsplit(o[i,'Alleles'],'')[[1]]
        }
        ped[o[i,-(1:3)] == 1,i+6] <- paste(alleles[1], alleles[1], sep=' ')
        ped[o[i,-(1:3)] == 2,i+6] <- paste(alleles[1], alleles[2], sep=' ')
        ped[o[i,-(1:3)] == 3,i+6] <- paste(alleles[2], alleles[2], sep=' ')
    }
    if(verbose && nrow(o) >= 1000) cat('\n')
    
    # header starts with FID IID PAT MAT SEX PHENOTYPE    
    # FID is IID
    ped[,1] <- colnames(o)[-(1:3)]
    # IID
    ped[,2] <- ped[,1]
    #PAT + MAT null
    ped[,3:4] <- '0'
    ped[,5] <- '1'
    ped[,6] <- '-9'
    colnames(ped) <-  c('FID','IID','PAT','MAT','SEX','PHENOTYPE', rownames(o))
    
    # write/print output
    if(!is.null(output.fp)) sink(output.fp)
    cat('#')
    write.table(ped, col=T, row=F, sep='\t', quote=F)
    if(!is.null(output.fp)) sink(NULL)
}
  
args <- commandArgs(T)
opticall2ped(args[1], paste(file.basename(args[1]), '.ped', sep=''))