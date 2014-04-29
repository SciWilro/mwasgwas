args <- commandArgs(TRUE)
ld <- read.table(args[1],head=T)
newld <- ld
allsnps <- c(ld[,3],ld[,6])
keep <- NULL
to.drop <- NULL
counts <- table(allsnps)
while(nrow(newld) > 0){
    hit <- names(which.max(counts))
    keep <- c(keep, hit)
    # remove all snps matching this hit from table
    # all entries for each matching snp
    direct.matches <- unique(c(newld[newld[,3]==hit,6], newld[newld[,6]==hit,3]))
    to.drop <- c(to.drop, direct.matches)
    
    indirect.matches.ixl <- (newld[,3] %in% direct.matches | newld[,6] %in% direct.matches)
    indirect.match.counts <- table(c(newld[indirect.matches.ixl,3], newld[indirect.matches.ixl,6]))

    counts[names(indirect.match.counts)] <- counts[names(indirect.match.counts)] - indirect.match.counts
    counts <- counts[counts > 0]
#     cat(sprintf('%d direct matches and %d indirect matches removed, %d rows remain.\n',length(direct.matches),sum(indirect.matches.ixl)-length(direct.matches),sum(!indirect.matches.ixl)))
    newld <- newld[!indirect.matches.ixl,,drop=F]
}

cat(to.drop,sep='\n')