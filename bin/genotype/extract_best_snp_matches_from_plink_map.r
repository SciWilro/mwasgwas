# run with
# Rscript extract_best.... snptable plink.map > snp-ichip-mapping.txt
argv <- commandArgs(TRUE)
snp.table.fp <- argv[1]
plink.map.fp <- argv[2]
m <- read.table(plink.map.fp,sep='\t',head=F)
snps <- read.table(snp.table.fp,sep='\t',head=T,row=1,check=F)

# first match by rsID
matches <- NULL

direct.match.rsID <- rownames(snps)[rownames(snps) %in% m[,2]]
matches <- rbind(matches, cbind(direct.match.rsID, direct.match.rsID))
snps <- snps[setdiff(rownames(snps), direct.match.rsID),]

# match remaining by position in build 37
mapnames <- paste(as.character(m[,1]), as.character(m[,4]),sep='__')
snpnames <- paste(as.character(snps$chr),as.character(snps$pos),sep='__')
match.by.pos <- snpnames[snpnames %in% mapnames]
rsIDs <- rownames(snps)[match(match.by.pos, snpnames)]
icIDs <- m[match(match.by.pos, mapnames),2]
matches <- rbind(matches, cbind(rsIDs, icIDs))

# match remaining by position in build 36
mapnames <- paste(as.character(m[,1]), as.character(m[,4]),sep='__')
snpnames <- paste(as.character(snps$chr36),as.character(snps$pos36),sep='__')
match.by.pos <- snpnames[snpnames %in% mapnames]
rsIDs <- rownames(snps)[match(match.by.pos, snpnames)]
icIDs <- m[match(match.by.pos, mapnames),2]
matches <- rbind(matches, cbind(rsIDs, icIDs))

write.table(matches[,2:1],sep='\t',quote=F,col=F,row=F)