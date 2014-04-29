# load data
p <- read.table('~/drive/research/prism/data/genetic/prism/preprocessed.raw',head=T,sep='\t')
m <- read.table('~/drive/research/prism/data/genetic/msh/preprocessed.raw',head=T,sep='\t')
n <- read.table('~/drive/research/prism/data/genetic/nl/preprocessed.raw',head=T,sep='\t',comment='')

rownames(p) <- p[,2]; p <- p[,-(1:6)]
rownames(m) <- m[,2]; m <- m[,-(1:6)]
rownames(n) <- n[,2]; n <- n[,-(1:6)]

cat('WARNING: manually deleting duplicate snp rs7608910 from prism data.\n')
p <- p[,-grep("\\.1", colnames(p))]

all.cols <- unique(c(colnames(m), colnames(p), colnames(n)))
gx <- matrix(NA,nrow(p)+nrow(m)+nrow(n), length(all.cols))
rownames(gx) <- c(rownames(p), rownames(m), rownames(n))
colnames(gx) <- sort(all.cols)
gx[rownames(p),colnames(p)] <- as.matrix(p)
gx[rownames(m),colnames(m)] <- as.matrix(m)
gx[rownames(n),colnames(n)] <- as.matrix(n)

# drop the allele from plink columns
colnames(gx) <- sapply(strsplit(colnames(gx),'_'),function(xx) xx[1])

# save merged table
write.table(gx,
            file='~/drive/research/prism/data/genetic/merged-prism-msh-nl.raw',
            sep='\t',
            quote=F)