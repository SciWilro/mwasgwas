basedir <- Sys.getenv('MWAS_GWAS_DIR')
source(sprintf('%s/lib/rpackage/load2.R',basedir))
IN.TABLE <- '~/drive/research/prism/data/genetic/merged-prism-msh-nl.raw'
OUT.TABLE <- '~/drive/research/prism/data/genetic/merged-prism-msh-nl-with-NOD2.raw'
IN.ANNOT <- '~/drive/research/prism/data/genetic/2013-01-09-generate-snps/all-snps-metadata-with-pos-with-categories.txt'
OUT.ANNOT <- '~/drive/research/prism/data/genetic/2013-01-09-generate-snps/all-snps-metadata-with-pos-with-categories-with-NOD2.txt'

gx <- load.genotype.data(IN.TABLE,min.maf=0, annotations.fp=IN.ANNOT)
nod2 <- rowSums(gx$x[,gx$annot$nod2.snps])
nod2[is.na(nod2)] <- 0
nod2[nod2 > 2] <- 2
X <- data.frame(NOD2=nod2)
X <- cbind(X,gx$x[,gx$annot$jostins.snps])

annot <- gx$annot
annot <- rbind(annot,annot['rs5743293',])
rownames(annot)[nrow(annot)] <- 'NOD2'
annot['NOD2','jostins.snps'] <- TRUE
annot <- annot[colnames(X),]

sink(OUT.TABLE)
cat('Subject\t')
write.table(X,quote=F,sep='\t')
sink(NULL)

sink(OUT.ANNOT)
cat('SNP\t')
write.table(annot, quote=F,sep='\t')
sink(NULL)