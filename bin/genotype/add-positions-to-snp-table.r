ichipmap.37 <- read.table('~/drive/research/prism/data/genetic/prism/2012-11-08-convert-SNPs/ichip2dbsnp_map-full.txt',sep='\t',head=T,check=F,comment='')
ichipmap.36 <- read.table('~/drive/research/prism/data/genetic/prism/2012-11-08-convert-SNPs/ichipmap-full-ref36.txt',sep='\t',head=T,check=F,comment='')

snps <- read.table('all-snps-metadata.txt',sep='\t',head=T,row=1,check=F)
ichipmap.36 <- ichipmap.36[ichipmap.36$distance <= 5,]
ichipmap.36 <- ichipmap.36[order(ichipmap.36$distance),]
ichipmap.37 <- ichipmap.37[ichipmap.37$distance <= 5,]
ichipmap.37 <- ichipmap.37[order(ichipmap.37$distance),]


chr.36 <- ichipmap.36$chromosome[match(rownames(snps),ichipmap.36$rsID)]
pos.36 <- ichipmap.36$ichip_build37_pos[match(rownames(snps),ichipmap.36$rsID)]
chr.37 <- ichipmap.37$chromosome[match(rownames(snps),ichipmap.37$rsID)]
pos.37 <- ichipmap.37$ichip_build37_pos[match(rownames(snps),ichipmap.37$rsID)]
snps$chr <- chr.37
snps$pos <- pos.37
snps$chr36 <- chr.36
snps$pos36 <- pos.36

sink('all-snps-metadata-with-pos.txt'); cat('IC_SNP\t'); write.table(snps,quot=F,sep='\t'); sink(NULL)



