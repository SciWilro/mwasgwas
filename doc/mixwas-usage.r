taxa <- load.microbiome.data(fp,type='taxa')
pc <- load.qiime.
pathways <- load.qiime.taxon()

gx <- loadgx
null.gx <- loadgx
md <- loadmetadata

mb <- list(taxa, pc, pathways, ...)[[1]]
mx <- mxwas.dataset(mb, gx, sample.md, taxa.md, x.md)
mx <- impute.missing.data(mb, gx)

res.pw <- pairwise.tests(mx, covariates=c('Age','Inflamed'), subset=
				mx.transform=,
				gx.transform=,
				test.type=,
				...)
				
res.snp <- snp.enrichment.tests(res.pw, )
res.path <- pathway.enrichment.tests(res.pw, )

# print/save output
top.hits(res.pw, n=NA, alpha=NA, p.adjust='fdr', file=NULL)
top.hits(res.snp, n=NA, alpha=NA, p.adjust='fdr', file=NULL)
top.hits(res.path, n=NA, alpha=NA, p.adjust='fdr', file=NULL)

# plot output
plot(res.pw)
plot(res.snp)
plot(res.path)
