source('~/drive/research/mwas/src/lib/stats.r')
source('~/drive/research/mwas/src/lib/util.load.r')
source('~/drive/research/mwas/src/lib/util.r')
source('~/drive/research/prism/src/lib/rpackage/util.r')
source('~/drive/research/prism/src/lib/rpackage/util.load.r')
source('~/drive/research/prism/src/lib/rpackage/general-utilities.r')
library('RColorBrewer')
library('beeswarm')

min.prevalence <- .1

cat('Loading taxa...\n')
taxon.fp <- 'input.v4/merged-taxa-subset.txt'
taxon.fp <- 'input.v4/merged-taxa.txt'
pathway2.fp <- 'input.v4/kegg_L2.txt'
pathway3.fp <- 'input.v4/kegg_L3.txt'

# taxon.fp <- '/Users/danknights/drive/research/prism/data/microbiome/gg11.2/merged-taxa-min500.txt'

map.fp <- 'input.v4/map-subset-imputed-with-NOD2.txt'
map.fp <- 'input.v4/map-imputed.txt'
map.fp <- 'input.v4/map-merged-no-filter.txt'

otu.fp <- 'input.v4/otutable-min500-with-taxa.txt'
# otu.fp <- '/Users/danknights/drive/research/prism/data/microbiome/gg11.2/otutable-merged-taxa-min500.txt'

m <- read.table(map.fp,sep='\t',head=T,row=1,check=F,comment='')
x <- t(read.table(taxon.fp,sep='\t',head=T,row=1,check=F))

pw2 <- t(read.table(pathway2.fp,sep='\t',head=T,row=1,check=F,comment='',skip=1))
pw3 <- t(read.table(pathway3.fp,sep='\t',head=T,row=1,check=F,comment='',skip=1,quote='"'))
outlier.pathways <- c('Ubiquinone and other terpenoid-quinone biosynthesis','Drug metabolism - other enzymes','Carbon fixation in photosynthetic organisms','Arachidonic acid metabolism','D-Glutamine and D-glutamate metabolism','Phosphotransferase system (PTS)','Vibrio cholerae pathogenic cycle','Aminobenzoate degradation','ABC transporters','Glutathione metabolism','Two-component system','beta-Alanine metabolism','Lipoic acid metabolism','Methane metabolism','Pathways in cancer','Benzoate degradation','Bacterial chemotaxis','D-Alanine metabolism','Alzheimers disease','Zeatin biosynthesis')
keep.pathways3 <- scan('~/drive/research/prism/data/microbiome/picrust_processing/kegg_l3_include.txt',w='c',sep='\n')
keep.pathways3 <- setdiff(keep.pathways3, outlier.pathways)
pw3 <- pw3[,colnames(pw3) %in% keep.pathways3]
keep.pathways2 <- scan('~/drive/research/prism/data/microbiome/picrust_processing/kegg_l2_include.txt',w='c',sep='\n')
pw2 <- pw2[,colnames(pw2) %in% keep.pathways2]

m <- droplevels(m[intersect(rownames(m),rownames(x)),])
m <- droplevels(m[!is.na(m$Gx_Subject_ID),])
m <- droplevels(m[m$Disease %in% c('CD','UC'),])
m <- droplevels(m[!is.na(m$Inflamed),])
x <- x[rownames(m),]
pw2 <- pw2[rownames(m),]
pw3 <- pw3[rownames(m),]

colnames(x) <- shorten.taxonomy(colnames(x),must.include.level=6)

cat('Loading otus...\n')
otus <- load.qiime.otu.table(otu.fp,include.lineages=TRUE)
lineages <- otus$lineages
otus <- otus$otus
rownames(otus) <- gsub('-','.',rownames(otus))
otus <- otus[rownames(m),]
otus <- sweep(otus,1,rowSums(otus),'/')
names(lineages) <- colnames(otus)


cat('Loading genetics...\n')
maf <- .1
gx <- read.table('input.v5/genotype.txt',sep='\t',head=T,row=1,check=F)
m <- droplevels(m[m$Gx_Subject_ID %in% rownames(gx),])
gx <- gx[match(m$Gx_Subject_ID,rownames(gx)),]
rownames(gx) <- rownames(m)
gx <- gx[,apply(gx,2,function(xx) mean(xx[!is.na(xx)]))/2 > maf]
otus <- otus[rownames(m),]
x <- x[rownames(m),]
pw2 <- pw2[rownames(m),]
pw3 <- pw3[rownames(m),]

class(pw2) <- 'numeric'
class(pw3) <- 'numeric'
pw2 <- pw2[,colMeans(pw2 > 0) >= min.prevalence]
pw3 <- pw3[,colMeans(pw3 > 0) >= min.prevalence]
x <- x[,colMeans(x > 0) >= min.prevalence]
lineages <- lineages[colSums(otus>0)>1]
otus <- otus[,colSums(otus>0)>1]

bdiv.uuf <- load.qiime.distance.matrix('input.v4/bdiv/unweighted_unifrac_otutable-min500.txt')
bdiv.uuf <- bdiv.uuf[rownames(m),rownames(m)]
bdiv.wuf <- load.qiime.distance.matrix('input.v4/bdiv/weighted_unifrac_otutable-min500.txt')
bdiv.wuf <- bdiv.wuf[rownames(m),rownames(m)]



LOAD.NULL <- FALSE
if(LOAD.NULL){
	cat('Loading null genetics...\n')
	gx.null <- read.table('~/drive//research/prism/data/genetic/merged/merge-prism-msh-neth-mli-filtered-pruned-null.raw',sep='\t',head=T,check=F)
	rownames(gx.null) <- gx.null[,2]
	rownames(gx.null) <- gsub('-','__',rownames(gx.null))
	gx.null <- gx.null[,-(1:6)]

	# drop samples not contained in null gx
	m <- droplevels(m[m$Gx_Subject_ID %in% rownames(gx.null),])
	gx.null <- gx.null[match(m$Gx_Subject_ID,rownames(gx.null)),]
	rownames(gx.null) <- rownames(m)
	otus <- otus[rownames(m),]
	x <- x[rownames(m),]
	gx <- gx[rownames(m),]
	gx.null <- gx.null[,apply(gx.null,2,function(xx) mean(xx[!is.na(xx)]))/2 > maf]
}

m$Immuno[m$Immuno == 'unknown'] <- 'No'
m$Immuno[is.na(m$Immuno)] <- 'No'
m$Immuno <- droplevels(m$Immuno)
m$Antibiotics[is.na(m$Antibiotics)] <- 'No'
m$Antibiotics <- droplevels(m$Antibiotics)

clin1 <- c('Age','Antibiotics','Biopsy_Location_General','Collection_Center','Disease','is.ileo','Gender','Immuno','Inflamed','Years_since_diagnosis')
clin2 <- c('Age','Antibiotics','Biopsy_Location_General','Collection_Center','Disease','is.ileo','Gender','Immuno','Inflamed')

snps <- read.table('~/drive/research/prism/data/genetic/2013-01-09-generate-snps/all-snps-metadata-with-pos-with-categories-with-NOD2.txt',sep='\t',head=T,row=1,check=F)

sample.site.order=c('Colon','Ileum','Rectum','Pouch','Pre_Pouch_Ileum')
inflamed.order=c('Yes','No')
sample.site.order=c('Ileum','Colon','Pre_Pouch_Ileum','Rectum','Pouch')
# inflamed.order=c('No','Yes')

clin <- clin2

subset.ix <- pick.unique.samples(m, sample.site.order=sample.site.order, inflamed.order=inflamed.order)

# subset.ix <- subset.ix & m$Age <=75 & m$Age >= 18 & m$Disease %in% c('CD','UC')
subset.ix <- subset.ix & m$Age <=75 & m$Age >= 1 & m$Disease %in% c('CD','UC')

subset1 <- subset.ix 

subset2 <- subset.ix & m$Biopsy_Location_General %in% c('Colon','Ileum','Pouch','Pre_Pouch_Ileum')

subset3 <- subset.ix & m$Biopsy_Location_General %in% c('Colon','Ileum')

subset4 <- subset.ix & m$Biopsy_Location_General %in% c('Colon','Ileum','Pre_Pouch_Ileum')

subset5 <- subset.ix & m$Biopsy_Location_General %in% c('Colon','Ileum','Rectum')

subset5.5 <- subset.ix & m$Biopsy_Location_General %in% c('Colon','Ileum','Pre_Pouch_Ileum','Rectum')


# drop abx
subset.ix <- subset.ix & m$Antibiotics == 'No'

subset6 <- subset.ix 

subset7 <- subset.ix & m$Biopsy_Location_General %in% c('Colon','Ileum','Pouch','Pre_Pouch_Ileum')

subset8 <- subset.ix & m$Biopsy_Location_General %in% c('Colon','Ileum')

subset9 <- subset.ix & m$Biopsy_Location_General %in% c('Colon','Ileum','Pre_Pouch_Ileum')

subset10 <- subset.ix & m$Biopsy_Location_General %in% c('Colon','Ileum','Rectum')

subset11 <- subset.ix & m$Biopsy_Location_General %in% c('Colon','Ileum','Pre_Pouch_Ileum','Rectum')

# 4, 7, 9
subsetx <- subset9

ix <- subsetx
gxx <- gx[,1,drop=F]
xx <- x
