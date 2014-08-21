
clin1 <- c('Age','Antibiotics','Biopsy_Location_General','Collection_Center','Disease','is.ileo','Gender','Immuno','Inflamed','Years_since_diagnosis')
clin2 <- c('Age','Antibiotics','Biopsy_Location_General','Collection_Center','Disease','is.ileo','Gender','Immuno','Inflamed')
clin <- clin1

snps <- read.table('~/drive/research/prism/data/genetic/2013-01-09-generate-snps/all-snps-metadata-with-pos-with-categories-with-NOD2.txt',sep='\t',head=T,row=1,check=F)

sample.site.order=c('Colon','Ileum','Rectum','Pouch','Pre_Pouch_Ileum')
sample.site.order=c('Ileum','Colon','Pre_Pouch_Ileum','Rectum','Pouch')
inflamed.order=c('Yes','No')
# inflamed.order=c('No','Yes')

subset.ix <- pick.unique.samples(m, sample.site.order=sample.site.order, inflamed.order=inflamed.order)

subset.ix <- subset.ix & m$Age <=75 & m$Age >= 18 & m$Disease %in% c('CD','UC')
# subset.ix <- subset.ix & m$Age <=75 & m$Age >= 1 & m$Disease %in% c('CD','UC')

subset.ix <- subset.ix & m$Collection_Center != 'RISK'

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
