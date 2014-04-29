# impute missing data
# Rscript impute.missing.metadata.r mappingfile taxonfile outfile
basedir <- Sys.getenv('MWAS_GWAS_DIR')
source(sprintf('%s/src/lib/rpackage/util.load.r',basedir))

covariate.names <- c('Antibiotics','Immuno','Inflamed','Age','Gender',
                     'Biopsy_Location_General','Disease','Disease_Location',
                     'Years_since_diagnosis','Collection_Center') 
args <- commandArgs(trail=TRUE)
md.fp <- args[1]
taxa.fp <- args[2]
out.fp <- args[3]

md <- load.metadata(md.fp)$x
taxa <- load.microbiome.data(taxa.fp,type='taxa')$x

# extract relevant covariates
# estimate years since diagnosis from "Immuno_close" date
cat('WARNING: several years since diagnosis undefined\n')
# set NA inflamed/Immuno to No, this is mostly HC
md$Inflamed[is.na(md$Inflamed)] <- 'No'
md$Immuno[is.na(md$Immuno)] <- 'No'

ysd <- md$Years_since_diagnosis
age <- md$Age
na.ix <- is.na(ysd)
if(sum(na.ix) > 0){
  mlm <- lm(ysd ~ age,subset=!na.ix)
  md$Years_since_diagnosis[na.ix] <- predict(mlm,data.frame(age=age[na.ix]))
  md <- droplevels(md)
}  

# necessary hacking of covariates
#     covariates$Antibiotics[is.na(covariates$Antibiotics)] <- 'unknown'
#     covariates$Antibiotics[covariates$Antibiotics == 'unknown'] <- 'unknown'
# predict missing antibiotics data
library('randomForest')
# only predict those rows for which we have microbiome data
common <- intersect(rownames(md),rownames(taxa))
orig.md <- md
md <- md[common,]
taxa <- taxa[common,]

# Impute Antibiotics
train.ix <- (md$is.msh | md$is.mgh) & !is.na(md$Antibiotics) & md$Antibiotics != 'unknown' & md$Disease != 'HC'
predict.ix <- is.na(md$Antibiotics) | md$Antibiotics == 'unknown'


if(sum(predict.ix) > 0){
  x <- as.data.frame(taxa)
  x <- cbind(x,md[,covariate.names])
  x <- x[,-match('Antibiotics',colnames(x))]
  x <- x[,-match('Disease_Location',colnames(x))]
  abx <- droplevels(factor(md$Antibiotics[train.ix]))
  for(i in 1:ncol(x)) if(class(x[,i]) == 'character') x[,i] <- droplevels(factor(x[,i]))
  mrf <- randomForest(x[train.ix,], abx)
  md$Antibiotics[predict.ix] <- predict(mrf,x[predict.ix,])
  md$Antibiotics <- droplevels(md$Antibiotics)
}

# single subject with unknown Immuno, set to "No"
md$Immuno[md$Immuno == 'unknown'] <- 'No'

# Assume most common for missing Disease Location
md$Disease_Location[md$has.cd & is.na(md$Disease_Location)] <- 'L3'
md$Disease_Location[md$has.uc & is.na(md$Disease_Location)] <- 'E3'

orig.md[rownames(md),'Disease_Location'] <- md$Disease_Location
orig.md[rownames(md),'Antibiotics'] <- md$Antibiotics
orig.md$is.inflamed <- orig.md$Inflamed == 'Yes'

sink(out.fp)
cat('#SampleID\t')
write.table(orig.md,sep='\t',quote=F)
sink(NULL)
