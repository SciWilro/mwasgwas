srcdir=Sys.getenv('MWAS_GWAS_DIR')
source(sprintf('%s/src/lib/rpackage/linear-tests.r',srcdir))
source(sprintf('%s/src/lib/rpackage/util.r',srcdir))
source(sprintf('%s/src/lib/rpackage/util.load.r',srcdir))
source(sprintf('%s/src/lib/rpackage/general-utilities.r',srcdir))
source(sprintf('%s/src/lib/rpackage/plink.r',srcdir))
#source('~/drive/enterotypes/scripts/package/clustering.r',srcdir))
library('vegan')

"load.MG.data" <- function(
  load.mbf=TRUE,              # load microbiome function data
  normalize=FALSE,            # pre-normalize within study (bad idea)
  mbf.merge.levels=FALSE,     # merge all levels of microbiome function data (bad idea)
  min.prevalence.mb=0.01,    # mb features must be present in at least this % samples
  min.prevalence.otus=0.01,   # mb features must be present in at least this % samples
  min.prevalence.mbf=0.01,   # mb features must be present in at least this % samples
  min.prevalence.gx=0.01,    # gx features must be present in at least this % samples
  min.mean.abundance.mb=0.01,    # mb features must be present in at least this % samples
  min.mean.abundance.otus=0.001,    # mb features must be present in at least this % samples
  min.mean.abundance.mbf=0.001,    # mb features must be present in at least this % samples
  min.mean.abundance.ko=0.0001,    # mb features must be present in at least this % samples
  collapse.at=0.95,            # collapse features at this level of spearman's correlation
  top.n.by.coefvar=NULL,      # keep only the top n features by coef of var
  load.ichip=TRUE,
  load.ichip.null=TRUE,
  load.pathways=TRUE,
  exclude.risk=TRUE,
  exclude.mli=FALSE,
  prefer.inflamed=FALSE,
  rarefaction.depth=c(1000,1500,2000,2500)[1],
  #             allowed.biopsy.sites=c('Colon','Ileum','Pre_Pouch_Ileum'), # keep only samples taken at these sites (NULL: keep all)
  allowed.biopsy.sites=NULL,
  covariate.names=c('Antibiotics','Immuno','Inflamed','Age','Gender',
                    'Biopsy_Location_General','Disease','Disease_Location',
                    'Years_since_diagnosis','Collection_Center')   # names of default covariates
  #             covariate.names=c('Biopsy_Location_General', 'Disease')   # names of default covariates
  ,disease.specific.snps=FALSE # use only CD-associated snps when calculating the CD ratio
  ,pathway.type=c('Reactome','Biocarta','KEGG')[1]
){
  # load data
  p <- read.table('~/drive/research/prism/data/genetic/prism/preprocessed.raw',head=T,sep='\t')
  #     p2 <- read.table('~/drive/research/prism/data/genetic/prism/2013-09-09-ashwin/preprocessed.raw',head=T,sep='\t')
  m <- read.table('~/drive/research/prism/data/genetic/msh/preprocessed.raw',head=T,sep='\t')
  mli <- read.table('~/drive/research/prism/data/genetic/mli/raw/preprocessed.raw',head=T,sep='\t')
  r <- read.table('~/drive/research/prism/data/genetic/risk/preprocessed.raw',head=T,sep='\t',comment='')
  n <- read.table('~/drive/research/prism/data/genetic/nl/preprocessed.raw',head=T,sep='\t',comment='')
  
  inflamed.suffix <- ''
  if(prefer.inflamed) inflamed.suffix <- '-prefer-inflamed'
  # load full gx
  if(load.ichip){
    gxx <- fast.read.table('~/drive/research/prism/data/genetic/merged/merge-prism-msh-neth-mli-filtered-pruned.raw',header=T,includes.rownames=T,type='character',verbose=T)
    rownames(gxx) <- gsub('-','__',gxx[,1])
    gxx <- gxx[,-(1:5)]
    class(gxx) <- 'numeric'
  } else {
    gxx <- NULL
  }
  # load just null portion of ichip
  if(load.ichip.null){
    gxx.null <- fast.read.table('~/drive/research/prism/data/genetic/merged/merge-prism-msh-neth-mli-filtered-pruned-null.raw',header=T,includes.rownames=T,type='character',verbose=T)
    rownames(gxx.null) <- gsub('-','__',gxx.null[,1])
    gxx.null <- gxx.null[,-(1:5)]
    class(gxx.null) <- 'numeric'
  } else {
    gxx.null <- NULL
  }
  
  cat('WARNING: manually deleting duplicate snp rs7608910 from prism data.\n')
  p <- p[,-grep("\\.1", colnames(p))]
  
  mb <- load.microbiome.data(sprintf('~/drive/research/prism/data/microbiome/gg11.2/taxa-rare-%d/merged.txt', rarefaction.depth))$x
  otus <- load.microbiome.data(sprintf('~/drive/research/prism/data/microbiome/gg11.2/otutable-merged-taxa-rare-%d.txt',rarefaction.depth))$x

  if(load.mbf){
    mbf <- list()
    mbf[['module']]   <- load.qiime.otu.table(sprintf('~/drive/research/prism/data/microbiome/gg12.10/otutable-merged-rare-%d-normalized-modules.txt',rarefaction.depth))
    module.desc <- read.table('~/drive/usr/ref/kegg/data/meta/kegg-modules-descriptions.txt',sep='\t',row=1)
    colnames(mbf[['module']]) <- sprintf('%s__%s',colnames(mbf[['module']]),module.desc[colnames(mbf[['module']]),1])
    mbf[['KO']] <- load.qiime.otu.table(sprintf('~/drive/research/prism/data/microbiome/gg12.10/otutable-merged-rare-%d-normalized-ko.txt',rarefaction.depth))
    ko.desc <- read.table('~/drive/usr/ref/kegg/data/pathways/ko-all-xl-nsc.txt',sep='\t',row=1,quote='"')
    colnames(mbf[['KO']])[colnames(mbf[['KO']]) %in% rownames(ko.desc)] <- sprintf('%s__%s',colnames(mbf[['KO']])[colnames(mbf[['KO']]) %in% rownames(ko.desc)],ko.desc[colnames(mbf[['KO']])[colnames(mbf[['KO']]) %in% rownames(ko.desc)],1])
    mbf[['pathway2']] <- load.qiime.otu.table(sprintf('~/drive/research/prism/data/microbiome/gg12.10/otutable-merged-rare-%d-normalized-pathways-L2.txt',rarefaction.depth))
    mbf[['pathway3']] <- load.qiime.otu.table(sprintf('~/drive/research/prism/data/microbiome/gg12.10/otutable-merged-rare-%d-normalized-pathways-L3.txt',rarefaction.depth))
    # sort by column names
    for(i in 1:length(mbf)) mbf[[i]] <- mbf[[i]][,sort(colnames(mbf[[i]]))]
  }
  # normalize to 1
  for(i in 1:length(mbf)) mbf[[i]] <- sweep(mbf[[i]], 1, rowSums(mbf[[i]]), '/')
  otus <- sweep(otus, 1, rowSums(otus),'/')
  mb <- sweep(mb, 1, rowSums(mb),'/')
  
  # drop "None" from functions
  for(i in 1:length(mbf)) mbf[[i]] <- mbf[[i]][,!grepl('None|Unclassified',colnames(mbf[[i]]))]
  
  snps <- read.table('~/drive/research/prism/data/genetic/2013-01-09-generate-snps/all-snps-metadata-with-pos.txt', sep='\t',head=T,comment='',check=F,row=1, stringsAsFactors=FALSE)
  snps.full <- snps
  map <- load.qiime.mapping.file(sprintf('~/drive/research/prism/data/metadata/map-merged%s.txt',inflamed.suffix))
  map$Microbiome_Sample_ID <- rownames(map)
  
  dash.ix <- grep('-',rownames(mbf[[1]]))
  nodash <- gsub('-','.',rownames(mbf[[1]]))[dash.ix]
  rownames(mb)[rownames(mb) %in% nodash] <- gsub('\\.','-', rownames(mb)[rownames(mb) %in% nodash])
  dash.ix <- grep('-',rownames(otus))
  nodash <- gsub('-','.',rownames(otus))[dash.ix]
  rownames(otus)[rownames(otus) %in% nodash] <- gsub('\\.','-', rownames(otus)[rownames(otus) %in% nodash])
  for(i in 1:length(mbf)) {
    dash.ix <- grep('-',rownames(mbf[[i]]))
    nodash <- gsub('-','.',rownames(mbf[[i]]))[dash.ix]
    rownames(mbf[[i]])[rownames(mbf[[i]]) %in% nodash] <- gsub('\\.','-', rownames(mbf[[i]])[rownames(mbf[[i]]) %in% nodash])
  }
  
  # keep only mb samples in map
  mb <- mb[rownames(mb) %in% map$Microbiome_Sample_ID,]
  otus <- otus[rownames(otus) %in% map$Microbiome_Sample_ID,]
  for(i in 1:length(mbf)) mbf[[i]] <- mbf[[i]][rownames(mbf[[i]]) %in% map$Microbiome_Sample_ID,]
  
  # subset map to match microbiome
  map <- droplevels(map[map$Microbiome_Sample_ID %in% rownames(mb),])
  
  # replace microbiome sample IDs with subject IDs
  rownames(mb) <- map$Gx_Subject_ID[match(rownames(mb),map$Microbiome_Sample_ID)]
  rownames(otus) <- map$Gx_Subject_ID[match(rownames(otus),map$Microbiome_Sample_ID)]
  for(i in 1:length(mbf)) rownames(mbf[[i]]) <- map$Gx_Subject_ID[match(rownames(mbf[[i]]),map$Microbiome_Sample_ID)]
  
  # replace map row names with Gx_Subject_ID
  rownames(map) <- map$Gx_Subject_ID
  
  map$Years_since_diagnosis <- map$Age - map$Age_at_IBD_diagnosis
  map$Recent_Diagnosis <- 'No'
  map$Recent_Diagnosis[map$Years_since_diagnosis <= 1/12] <- 'Yes'
  map$Biopsy_Location_General <- as.factor(map$Biopsy_Location_General)
  # merge genomics into common table
  #     common.cols <- intersect(intersect(colnames(m), colnames(p)), colnames(r))
  #     gx <- rbind(p[,common.cols],m[,common.cols],r[,common.cols])
  #     rownames(gx) <- gx[,2]
  #     gx <- as.matrix(gx[,-(1:6)])
  rownames(p) <- p[,2]; p <- p[,-(1:6)]
  #     rownames(p2) <- sapply(strsplit(p2[,1],'_'),'[',2); p2 <- p2[,-(1:6)]
  rownames(m) <- m[,2]; m <- m[,-(1:6)]
  rownames(mli) <- sprintf('%07d',as.numeric(mli[,1])); mli <- mli[,-(1:6)]
  # hac: remove stray duplicate snp from mli
  
  mli[,'rs7608910_A'] <- mli[,'rs7608910_A.1']
  mli[,'rs7608910_A.1'] <- NULL
  rownames(r) <- r[,2]; r <- r[,-(1:6)]
  rownames(n) <- n[,2]; n <- n[,-(1:6)]
  #     all.cols <- unique(c(colnames(m), colnames(mli), colnames(p), colnames(p2), colnames(r), colnames(n)))
  #     gx <- matrix(NA,nrow(p)+nrow(p2)+nrow(m)+nrow(mli)+nrow(r)+nrow(n), length(all.cols))
  #     rownames(gx) <- c(rownames(p), rownames(p2), rownames(m), rownames(mli), rownames(r), rownames(n))
  all.cols <- unique(c(colnames(m), colnames(mli), colnames(p), colnames(r), colnames(n)))
  gx <- matrix(NA,nrow(p)+nrow(m)+nrow(mli)+nrow(r)+nrow(n), length(all.cols))
  rownames(gx) <- c(rownames(p), rownames(m), rownames(mli), rownames(r), rownames(n))
  colnames(gx) <- sort(all.cols)
  gx[rownames(p),colnames(p)] <- as.matrix(p)
  #     gx[rownames(p2),colnames(p2)] <- as.matrix(p2)
  gx[rownames(m),colnames(m)] <- as.matrix(m)
  gx[rownames(mli),colnames(mli)] <- as.matrix(mli)
  gx[rownames(r),colnames(r)] <- as.matrix(r)
  gx[rownames(n),colnames(n)] <- as.matrix(n)
  
  # drop the allele from plink columns
  colnames(gx) <- sapply(strsplit(colnames(gx),'_'),function(xx) xx[1])
  snps <- snps[colnames(gx),]
  snps$AssociationType <- as.character(snps$AssociationType)
  snps$AssociationType[is.na(snps$AssociationType)] <- 'Unknown'
  
  # get common rows of gx and mb and map
  common.rows <- intersect(rownames(mb), rownames(gx))
  common.rows <- intersect(common.rows, rownames(otus))
  common.rows <- intersect(common.rows, rownames(map))
  if(load.ichip) common.rows <- intersect(common.rows, rownames(gxx))
  if(load.ichip.null) common.rows <- intersect(common.rows, rownames(gxx.null))
  if(load.mbf) common.rows <- intersect(common.rows, rownames(mbf[[1]]))
  mb <- mb[common.rows,]
  otus <- otus[common.rows,]
  if(load.mbf) for(i in 1:length(mbf)) mbf[[i]] <- mbf[[i]][common.rows,]
  if(load.ichip) gxx <- gxx[common.rows,]
  if(load.ichip.null) gxx.null <- gxx.null[common.rows,]
  gx <- gx[common.rows,]
  map <- map[common.rows,]
  # set certain metavariables as factors
  metavars <- c('Gender','Antibiotics','Study')
  for(i in 1:length(metavars)){
    map[,metavars[i]] <- as.factor(map[,metavars[i]])
  }
  map$Gender <- as.factor(map$Gender)
  location <- factor(map$Montreal_Classification_disease_location,levels=c('L1','L2','L3','E1','E2','E3','Lx'))
  location[grep('L1',map$Montreal_Classification_disease_location)] <- 'L1'
  location[grep('L2',map$Montreal_Classification_disease_location)] <- 'L2'
  location[grep('L3',map$Montreal_Classification_disease_location)] <- 'L3'
  location[map$Disease == 'CD'] <- gsub('E','L',location[map$Disease == 'CD'])
  location[map$Disease == 'UC'] <- gsub('L','E',location[map$Disease == 'UC'])
  location[map$Disease %in% c('IC','HC')] <- 'Lx'
  location <- as.character(location)
  map$Ileal <- grepl('1|3',location)
  map$Ileal[location == 'Lx'] <- NA
  map$Disease_Location <- location
  
  map$Immuno[map$Immuno == 'N'] <- 'No'
  map$Immuno[map$Immuno == 'Y'] <- 'Yes'
  
  # missing.mb <- setdiff(rownames(mb),common.rows)
  # missing.gx <- setdiff(rownames(gx),common.rows)
  # missing.snps <- setdiff(rownames(snps), sapply(strsplit(common,'_'),'[[',1))
  # MANDATORY SUBSETS
  # keep only subjects with abx info (include all now--risk has NA for abx)
  #     keep.subj <- !is.na(map$Antibiotics) & !is.na(map$Age)
  
  
  keep.subj <- rep(TRUE,nrow(map))	
  #     keep.subj <- keep.subj & !is.na(map$Age)
  #     keep.subj <- keep.subj & map$Disease != 'HC'
  #     keep.subj <- keep.subj & map$Disease != 'IC'
  #     keep.subj <- keep.subj & (map$Disease == 'CD' | map$Disease == 'UC')
  if(exclude.risk) keep.subj <- keep.subj & map$Study != 'RISK'
  if(exclude.mli) keep.subj <- keep.subj & map$Study != 'MLI'
  # keep biopsies only at Colon, Ileum?
  if(!is.null(allowed.biopsy.sites)){
    keep.subj <- keep.subj & map$Biopsy_Location_General %in% allowed.biopsy.sites
  }
  
  # drop NA's at Immuno, Disease_Location, Antibiotics, Years_since_IBD_diagnosis
  # 	keep.subj <- keep.subj & !is.na(map$Antibiotics) & !(map$Antibiotics == 'unknown')
  # 	keep.subj <- keep.subj & !is.na(map$Immuno) & !(map$Immuno == 'unknown')
  # 	keep.subj <- keep.subj & !is.na(map$Disease_Location) & !(map$Disease_Location == 'unknown')
  # 	keep.subj <- keep.subj & !is.na(map$Years_since_diagnosis) & !(map$Years_since_diagnosis == 'unknown')
  
  map <- map[keep.subj,]
  gx <- gx[keep.subj,]
  mb <- mb[keep.subj,]
  otus <- otus[keep.subj,]
  if(load.mbf) for(i in 1:length(mbf)) mbf[[i]] <- mbf[[i]][keep.subj,]
  if(load.ichip) gxx <- gxx[keep.subj,]
  if(load.ichip.null) gxx.null <- gxx.null[keep.subj,]
  
  # drop taxa in > min.prevalence samples
  keep.taxa <- colMeans(mb>0) >= min.prevalence.mb & colMeans(mb) >= min.mean.abundance.mb
  mb <- mb[, keep.taxa]
  keep.otus <- colMeans(otus>0) >= min.prevalence.otus & colMeans(otus) >= min.mean.abundance.otus
  otus <- otus[, keep.otus]
  if(load.mbf) {
    for(i in 1:length(mbf)){
      if(names(mbf)[i] == 'KO'){
        keep.mbf <- colMeans(mbf[[i]]>0) >= min.prevalence.mbf & colMeans(mbf[[i]]) >= min.mean.abundance.ko
      } else {
        keep.mbf <- colMeans(mbf[[i]]>0) >= min.prevalence.mbf & colMeans(mbf[[i]]) >= min.mean.abundance.mbf
      } 
      mbf[[i]] <- mbf[[i]][, keep.mbf]
      if(!is.null(top.n.by.coefvar)){
        if(ncol(mbf[[i]]) > top.n.by.coefvar){
          coefvar <- apply(mbf[[i]],2,function(xx) sd(xx)/mean(xx))
          mbf[[i]] <- mbf[[i]][,order(coefvar,decreasing=TRUE)[1:top.n.by.coefvar]]
        }
      }
    }
  }
  # keep only genera in > .1 samples
  # genus.ix <- which(grepl('g__',colnames(mb)) & !grepl('s__',colnames(mb)))
  # mb <- mb[,genus.ix]
  # keep only snps in jostins et al
  #     snp.ix <- snps$Reference == '23128233'
  #     snps <- snps[snp.ix,]
  
  gx <- gx[,rownames(snps)]
  good.snps <- grepl('Daly',snps$Reference) | grepl('23128233',snps$Reference)
  # keep only snps called in 99% of subjects
  #     good.snps <- apply(is.na(gx),2,mean) < min.prevalence.gx
  #     good.snps <- colMeans(gx>0)
  snps <- snps[good.snps,]
  gx <- gx[,good.snps]
  
  mb.collapse <- NULL
  otus.collapse <- NULL
  mbf.collapse <- NULL
  # collapse once at top level
  if(!is.null(collapse.at)){
    mb.collapse <- collapse.by.correlation(mb, collapse.at, "mean", verbose=TRUE)
    mb <- mb[,mb.collapse$reps]
    otus.collapse <- collapse.by.correlation(otus, collapse.at, "mean", verbose=TRUE)
    otus <- otus[,otus.collapse$reps]
    
    if(load.mbf) {
      mbf.collapse <- list()
      for(i in 1:length(mbf)){
        mbf.collapse[[i]] <- collapse.by.correlation(mbf[[i]], collapse.at, "mean", verbose=TRUE)
        mbf[[i]] <- mbf[[i]][,mbf.collapse[[i]]$reps]
      }
      names(mbf.collapse) <- names(mbf)
    }
    # collapse gx
    #         gx.collapse <- collapse.by.correlation(gx, collapse.at, "mean", verbose=TRUE)
    #         gx <- gx[,gx.collapse$reps]
    #         snps <- snps[gx.collapse$reps,]
  }
  # OPTIONAL SUBSETS
  # SUBJECT SUBSETS
  has.cd <- map$Disease == 'CD'
  has.uc <- map$Disease == 'UC'
  is.prism <- map$Study == 'MGH' | map$Study == 'MSH'
  is.msh <- map$Collection_Center == 'MSH'
  is.mgh <- map$Collection_Center == 'MGH'
  is.mli <- map$Collection_Center == 'MLI'
  is.risk <- map$Study == 'RISK'
  is.neth <- map$Study == 'NETH'
  is.colon <- map$Biopsy_Location_General == 'Colon'
  is.ileum <- map$Biopsy_Location_General == 'Ileum'
  is.pouch <- map$Biopsy_Location_General == 'Pouch'
  is.prepouch <- map$Biopsy_Location_General == 'Pre_Pouch_Ileum'
  is.rectum <- map$Biopsy_Location_General == 'Rectum'
  is.ileo <- grepl('L1',map$Montreal_Classification_disease_location)
  if(!exclude.risk){
    map$Immuno[is.risk] <- 'No'
    map$Antibiotics[is.risk] <- 'No'
  }
  map$Recent_Diagnosis[is.prepouch] <- 'No'
  
  # GENETIC SUBSETS
  # subset of snps on important genes
  jostins.snps <- grepl('23128233',snps$Reference)
  mb.snps <- grepl('ATG16L1|IL23R|NOD2|IRGM',snps$`All Genes`) & grepl('23128233',snps$Reference)
  mb.snps <- grepl('ATG16L1|IL23R|NOD2|IRGM|CARD9|FUT2',snps$`All Genes`) & grepl('23128233',snps$Reference)
  #     atg.snps <- grepl('ATG16L1|IRGM',snps$`All Genes`) & grepl('23128233',snps$Reference)
  atg.snps <- grepl('ATG16L1|NOD2|IRGM',snps$`All Genes`) & grepl('23128233',snps$Reference)
  atg.snps2 <- grepl('ATG16L1|NOD2|IRGM|ORMDL3',snps$`All Genes`) & grepl('23128233',snps$Reference)
  fut.snps <- grepl('FUT2',snps$`All Genes`) & grepl('23128233',snps$Reference)
  jostins.nw.snps <- grepl('NOD2|IL10|HCK|DOK3|VDR|SLC11A1|CARD9|LGALS9',snps$`All Genes`)
  jakstat.snps <- grepl('JAK|STAT',snps$`All Genes`) & grepl('23128233',snps$Reference)
  mb.snps.long <- grepl('ATG16L1|IL23R|IL12B|CARD9|MUC19|TAB1|FUT2|NOD2|IRGM',snps$`All Genes`)  & grepl('23128233',snps$Reference)
  cd.snps <- snps$AssociationType == 'CD' | snps$AssociationType == 'IBD'
  uc.snps <- snps$AssociationType == 'UC' | snps$AssociationType == 'IBD'
  ibd.snps <- snps$AssociationType == 'IBD'
  nod2.snps.long <- grepl('NOD2',snps$`All Genes`)
  #     nod2.snps <- grepl('NOD2',snps$`All Genes`) & grepl('23128233|Daly',snps$Reference)
  daly.snps <- grepl('Daly',snps$Reference)
  nod2.snps <- grepl('NOD2',snps$`All Genes`) & daly.snps
  sensing.snps <- nod2.snps | (grepl('CARD9',snps$`All Genes`) & grepl('23128233',snps$Reference))
  #     rx1.snps <- grepl('SP110',snps$`All Genes`)
  #     rx2.snps <- grepl('ATG16L1|ATG5|IRGM|ANKRD33B,DAP',snps$`All Genes`)
  #     rx3.snps <- grepl('FUT2|PTGER4|MUC19',snps$`All Genes`)
  #     rx4.snps <- grepl('IL23R|IL12B|STAT3|JAK2|CCR6',snps$`All Genes`)
  #     rx5.snps <- grepl('C1orf106',snps$`All Genes`)
  #     rx6.snps <- grepl('IL18RAP',snps$`All Genes`)
  #     rx7.snps <- grepl('MST1',snps$`All Genes`)
  #     rx8.snps <- grepl('CARD9',snps$`All Genes`)
  
  # drop very rare NOD2 snps
  daly.snps[daly.snps] <- apply(gx[,daly.snps],2,function(xx) mean(drop.na(xx)/2) > 0.01)
  
  topgenes <- as.character(snps$`All Genes`)
  names(topgenes) <- rownames(snps)
  topgenes[topgenes == ''] <- rownames(snps)[topgenes=='']
  topgenes['rs7134599'] <- 'IL22'
  topgenes['rs12942547'] <- 'IL22'
  topgenes['rs12942547'] <- 'STAT3'
  topgenes['rs10758669'] <- 'JAK2'
  topgenes['rs4246905'] <- 'TNFSF15'
  topgenes['rs1569328'] <- 'FOS'
  topgenes['rs1893217'] <- 'PTPN2'
  topgenes['rs4845604'] <- 'RORC'
  topgenes['rs11168249'] <- 'VDR'
  topgenes['rs3024505'] <- 'IL10'
  topgenes['rs4243971'] <- 'HCK'
  topgenes['rs4976646'] <- 'DOK3'
  topgenes['rs2382817'] <- 'SLC11A1'
  topgenes['rs2945412'] <- 'LGALS9'	
  topgenes[grep('SP110',topgenes)] <- 'SP110'
  topgenes[grep('ANKRD33B,DAP',topgenes)] <- 'DAP'
  topgenes[grep('ANKRD33B,DAP',topgenes)] <- 'DAP'
  topgenes[grep('ATG5',topgenes)] <- 'ATG5'
  topgenes[grep('PTGER4',topgenes)] <- 'PTGER4'
  topgenes[grep('MUC19',topgenes)] <- 'MUC19'
  topgenes[grep('IL12B',topgenes)] <- 'IL12B'
  topgenes[grep('CCR6',topgenes)] <- 'CCR6'
  topgenes[grep('IL18RAP',topgenes)] <- 'IL18RAP'
  topgenes[grep('C1orf106',topgenes)] <- 'C1orf106'
  topgenes[grep('MST1',topgenes)] <- 'MST1'
  topgenes[grep('TAB1',topgenes)] <- 'TAB1'
  topgenes[grep('CARD9',topgenes)] <- 'CARD9'
  topgenes[grep('IL23R',topgenes)] <- 'IL23R'
  topgenes[grep('IRGM',topgenes)] <- 'IRGM'
  topgenes[grep('ATG16L1',topgenes)] <- 'ATG16L1'
  topgenes[grep('FUT2',topgenes)] <- 'FUT2'
  topgenes['rs7134599'] <- 'IFNG'
  topgenes['rs913678'] <- 'CEBPB'
  topgenes['rs727088'] <- 'CD226'
  topgenes['rs921720'] <- 'TRIB1'
  topgenes['rs6740462'] <- 'SPRED2'
  topgenes['rs11739663'] <- 'SLC9A3'
  topgenes['rs2227551'] <- 'rs2227551'
  topgenes['rs10865331'] <- 'rs10865331'
  topgenes['rs11168249'] <- 'HDAC7/VDR'
  topgenes['rs2266959'] <- 'rs2266959'
  topgenes['rs2231884'] <- 'RELA'
  topgenes['rs9847710'] <- 'PRKCD'
  topgenes['rs670523'] <- 'rs670523'
  topgenes['rs4743820'] <- 'NFIL3'
  topgenes['rs6716753'] <- 'SP140'
  topgenes['rs1728785'] <- 'ZFP90'
  topgenes['rs1847472'] <- 'rs1847472'
  topgenes[''] <- ''
  topgenes[''] <- ''
  topgenes[''] <- ''
  topgenes[''] <- ''
  topgenes[''] <- ''
  topgenes[''] <- ''
  topgenes[''] <- ''
  topgenes[''] <- ''
  topgenes[''] <- ''
  topgenes[''] <- ''
  topgenes[''] <- ''
  
  topgenes <- c(topgenes,'NOD2')
  names(topgenes)[length(topgenes)] <- 'NOD2'
  topgenes <- sapply(strsplit(topgenes,','),'[',1)
  
  # MICROBIOME SUBSETS
  known.taxa <- scan('~/drive/research/prism/data/metadata/taxon-assoc-all.txt',what='c',quiet=T)
  known.taxa <- colnames(mb) %in% known.taxa
  cd.taxa <- colnames(mb) %in% scan('~/drive/research/prism/data/metadata/taxon-assoc-cd.txt',what='c',quiet=T)
  uc.taxa <- colnames(mb) %in% scan('~/drive/research/prism/data/metadata/taxon-assoc-uc.txt',what='c',quiet=T)
  ibd.taxa <- colnames(mb) %in% scan('~/drive/research/prism/data/metadata/taxon-assoc-ibd.txt',what='c',quiet=T)
  cd.taxa <- cd.taxa | colnames(mb) == "pBctrdts;cBctrd;oBctrdls;fPrvtllc"
  uc.taxa <- uc.taxa | colnames(mb) == "pBctrdts;cBctrd;oBctrdls;fPrvtllc"
  ibd.taxa <- ibd.taxa | colnames(mb) == "pBctrdts;cBctrd;oBctrdls;fPrvtllc"
  
  
  # gxx PC
  if(load.ichip){
    cat('Adding genetic principal coordinates...\n')
    gxx.pc <- replicate(10,{
      tmp <- gxx
      tmp[is.na(tmp)] <- 0
      princomp(tmp[,sample(ncol(tmp),size=ceiling(.99*nrow(tmp)))])$scores[,1:3]
    })
    gxx.pc <- apply(gxx.pc,1:2,mean)
  } else {
    gxx.pc <- NULL
  }
  
  # rename taxa by abbreviated name
  abbrev.taxon.names <- sapply(strsplit(colnames(mb),';'),function(xx) {xx <- xx[-1]; paste(gsub('~','o',gsub('__','',gsub('[aeiouy]','',gsub('o__','~__',xx)))),collapse=';')})
  full.taxon.names <- colnames(mb)
  colnames(mb) <- abbrev.taxon.names
  
  # create log version of mb
  lmb <- mb; lmb[lmb==0] <- min(lmb[lmb > 0]) / 2; lmb <- log10(lmb)
  #     lmb <- asin(sqrt(mb));
  if(load.mbf) {
    lmbf <- list()
    for(i in 1:length(mbf)){
      lmbf[[i]] <- mbf[[i]]
      lmbf[[i]][lmbf[[i]]==0] <- min(lmbf[[i]][lmbf[[i]] > 0]) / 2
      lmbf[[i]] <- log10(lmbf[[i]])
    }
    names(lmbf) <- names(mbf)
  }
  
  # DERIVED FEATURES
  # load alpha, beta div
  LOAD.DIV <- FALSE
  LOAD.DIV <- TRUE
  bdiv <- NULL
  dx <- NULL
  if(LOAD.DIV){
    adiv <- read.table(sprintf('~/drive/research/prism/data/microbiome/gg11.2/alpha-rare-%d.txt',rarefaction.depth),sep='\t',head=T,check=F,row=1)
    bdiv.wuf <- load.qiime.distance.matrix(sprintf('~/drive/research/prism/data/microbiome/gg11.2/beta-rare-%d/weighted_unifrac_otutable-merged-taxa-rare-%d.txt',rarefaction.depth,rarefaction.depth))
    bdiv.uuf <- load.qiime.distance.matrix(sprintf('~/drive/research/prism/data/microbiome/gg11.2/beta-rare-%d/unweighted_unifrac_otutable-merged-taxa-rare-%d.txt',rarefaction.depth,rarefaction.depth))
    
    adiv <- adiv[rownames(adiv) %in% map$Microbiome_Sample_ID,]
    rownames(adiv) <- map$Gx_Subject_ID[match(rownames(adiv),map$Microbiome_Sample_ID)]
    bdiv.wuf <- bdiv.wuf[rownames(bdiv.wuf) %in% map$Microbiome_Sample_ID,rownames(bdiv.wuf) %in% map$Microbiome_Sample_ID]
    rownames(bdiv.wuf) <- map$Gx_Subject_ID[match(rownames(bdiv.wuf),map$Microbiome_Sample_ID)]
    colnames(bdiv.wuf) <- map$Gx_Subject_ID[match(rownames(bdiv.wuf),map$Microbiome_Sample_ID)]
    bdiv.uuf <- bdiv.uuf[rownames(bdiv.uuf) %in% map$Microbiome_Sample_ID,rownames(bdiv.uuf) %in% map$Microbiome_Sample_ID]
    rownames(bdiv.uuf) <- map$Gx_Subject_ID[match(rownames(bdiv.uuf),map$Microbiome_Sample_ID)]
    colnames(bdiv.uuf) <- map$Gx_Subject_ID[match(rownames(bdiv.uuf),map$Microbiome_Sample_ID)]
    
    pc.wuf <- cmdscale(bdiv.wuf,2)
    pc.uuf <- cmdscale(bdiv.uuf,2)
    pc1.wuf <- pc.wuf[,1]
    pc2.wuf <- pc.wuf[,2]
    pc1.uuf <- pc.uuf[,1]
    pc2.uuf <- pc.uuf[,2]
    dx <- cbind(adiv,pc1.wuf,pc2.wuf,pc1.uuf,pc2.uuf)
    colnames(dx)[ncol(dx) -(3:0)] <- c('W. UniFrac PC1','W. UniFrac PC2','U. UniFrac PC1','U. UniFrac PC2')
    
    d.js <- jsdist(cbind(mbf[['pathway3']],1-rowSums(mbf[['pathway3']])))
    pc.js <- cmdscale(d.js,2)
    d.bc <- as.matrix(vegdist(mbf[['pathway3']]))
    pc.bc <- cmdscale(d.bc,2)
    dx[['Pathway-3 Bray Curtis PC1']] <- pc.bc[,1]
    dx[['Pathway-3 Bray Curtis PC2']] <- pc.bc[,2]
    dx[['Pathway-3 Jensen Shannon PC1']] <- pc.js[,1]
    dx[['Pathway-3 Jensen Shannon PC2']] <- pc.js[,2]
    
    # also get module-level PC
    d.js.m <- jsdist(cbind(mbf[['module']],1-rowSums(mbf[['module']])))
    pc.js.m <- cmdscale(d.js.m,2)
    d.bc.m <- as.matrix(vegdist(mbf[['module']]))
    pc.bc.m <- cmdscale(d.bc.m,2)
    dx[['Module Bray Curtis PC1']] <- pc.bc.m[,1]
    dx[['Module Bray Curtis PC2']] <- pc.bc.m[,2]
    dx[['Module Jensen Shannon PC1']] <- pc.js.m[,1]
    dx[['Module Jensen Shannon PC2']] <- pc.js.m[,2]
    bdiv <- list('wuf'=bdiv.wuf,'uuf'=bdiv.uuf,
                 'pathway3-js'=as.matrix(d.js),'pathway3-bc'=as.matrix(d.bc),
                 'module-js'=as.matrix(d.js.m),'module-bc'=as.matrix(d.bc.m))
  }
  
  # correct for allele inversion in several prism snps
  #     mgh.means <- apply(gx[map$Collection_Center == 'MGH',],2,function(xx) mean(drop.na(xx))/2)
  #     msh.means <- apply(gx[map$Collection_Center == 'MSH',],2,function(xx) mean(drop.na(xx))/2)
  #     if(!exclude.risk){
  #         risk.means <- apply(gx[map$Collection_Center == 'RISK',],2,function(xx) mean(drop.na(xx))/2)
  #     }
  #     disagreement <- apply(cbind(mgh.means, msh.means, risk.means),1,var)
  #     gx[which(map$Collection_Center == 'MGH'),disagreement > 0.01] <- 2 - gx[which(map$Collection_Center == 'MGH'),disagreement > 0.01]
  
  # calculate risk factors for cd, uc, ibd using microbiome-related genes
  # note: some NA's
  # for microbiome-related snps, drop subjects with NA
  or.snp.ix <- !is.na(snps$OR_IBD) | !is.na(snps$OR_CD)
  or.subj.ix <- rowSums(is.na(gx[,or.snp.ix])) == 0
  
  # compute risk scores
  diseases <- c('CD','UC','IBD')
  studies <- c('ALL')
  odds.ratio.names <- c(CD='OR_CD',UC='OR_UC',IBD='OR_IBD')
  subject.subsets <- list(ALL=list(CD=has.cd,
                                   UC=has.uc,
                                   IBD=has.cd | has.uc)
  )
  snp.subsets <- list(CD= cd.snps | ibd.snps,
                      UC= uc.snps | ibd.snps,
                      IBD=ibd.snps | cd.snps | uc.snps)
  if(!load.pathways){
    cd.scores=NULL
    uc.scores=NULL
    ibd.scores=NULL
    cd.scores.pathway=NULL
    uc.scores.pathway=NULL
    ibd.scores.pathway=NULL
    disease.scores=NULL
    pathways=NULL
    pathways.gene=NULL
    pathway.sizes=NULL
  } else{
    cat('Loading pathways...\n')
    # list of list of matrices of sample x pathway risk score
    disease.scores <- list()
    # load pathways
    snp.genes <- sapply(snps[,'All Genes'],function(xx) strsplit(xx,',')[[1]])
    names(snp.genes) <- rownames(snps)
    if(pathway.type == 'Biocarta'){
      pathways.res <- pathway.membership(snp.genes,verbose=TRUE,pathway.fp='~/drive/usr/ref/MSigDB/c2.cp.biocarta.v3.1.symbols.gmt')
    } else if(pathway.type == 'KEGG'){
      pathways.res <- pathway.membership(snp.genes,verbose=TRUE,pathway.fp='~/drive/usr/ref/MSigDB/c2.cp.kegg.v3.1.symbols.gmt')
    } else if(pathway.type == 'Reactome'){
      pathways.res <- pathway.membership(snp.genes,verbose=TRUE,pathway.fp='~/drive/usr/ref/MSigDB/c2.cp.reactome.v3.1.symbols.gmt')
    } else {
      stop('Unknown pathway type:',pathway.type,'\n')
    }
    
    pathways <- pathways.res$pathways.snp
    pathways.gene <- pathways.res$pathways.gene
    pathway.sizes <- pathways.res$pathway.sizes
    # add a few custom pathways
    # 		pathways <- rbind(mb.snps.long,pathways)
    # 		rownames(pathways)[1] <- 'MICROBIOME_DYSREGULATION_LONG'
    # 		pathways <- rbind(mb.snps, pathways)
    # 		rownames(pathways)[1] <- 'MICROBIOME_DYSREGULATION'
    pathways <- rbind(atg.snps, pathways)
    rownames(pathways)[1] <- 'AUTOPHAGY'
    # 		pathways <- rbind(atg.snps2, pathways)
    # 		rownames(pathways)[1] <- 'AUTOPHAGY2'
    # 		pathways <- rbind(sensing.snps, pathways)
    # 		rownames(pathways)[1] <- 'LIGAND_SENSING'
    pathways <- rbind(jakstat.snps, pathways)
    rownames(pathways)[1] <- 'JAK_STAT'
    pathways <- rbind(rep(TRUE,ncol(pathways)), pathways)
    rownames(pathways)[1] <- 'FULL'
    # 		pathways <- rbind(fut.snps, pathways)
    # 		rownames(pathways)[1] <- 'FUT'
    pathways <- rbind(jostins.nw.snps, pathways)
    rownames(pathways)[1] <- 'JNW'
    
    # 		pathways <- rbind(rx1.snps, pathways)
    # 		rownames(pathways)[1] <- 'RX1'
    # 		pathways <- rbind(rx2.snps, pathways)
    # 		rownames(pathways)[1] <- 'RX2'
    # 		pathways <- rbind(rx3.snps, pathways)
    # 		rownames(pathways)[1] <- 'RX3'
    # 		pathways <- rbind(rx4.snps, pathways)
    # 		rownames(pathways)[1] <- 'RX4'
    # 		pathways <- rbind(rx5.snps, pathways)
    # 		rownames(pathways)[1] <- 'RX5'
    # 		pathways <- rbind(rx6.snps, pathways)
    # 		rownames(pathways)[1] <- 'RX6'
    # 		pathways <- rbind(rx7.snps, pathways)
    # 		rownames(pathways)[1] <- 'RX7'
    # 		pathways <- rbind(rx8.snps, pathways)
    # 		rownames(pathways)[1] <- 'RX8'
    # 		pathways <- pathways[-(1:16),]
    
    #  drop singleton pathways
    single.gene.pathways <- rowSums(pathways) == 1
    pathways <- pathways[!single.gene.pathways,]
    
    #         pathways <- pathways[c(grep('NOD1|IMMUNE|TCR',rownames(pathways)), 1:8),]
    cat('calculating disease scores...\n')
    disease.scores <- list()
    for(disease in diseases){
      disease.scores[[disease]] <- matrix(NA,nrow=nrow(gx),ncol=nrow(pathways))
      for(study in studies){
        rownames(disease.scores[[disease]]) <- rownames(gx)
        colnames(disease.scores[[disease]]) <- rownames(pathways)
        for(i in 1:nrow(pathways)){
          subject.ix <- subject.subsets[[study]][[disease]]
          or.snp.ix <- !is.na(snps[,odds.ratio.names[[disease]]])
          
          snp.ix <- or.snp.ix & pathways[i,]
          if(disease.specific.snps) snp.ix <- snps.ix & snp.subsets[[disease]]
          if(sum(snp.ix) >= 3){
            
            scores <- rep(NA,nrow(gx))
            drop.ix <- try(drop.na.2d(gx[subject.ix,snp.ix,drop=F],prefer.columns=2))
            if(class(drop.ix) == 'try-error') browser()
            # drop singleton pathways
            if(sum(!drop.ix$drop.col) > 1){
              B <- snps[snp.ix,][!drop.ix$drop.col,odds.ratio.names[[disease]]]
              X <- gx[subject.ix,snp.ix,drop=F][!drop.ix$drop.row, !drop.ix$drop.col,drop=F]
              if(any(snp.ix[!drop.ix$drop.col])){
                scores <- as.numeric(B %*% t(X))
                disease.scores[[disease]][subject.ix,i][!drop.ix$drop.row] <- scores
                if(all(is.na(scores))){
                  cat('Warning: pathway,',rownames(pathways)[i],'all NA\n')
                }
              }
              #                 } else if(sum(!drop.ix$drop.col) == 1){
              #                     disease.scores[[disease]][subject.ix,i][!drop.ix$drop.row] <- scores
            } else if(sum(!drop.ix$drop.col) == 1){
              X <- gx[subject.ix,snp.ix,drop=F][!drop.ix$drop.row, !drop.ix$drop.col,drop=F]
              disease.scores[[disease]][subject.ix,i][!drop.ix$drop.row] <- X
            }
          }
        }
        drop.ix <- colMeans(is.na(disease.scores[[disease]])) == 1
        disease.scores[[disease]][,!drop.ix]
      }
    }
    cat('Done with disease scores.\n')
    cd.scores <- uc.scores <- ibd.scores <- rep(NA,nrow(gx))
    cd.scores[or.subj.ix] <- as.numeric(snps[or.snp.ix & (cd.snps | ibd.snps),'OR_CD'] %*% t(gx[or.subj.ix,or.snp.ix & (cd.snps | ibd.snps)]))
    uc.scores[or.subj.ix] <- as.numeric(snps[or.snp.ix & (uc.snps | ibd.snps),'OR_UC'] %*% t(gx[or.subj.ix,or.snp.ix & (uc.snps | ibd.snps)]))
    ibd.scores[or.subj.ix] <- as.numeric(snps[or.snp.ix,'OR_IBD'] %*% t(gx[or.subj.ix,or.snp.ix]))
    
    # load risk scores on a per-pathway basis
    snp.genes <- sapply(snps[,'All Genes'],function(xx) strsplit(xx,',')[[1]])
    snp.genes <- apply(snps[,c('GRAIL genes','DAPPLE genes','cSNP','eQTL Genes','Co-expression net. Genes')],1, function(xx) {xx <-strsplit(paste(xx,collapse=','),',')[[1]]; xx[xx != '']})
    names(snp.genes) <- rownames(snps)
    
    # load per-pathway risk scores
    cd.scores.pathway <- matrix(NA,nrow=nrow(gx),ncol=nrow(pathways))
    uc.scores.pathway <- matrix(NA,nrow=nrow(gx),ncol=nrow(pathways))
    ibd.scores.pathway <- matrix(NA,nrow=nrow(gx),ncol=nrow(pathways))
    rownames(cd.scores.pathway) <- rownames(gx)
    rownames(uc.scores.pathway) <- rownames(gx)
    rownames(ibd.scores.pathway) <- rownames(gx)
    colnames(cd.scores.pathway) <- rownames(pathways)
    colnames(uc.scores.pathway) <- rownames(pathways)
    colnames(ibd.scores.pathway) <- rownames(pathways)
    for(i in 1:nrow(pathways)){
      in.pathway <- pathways[i,]
      ix <- in.pathway & or.snp.ix & (cd.snps | ibd.snps)
      if(any(ix)) cd.scores.pathway[or.subj.ix,i] <- as.numeric(snps[ix,'OR_CD'] %*% t(gx[or.subj.ix,ix]))
      ix <- in.pathway & or.snp.ix & (uc.snps | ibd.snps)
      if(any(ix)) uc.scores.pathway[or.subj.ix,i] <- as.numeric(snps[ix,'OR_UC'] %*% t(gx[or.subj.ix,ix]))
      ix <- in.pathway & or.snp.ix
      if(any(ix)) ibd.scores.pathway[or.subj.ix,i] <- as.numeric(snps[ix,'OR_IBD'] %*% t(gx[or.subj.ix,ix]))
    }
    cd.scores.pathway <- cd.scores.pathway[,colMeans(is.na(cd.scores.pathway)) < 1]
    uc.scores.pathway <- uc.scores.pathway[,colMeans(is.na(uc.scores.pathway)) < 1]
    ibd.scores.pathway <- ibd.scores.pathway[,colMeans(is.na(ibd.scores.pathway)) < 1]
    
    # add 'NOD2' column to pathways
    pathways <- cbind(pathways,rep(FALSE,nrow(pathways)))
    colnames(pathways)[ncol(pathways)] <- 'NOD2'
    for(i in 1:nrow(pathways)){
      
      matches <- colnames(pathways)[pathways[i,]] %in% topgenes[names(topgenes == 'NOD2')]
      if(any(matches)){
        pathways[i,'NOD2'] <- TRUE
      }
    }
    
  }
  
  
  # normalize within each disease and study
  normalize <- FALSE
  if(normalize){
    disease.list <- list(cd=has.cd, uc=has.uc, ibd=has.cd | has.uc)[3]
    study.list <- list(prism=rep(TRUE,nrow(map)))
    for(k in 1:length(disease.list)){
      for(l in 1:length(study.list)){
        cat(names(disease.list)[k],'and',names(study.list)[l],'...\n')
        disease <- disease.list[[k]]
        study <- study.list[[l]]
        for(j in 1:ncol(mb)){
          lmb[disease & study,j] <- lmb[disease & study,j] - mean(lmb[disease & study,j])
          lmb[disease & study,j] <- lmb[disease & study,j] / sd(lmb[disease & study,j])
        }
        if(load.mbf){
          for(i in 1:length(lmbf)){
            for(j in 1:ncol(lmbf[[i]])){
              lmbf[[i]][disease & study,j] <- lmbf[[i]][disease & study,j] - mean(lmbf[[i]][disease & study,j])
              lmbf[[i]][disease & study,j] <- lmbf[[i]][disease & study,j] / sd(lmbf[[i]][disease & study,j])
            }
          }
        }
      }
    }
  }
  
  
  # extract relevant covariates
  # estimate years since diagnosis from "Immuno_close" date
  cat('WARNING: several years since diagnosis undefined\n')
  IMPUTE <- TRUE
  if(IMPUTE){
    # set NA inflamed/Immuno to No, this is mostly HC
    map$Inflamed[is.na(map$Inflamed)] <- 'No'
    map$Immuno[is.na(map$Immuno)] <- 'No'
    
    ysd <- map$Years_since_diagnosis
    age <- map$Age
    na.ix <- is.na(ysd)
    if(sum(na.ix) > 0){
      mlm <- lm(ysd ~ age,subset=!na.ix)
      map$Years_since_diagnosis[na.ix] <- predict(mlm,data.frame(age=age[na.ix]))
      map <- droplevels(map)
    }	
    
    # necessary hacking of covariates
    #     covariates$Antibiotics[is.na(covariates$Antibiotics)] <- 'unknown'
    #     covariates$Antibiotics[covariates$Antibiotics == 'unknown'] <- 'unknown'
    # predict missing antibiotics data
    library('randomForest')
    # Impute Antibiotics
    train.ix <- (is.msh | is.mgh) & !is.na(map$Antibiotics) & map$Antibiotics != 'unknown' & map$Disease != 'HC'
    predict.ix <- is.na(map$Antibiotics) | map$Antibiotics == 'unknown'
    
    if(sum(predict.ix) > 0){
      x <- as.data.frame(mb)
      x <- cbind(x,mbf[['module']])
      x <- cbind(x,map[,covariate.names,drop=F])
      x <- x[,-match('Antibiotics',colnames(x))]
      x <- x[,-match('Disease_Location',colnames(x))]
      abx <- droplevels(factor(map$Antibiotics[train.ix]))
      for(i in 1:ncol(x)) if(class(x[,i]) == 'character') x[,i] <- droplevels(factor(x[,i]))
      mrf <- randomForest(x[train.ix,], abx)
      map$Antibiotics[predict.ix] <- predict(mrf,x[predict.ix,])
      map$Antibiotics <- droplevels(map$Antibiotics)
    }
    
    # single subject with unknown Immuno, set to "No"
    map$Immuno[map$Immuno == 'unknown'] <- 'No'
    
    # Assume most common for missing Disease Location
    map$Disease_Location[has.cd & is.na(map$Disease_Location)] <- 'L3'
    map$Disease_Location[has.uc & is.na(map$Disease_Location)] <- 'E3'
  }
  covariates <- map[,covariate.names,drop=F]
  for(i in 1:ncol(covariates)) if(class(covariates[,i]) == 'character') covariates[,i] <- factor(covariates[,i])
  is.inflamed <- map$Inflamed == 'Yes'
  
  
  colnames(lmbf[['module']]) <- sapply(strsplit(colnames(lmbf[['module']]),' \\['),'[',1)
  colnames(mbf[['module']]) <- sapply(strsplit(colnames(mbf[['module']]),' \\['),'[',1)
  if(!is.null(collapse.at)){
    names(mbf.collapse$module$groups) <- sapply(strsplit(names(mbf.collapse$module$groups),' \\['),'[',1)
  }
  
  ignore.pathways <- c("Tuberculosis","Transcription machinery","Translation proteins","DNA replication","Function unknown","General function prediction only","Homologous recombination","Membrane and intracellular structural molecules","Transcription machinery","Translation proteins","Base excision repair","Cell cycle - Caulobacter","Chaperones and folding catalysts","Chromosome","DNA repair and recombination proteins","DNA replication proteins","Mismatch repair","Others","Plant-pathogen interaction","Replication, recombination and repair proteins","Ribosome Biogenesis","RNA polymerase","Transcription factors","Translation factors",'Ribosome',"Protein folding and associated processing","RNA degradation","Drug metabolism - other enzymes","Drug metabolism - other enzymes","Nucleotide excision repair","Photosynthesis proteins", "Carbon fixation in photosynthetic organisms","Photosynthesis proteins")
  # 	outlier.pathways <- c('Methane metabolism','Glycolysis / Gluconeogenesis','Phosphotransferase system (PTS)',    'Porphyrin and chlorophyll metabolism',    'Selenocompound metabolism',    'Vitamin B6 metabolism',    'Fatty acid biosynthesis',    'Cysteine and methionine metabolism',    'Riboflavin metabolism',    'ABC transporters',    'Ubiquinone and other terpenoid-quinone biosynthesis',    'Thiamine metabolism',    'Aminoacyl-tRNA biosynthesis',    'Ribosome',    'Amino sugar and nucleotide sugar metabolism',    'Streptomycin biosynthesis',    'Glutathione metabolism',    'Two-component system',    'beta-Alanine metabolism',    'Carbon fixation in photosynthetic organisms',    'Biotin metabolism',    'Purine metabolism',    'Arginine and proline metabolism',    'D-Glutamine and D-glutamate metabolism',    'D-Alanine metabolism',    'Pyrimidine metabolism',    'Pantothenate and CoA biosynthesis',    'One carbon pool by folate',    'Lipoic acid metabolism',    'C5-Branched dibasic acid metabolism')
  outlier.pathways <- c('Ubiquinone and other terpenoid-quinone biosynthesis','Drug metabolism - other enzymes','Carbon fixation in photosynthetic organisms','Arachidonic acid metabolism','D-Glutamine and D-glutamate metabolism','Phosphotransferase system (PTS)','Vibrio cholerae pathogenic cycle','Aminobenzoate degradation','ABC transporters','Glutathione metabolism','Two-component system','beta-Alanine metabolism','Lipoic acid metabolism','Methane metabolism','Pathways in cancer','Benzoate degradation','Bacterial chemotaxis','D-Alanine metabolism','Alzheimers disease','Zeatin biosynthesis')
  exclude.pathways <- sort(unique(c(ignore.pathways, outlier.pathways)))
  
  names(full.taxon.names) <- colnames(mb)
  
  
  # upon return, call "attach" on result to access all elements directly
  return(list(
    # data tables
    map=map,
    gx=gx,
    gxx=gxx,
    gxx.null=gxx.null,
    gxx.pc=gxx.pc,
    snps=snps,
    mb=mb,
    otus=otus,
    lmb=lmb,
    mbf=mbf,
    lmbf=lmbf,
    dx=dx, # diversity
    bdiv=bdiv,
    covariates=covariates,
    snps.full=snps.full,
    mb.collapse=mb.collapse,
    mbf.collapse=mbf.collapse,
    
    
    # subject subsets
    has.cd=has.cd,
    has.uc=has.uc,
    has.ibd=has.cd | has.uc,
    is.prism=is.prism,
    is.msh=is.msh,
    is.mgh=is.mgh,
    is.mli=is.mli,
    is.risk=is.risk,
    is.neth=is.neth,
    is.colon=is.colon,
    is.ileum=is.ileum,
    is.pouch=is.pouch,
    is.prepouch=is.prepouch,
    is.rectum=is.rectum,
    is.inflamed=is.inflamed,
    
    # GENETIC SUBSETS
    # subset of snps on important genes
    mb.snps.long=mb.snps.long,
    mb.snps=mb.snps,
    cd.snps=cd.snps,
    uc.snps=uc.snps,
    ibd.snps=ibd.snps,
    nod2.snps.long=nod2.snps.long,
    nod2.snps=nod2.snps,
    fut.snps=fut.snps,
    atg.snps=atg.snps,
    jostins.snps=jostins.snps,
    jakstat.snps=jakstat.snps,
    sensing.snps=sensing.snps,
    daly.snps=daly.snps,
    topgenes=topgenes,
    
    # MICROBIOME SUBSETS
    known.taxa=known.taxa,
    cd.taxa=cd.taxa,
    uc.taxa=uc.taxa,
    ibd.taxa=ibd.taxa,
    full.mb.names=full.taxon.names,
    
    # genetic burden info
    or.snp.ix=or.snp.ix,
    or.subj.ix=or.subj.ix,
    
    cd.scores=cd.scores,
    uc.scores=uc.scores,
    ibd.scores=ibd.scores,
    cd.scores.pathway=cd.scores.pathway,
    uc.scores.pathway=uc.scores.pathway,
    ibd.scores.pathway=ibd.scores.pathway,
    disease.scores=disease.scores,
    pathways=pathways,
    pathways.gene=pathways.gene,
    pathway.sizes=pathway.sizes,
    ignore.pathways=ignore.pathways,
    outlier.pathways=outlier.pathways,
    exclude.pathways=exclude.pathways        
  ))
}