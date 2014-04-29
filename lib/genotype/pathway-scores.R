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
  #   	pathways <- rbind(mb.snps.long,pathways)
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
