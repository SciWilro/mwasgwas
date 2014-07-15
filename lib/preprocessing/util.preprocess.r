# note: preference.rank does not work for multiple=TRUE
"unique.biopsy.location.subset" <- function(map,multiple=FALSE,preference.rank=NULL,
		inflamed.separate=TRUE,
        bodysite.order=c('Ileum','Colon','Rectum','Pre_Pouch_Ileum','Pouch'),
        infstates=c('No','Yes')){
    map <- as.data.frame(map)
     mb.ix <- map[,'Have_Microbiome'] == 'Yes'
     gx.ix <- map[,'Have_Genomics'] == 'Yes'
	keep <- rep(FALSE,nrow(map))
	names(keep) <- rownames(map)
	# sort by body.site, then split by person, take first one
	if(!multiple){
		if(inflamed.separate){
			for(infstate in infstates){
				m2 <- map[map$Inflamed==infstate & mb.ix & gx.ix,]
			
				m2$Biopsy_Location_General <- ordered(m2$Biopsy_Location_General,levels=bodysite.order)
				m2 <- m2[order(m2$Biopsy_Location_General),]
				if(!is.null(preference.rank)) m2 <- m2[order(preference.rank),]
				keep.ids <- sapply(split(rownames(m2),m2$Gx_Subject_ID),'[',1)
				keep[keep.ids] <- TRUE
			}
		} else {
			# sort by biopsy, then inflamed		
			m2 <- map[mb.ix & gx.ix,]
			m2$Biopsy_Location_General <- ordered(m2$Biopsy_Location_General,levels=bodysite.order)
			m2 <- m2[order(m2$Biopsy_Location_General),]
			m2$Inflamed <- ordered(m2$Inflamed,levels=infstates)
			m2 <- m2[order(m2$Inflamed),]
			if(!is.null(preference.rank)) m2 <- m2[order(preference.rank),]
			keep.ids <- sapply(split(rownames(m2),m2$Gx_Subject_ID),'[',1)
		}		
		keep[keep.ids] <- TRUE
		return(keep)
	}
	
    bodysite.ixs <- list()
    for(bs in bodysite.order){
        bodysite.ixs[[bs]] <- map[,'Biopsy_Location_General'] == bs
    }
    keep <- rep(FALSE,nrow(map))
    for(subj in unique(drop.na(map[,'Gx_Subject_ID']))){
        
        if(inflamed.separate){
            # keep one biopsy from inflamed and one from non-inflamed.
            for(infstate in infstates){
                ix <- map[,'Gx_Subject_ID'] == subj & map[,'Inflamed'] == infstate
                ix <- ix & mb.ix & gx.ix
                ix[is.na(ix)] <- FALSE
                if(any(ix)){
                    for(bs in bodysite.order){
                        if(any(bodysite.ixs[[bs]] & ix) & !any(keep[ix])){
                            keep[bodysite.ixs[[bs]] & ix] <- TRUE                            
                        }
                    }
                }
            }
        } else {
            # merge inflamed & non-inflamed
            # choose first listed in infstates first, if not present then second
            # default infstates is 'No','Yes', so non-inflamed chosen first
            already.chosen <- FALSE
            for(infstate in infstates){
                # keep one biopsy from inflamed and one from non-inflamed.
                ix <- map[,'Gx_Subject_ID'] == subj & map[,'Inflamed'] == infstate
                ix <- ix & mb.ix & gx.ix
                ix[is.na(ix)] <- FALSE
                if(any(ix) & !already.chosen){
                    for(bs in bodysite.order){
                        if(any(bodysite.ixs[[bs]] & ix) & !any(keep[ix]))
                            keep[bodysite.ixs[[bs]] & ix] <- TRUE
                    }
                    already.chosen <- TRUE
                }
            }
        }
    }
    return(keep)
}

