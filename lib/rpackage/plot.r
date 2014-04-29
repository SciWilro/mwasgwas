# plots an area chart of mixtures
# where mixtures are sorted by JS distance PC 1
# alpha is transparency: 00 to FF
"sorted.areachart" <- function(x,alpha='FF',add=FALSE, ylim=NULL, fill=TRUE,n.collapse=5){
    # note: if x has negatives, then it's log values; re-normalize between 0 and 1
    if(min(x) < 0){
        for(i in 1:nrow(x)){
            x[i, ] <- x[i,] - min(x[i,])
            x[i,] <- x[i,] / sum(x[i,])
        }
    }

    # collapse x by top 5 features
    largest.features <- order(colMeans(x),decreasing=TRUE)[1:min(n.collapse,ncol(x)-2)]
    x <- cbind(rowSums(x[,-largest.features]), x[,largest.features])
    x <- x[,rev(order(colMeans(x)))]
    cols <- c(bong.cols,brewer.pal(9,'Set1'))
    cols <- sprintf('%s%s',cols,alpha)
    plot.ix <- order(cmdscale(jsdist(cbind(x,1-rowSums(x))),1))
#     plot.ix <- order(x[,1])
#     barplot(t(x[plot.ix,]),beside=FALSE,
#         las=2,cex.names=.25,border=NA,space=0,col=cols,...)
    
    
    newx <- rbind(rep(0,nrow(x)),t(x[plot.ix,]))
    newx <- apply(newx,2,cumsum)
    # determine x-axis positions of samples from 0 to 1
    at <- seq(0,1,length=2 * nrow(x) + 1)[2 * (1:nrow(x))]
    if(is.null(ylim)) ylim <- c(0,max(newx))
    if(!add) plot(0,0,type='n',xlab='',ylab='',xlim=c(0,1),ylim=ylim)
#     browser()
    for(i in 2:nrow(newx)){
#         if(i==2) browser()
        if(fill){
            polygon(c(at,rev(at)),c(newx[i,],newx[i-1,ncol(newx):1]),col=cols[(i-1) %% length(cols) + 1],border=NA)
        } else {
            lines(at,newx[i,],col=cols[(i-1) %% length(cols) + 1])
        }
    }

}


"multi.stripchart" <- function(){
# abbreviate names
for(plot.type in c('up','down')){
    abbrev.module.names <- colnames(mbf[['module']])
    names(abbrev.module.names) <- abbrev.module.names
    names.to.change <-  c('M00027__GABA (gamma-Aminobutyrate) shunt','M00198__Putative sn-glycerol-phosphate transport system','M00234__Cystine transport system','M00185__Sulfate transport system','M00270__PTS system, trehalose-specific II component','M00267__PTS system, N-acetylglucosamine-specific II component','M00080__Lipopolysaccharide biosynthesis, inner core => outer core => O-antigen','M00176__Sulfur reduction, sulfate => H2S','M00223__Phosphonate transport system','M00134__Polyamine biosynthesis, arginine => ornithine => putrescine')
    abbrev.module.names[names.to.change] <- c('GABA shunt','sn-glycerol-phosphate tr.','Cystine tr.','Sulfate tr.','PTS sys, trehalose','PTS sys, N-acetylglcsmn','Lipopolysaccharides','Sulfate => H2S','Phosphonate tr.','Arg => putrescine')
    
    res <- res.joint.module$module
    if(plot.type == 'up'){
        signif <- res$qvalue <= .05 & res$direction > 0
    } else {
        signif <- res$qvalue <= .05 & res$direction < 0
    }
    hitlist.full.names <- res[signif,2]
    
    for(i in 1:nrow(res)) res[i,2] <- abbrev.module.names[res[i,2]]
    hitlist <- res[signif,2]
    
    ix <- map$Age >= 18 & map$Age <= 80
    resid <- apply(mbf[['module']][ix,hitlist.full.names], 2,
                get.residuals,x=covariates[ix,],
                              outlier.range=3,
                              do.power.transform=TRUE);
    
    # standardize residuals for display
    for(i in 1:ncol(resid)){
        resid[!is.na(resid[,i]),i] <- drop.na(resid[,i]) / sd(drop.na(resid[,i]))
    }
    
    x <- data.frame(
        abund=as.numeric(resid),
        NOD2=as.factor(rep(rowSums(gx[ix,daly.snps]),length(hitlist))),
         mwas=factor(rep(hitlist,each=nrow(mbf[['module']][ix,])),levels=hitlist),
        variant=as.factor(apply(gx[ix,daly.snps],1,function(xx) ifelse(any(xx > 0),which.max(xx),0)))
    #     updown=as.factor(rep(ifelse(res$direction[signif] > 0,'Up','Down'),each=nrow(resid)))
    )
    
    
    alphas <- c('33',rep('bb',6))
    colours <- sprintf('%s%s',
                      c('#444444','#21FFFF','#FB0007','#1AC805','#FFFF0B','#FA00FF','#0000FF'),
                      alphas)
    
    
    strip.font.size <- ifelse(plot.type=='up',.9,.3)
    plot.nrow <- ifelse(plot.type=='up',2,4)
    p <- ggplot(x, aes(x=NOD2,y=abund,colour=variant))
    p <- p + geom_jitter(position=position_jitter(width=1/4),
             na.rm=TRUE)
    p <- p + facet_wrap(~ mwas, nrow=plot.nrow,scales='fixed')
    p <- p + geom_violin(fill='#00000000',
            na.rm=TRUE,
            scale='width',adjust=2,
            trim=TRUE,colour='#00000066')
    p <- p + theme(strip.text.x = element_text(size = rel(strip.font.size)))
    p <- p + scale_colour_manual(values=colours,
                    labels = c('None',rownames(snps)[daly.snps]))
    p <- p + xlab('NOD2 risk allele dosage')
    p <- p + ylab('Residual relative abundance (standardized)')
    
    if(plot.type == 'up'){
        width <- 10
        height <- 6
    } else {
        width <- 11
        height <- 8.5
    }
#     pdf(sprintf('modules.%s.pdf',plot.type),width=width,height=height)
    print(p)
#     dev.off()
}







}