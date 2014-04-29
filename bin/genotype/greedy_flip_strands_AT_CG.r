# run with
# Rscript greedy_... my.flipscan fliplist.txt
# note: flipscan output from plink
# must be preprocessed with:
# cat my.flipscan | sed -e 's: [ ]*: :g' | sed -e 's:^ ::' > my.flipscan.txt

args <- commandArgs(TRUE)
# 

f <- read.table(args[1],head=T,sep=' ',stringsAsFactors=FALSE)

# sort by number of negative correlations
f <- droplevels(f[f[,'NEG'] > 0,])
f <- f[order(f[,'NEG'],decreasing=TRUE),]
is.at <- (f$A1 == 'A' & f$A2 == 'T') | (f$A1 == 'T' & f$A2 == 'A')
is.cg <- (f$A1 == 'C' & f$A2 == 'G') | (f$A1 == 'G' & f$A2 == 'C')
is.amb <- is.at | is.cg

# cumsum.at.cg <- cumsum(is.at | is.cg)
# plot(cumsum.at.cg,type='l',xlim=c(0,length(is.at)), ylim=c(0,length(is.at))); abline(0,1,lty=2,col='#00000099')
# lines(f$NEG / 18 * 120000,col='blue')
# axis(4,at=(0:18)/18 * 120000,labels=0:18,cex.axis=.7,col='blue')

# greedy strand flip of at and cg
flip <- rep(FALSE,nrow(f))
for(i in 1:nrow(f)){
    if(i %% 1000 == 0) cat(i,' ',sep='')
    if(flip[i]) next
    if(is.amb[i] & f$NEG[i] > 0){
        # if this is anti-correlated with at least one unambiguous snp, flip it
        negsnps <- strsplit(f$NEGSNPS,'\\|')[[1]]
        if(any(!is.amb[match(negsnps,f$SNP)])){
            flip[i] <- TRUE
            # if is ambiguous (AT/CG)
            # subtract this SNP from negative correlation counts of all SNPs
            matches <- grep(f$SNP[i], f$NEGSNPS)
            if(length(matches) > 0){
                f$NEG[matches] <- f$NEG[matches] - 1
            }
        }
    }
}
cat('\n')

sink(args[2])
cat(f$SNP[flip],sep='\n')
sink(NULL)