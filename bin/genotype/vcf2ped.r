# usage:
# R --slave --args vcffile
args <- commandArgs(trailing=TRUE)
prefix <- paste(strsplit(args[1],'\\.')[[1]][-length(strsplit(args[1],'\\.')[[1]])],collapse='')


fullx <- read.table(args[1],sep='\t',head=T,skip=34,comment='',check=F,row=NULL)
fullx <- apply(fullx,2,as.character)
rownames(fullx) <- fullx[,3]
x <- fullx[,-(1:9)]
newx <- x

# make a new table containing the genotype in ped format
for(i in 1:nrow(x)){
    ref <- fullx[i,4]
    alt <- strsplit(fullx[i,5],",")[[1]]
    chars <- c(ref, alt)
    for(j in 1:ncol(x)){
        xx <- strsplit(x[i,j],':')[[1]][1]
        xx <- strsplit(xx,'/')[[1]]
        if(length(xx)==1) xx <- strsplit(xx,'\\|')[[1]]
        
        # missing genotype
        if(xx[1] == '.'){
            newx[i,j] <- "0 0"
        } else {
            # non-missing genotype
            xx <- sort(as.numeric(xx)) # should be 0 for ref, 1 for first alt
            newx[i,j] <- paste(chars[xx+1],collapse=' ')
        }
    }
}

# add columsn for ped format: Family ID, Individual ID,
# Paternal ID, Maternal ID, Sex (1=male; 2=female; other=unknown), Phenotype
newx <- t(newx)
newx <- cbind(rownames(newx), rep('0',nrow(newx)), 
            rep('0',nrow(newx)), rep('other',nrow(newx)),
            rep('1',nrow(newx)), newx)
colnames(newx)[1:5] <- c('Individual ID','Paternal ID','Maternal ID','Sex','Phenotype')

sink(sprintf('%s.ped',prefix))
cat('# Analysis Set: HMP shotgun genomics\n')
cat('#\nFamily ID\t')
write.table(newx,sep='\t',quote=F)
sink(NULL)

# remove extra spaces from fullx columns
for(i in 1:3) fullx[,i] <- gsub(' ','', fullx[,i])
map <- cbind(fullx[,1], fullx[,3], rep('0',nrow(fullx)), fullx[,2])
write.table(map,sprintf('%s.map',prefix),sep='\t',quote=F,col.names=FALSE, row.names=FALSE)