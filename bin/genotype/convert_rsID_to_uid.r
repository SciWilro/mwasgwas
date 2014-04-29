# usage convert_rsID_to_uid.r rsfile.txt > outfile.txt
args <- commandArgs(trailing=TRUE)
icmap <- read.table('~/drive/research/prism/data/genetic/prism/2012-11-08-convert-SNPs/ichip2dbsnp_map-full.txt',sep='\t',head=T,comment='')
icmap <- icmap[icmap[,'distance'] == 0,]
null.list <- scan(args[1],w='c',quiet=TRUE)
matched <- null.list %in% icmap[,'rsID']

icmap <- icmap[match(null.list[matched],icmap[,'rsID']),]
cat(sprintf('%s\tuid_%s_%d',icmap[,'rsID'],icmap[,'chromosome'],icmap[,'ichip_build37_pos']),sep='\n')
#  also print any input IDS that already start with uid_
cat(sprintf('%s\t%s',
        grep('^uid_',null.list,value=TRUE),
        grep('^uid_',null.list,value=TRUE)),sep='\n')

# load custom mapping file
custom.mapping <- read.table('custom_rs2uid_mapping.txt',sep='\t',head=F,row=1)
null.list <- null.list[!matched]
null.list <- grep('^rs',null.list,value=TRUE)

for(x in null.list){
    if(x %in% rownames(custom.mapping)){
        cat(sprintf('%s\t%s\n',x,custom.mapping[x,]))
    }
}

null.list <- setdiff(null.list,rownames(custom.mapping))
if(length(null.list) > 0){
    # print failures
    cat(sprintf('%s\tFailed',null.list),sep='\n')
}
