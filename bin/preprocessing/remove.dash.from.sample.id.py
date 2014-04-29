# Rscript scriptname mappingfile
# because qiime merge_mapping_files.py converts "." to "-"
source(sprintf('%s/src/lib/rpackage/util.load.r',Sys.getenv('MWAS_GWAS_DIR')))

args <- commandArgs(trailing=TRUE)

map <- load.qiime.mapping.file(args[1])
rownames(map) <- gsub('-','.',rownames(map))
sink(args[1])
cat('#SampleID\t')
write.table(map,sep='\t',quote=F)
sink(NULL)