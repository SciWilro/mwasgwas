# Rscript add-snp-categories.r snp-annotation-table snp-annotation-table-out.txt
args <- commandArgs(trailing=TRUE)
fp.in <- args[1]
fp.out <- args[2]

snps <- read.table(fp.in,sep='\t',head=T,row=1,comment='',stringsAsFactors=FALSE,check=F)

# GENETIC SUBSETS
# subset of snps on important genes
snps$jostins.snps <- grepl('23128233',snps$Reference)
snps$mb.snps <- grepl('ATG16L1|IL23R|NOD2|IRGM',snps$`All Genes`) & grepl('23128233',snps$Reference)
snps$mb.snps <- grepl('ATG16L1|IL23R|NOD2|IRGM|CARD9|FUT2',snps$`All Genes`) & grepl('23128233',snps$Reference)
snps$atg.snps <- grepl('ATG16L1|NOD2|IRGM',snps$`All Genes`) & grepl('23128233',snps$Reference)
snps$atg.snps2 <- grepl('ATG16L1|NOD2|IRGM|ORMDL3',snps$`All Genes`) & grepl('23128233',snps$Reference)
snps$fut.snps <- grepl('FUT2',snps$`All Genes`) & grepl('23128233',snps$Reference)
snps$jostins.nw.snps <- grepl('NOD2|IL10|HCK|DOK3|VDR|SLC11A1|CARD9|LGALS9',snps$`All Genes`)
snps$jakstat.snps <- grepl('JAK|STAT',snps$`All Genes`) & grepl('23128233',snps$Reference)
snps$mb.snps.long <- grepl('ATG16L1|IL23R|IL12B|CARD9|MUC19|TAB1|FUT2|NOD2|IRGM',snps$`All Genes`)  & grepl('23128233',snps$Reference)
snps$cd.snps <- snps$AssociationType == 'CD' | snps$AssociationType == 'IBD'
snps$uc.snps <- snps$AssociationType == 'UC' | snps$AssociationType == 'IBD'
snps$ibd.snps <- snps$AssociationType == 'IBD'
snps$nod2.snps.long <- grepl('NOD2',snps$`All Genes`)
snps$daly.snps <- grepl('Daly',snps$Reference)
snps$nod2.snps <- grepl('NOD2',snps$`All Genes`) & snps$daly.snps 
snps$sensing.snps <- snps$nod2.snps | (grepl('CARD9',snps$`All Genes`) & grepl('23128233',snps$Reference))
snps$rx1.snps <- grepl('SP110',snps$`All Genes`)
snps$rx2.snps <- grepl('ATG16L1|ATG5|IRGM|ANKRD33B,DAP',snps$`All Genes`)
snps$rx3.snps <- grepl('FUT2|PTGER4|MUC19',snps$`All Genes`)
snps$rx4.snps <- grepl('IL23R|IL12B|STAT3|JAK2|CCR6',snps$`All Genes`)
snps$rx5.snps <- grepl('C1orf106',snps$`All Genes`)
snps$rx6.snps <- grepl('IL18RAP',snps$`All Genes`)
snps$rx7.snps <- grepl('MST1',snps$`All Genes`)
snps$rx8.snps <- grepl('CARD9',snps$`All Genes`)

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
topgenes['rs9847710'] <- 'PRKCD'
topgenes['rs670523'] <- 'rs670523'
topgenes['rs4743820'] <- 'NFIL3'
topgenes['rs6716753'] <- 'SP140'
topgenes['rs1728785'] <- 'ZFP90'
topgenes['rs1847472'] <- 'rs1847472'
topgenes <- sapply(strsplit(topgenes,','),'[',1)

snps$topgene <- topgenes

sink(fp.out)
cat('SNP\t')
write.table(snps,sep='\t',quote=F)
sink(NULL)

