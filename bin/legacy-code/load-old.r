source('~/drive/prism/code/plink.r')
source('~/drive/prism/code/util.load.r')
source('~/drive/prism/code/util.r')
source('~/drive/prism/code/non-parametric-tests.r')
library('optparse',quietly=TRUE)

# set up and parse command-line options
option_list <- list(
    make_option(c("-v", "--verbose"), action="store_true", default=FALSE,
        help="Print extra output [default %default]"),
    make_option(c("-L", "--load_gwas"), type='character', default=NULL,
        help="Preload gwas object from given file [default %default]"),
    make_option(c("-S", "--save_gwas"), type='character', default=NULL,
        help="Save gwas object to given file [default %default]"),
    make_option(c("-O", "--opticall"), type="character", default=NULL,
        help="Path to opticall file [default %default]",
        metavar="path"),
    make_option(c("-p", "--ped"), type="character", default=NULL,
        help="Path to ped file [default %default]",
        metavar="path"),
    make_option(c("-M", "--map"), type="character", default=NULL,
        help="Path to plink map file [default %default]",
        metavar="path"),
    make_option(c("-t", "--taxon_table"), type="character", default=NULL,
        help="Path to tab-delimited taxon table [default %default]",
        metavar="path"),
    make_option(c("-r", "--variant_subset"), type="character", default=NULL,
        help="Path to file listing subset of RS ID's [default %default]",
        metavar="path"),
    make_option(c("-T", "--taxon_subset"), type="character", default=NULL,
        help="Path to file listing subset of taxa [default %default]",
        metavar="path"),
    make_option(c("-m", "--metadata"), type="character", default=NULL,
        help="Path to metadata file [default %default]",
        metavar="path"),
    make_option(c("-s", "--states"), type="character", default=NULL,
        help="Only keep samples with these values in metadata, e.g. 'GENDER:Male;COUNTRY:*,!USA'; if NULL, includes only samples in metadata file [default %default]",
        metavar="state_list"),
    make_option(c("-G", "--groupfile"), type="character", default=NULL,
        help="Also test for grouping of snps in file [default %default]",
        metavar="path"),
    make_option(c("-o", "--outdir"), type="character", default=NULL,
        help=" [default %default]",
        metavar="path")  
)
args <- parse_args(OptionParser(option_list = option_list))

# create outdir
if (!file.exists(args$outdir)) dir.create(args$outdir)

# pre-load list of acceptable samples
if(!is.null(args$metadata)){
    metadata <- load.qiime.mapping.file(args$metadata)
    if(!is.null(args$states)) {
        sample.ids <- extract.samples.by.metadata(metadata, args$states)
        metadata <- droplevels(metadata[sample.ids,,drop=F])
    } else {
        sample.ids <- rownames(metadata)
    }
} else {
    sample.ids <- NULL
    metadata <- NULL
}

# load microbiome with subset taxa and samples
if(args$verbose) cat('Loading mwas table...\n')
mwas <- load.taxon.table(args$taxon_table, taxon.subset=args$taxon_subset, sample.subset=sample.ids)

# load gwas data with subset variants and samples
if(is.null(args$load_gwas)){
    if(args$verbose) cat('Loading gwas table from scratch...\n')
    if(is.null(args$ped)){
        # must be opticall
        gwas <- load.opticall.file(args$opticall, sample.subset=sample.ids, variant.subset=args$variant_subset, verbose=args$verbose)
    } else {
        gwas <- load.ped.file(args$ped, args$map, sample.subset=sample.ids, variant.subset=args$variant_subset, verbose=args$verbose)
    }
} else {
    if(args$verbose) cat(sprintf('Loading gwas table from %s...\n',args$load_gwas))
    load(args$load_gwas)
}

# Extract requested subset of samples
if(args$verbose) cat('Preprocess data...\n')
res <- remove.nonoverlapping.samples(gwas, mwas, metadata, sample.ids)
gwas <- res$gwas
mwas <- res$mwas
metadata <- res$metadata;

# Drop zero-variance features
mwas <- drop.empty.mwas.features(mwas)    
gwas <- drop.empty.gwas.features(gwas)
if(args$verbose) cat(sprintf('%d samples, %d MWAS features, %d GWAS features remain after dropping empty features.\n', nrow(mwas), ncol(mwas), ncol(gwas$ped)-6))

# preprocess
mwas.keepix <- which(colSums(mwas>0) >= 5 & apply(mwas,2,max) > .01)
gwas.keepix <- which(colSums(gwas$counts > 0) >= 5)

pca <- princomp(mwas[,order(colSums(mwas>0),decreasing=TRUE)[1:min(c(ncol(mwas),(.8 * nrow(mwas))))]])
lmwas <- mwas; lmwas[lmwas==0] <- min(lmwas[lmwas > 0])/2; lmwas <- log10(lmwas)

ichipmap <- read.table('~/drive/research/prism/data/genetic/2012-11-08-convert-SNPs/ichip2dbsnp_map.xls',sep='\t',head=T)

LOAD.DISTANCES <- FALSE
LOAD.DIVERSITY <- FALSE
RUN.GX.BURDEN.ANALYSIS <- FALSE

if(LOAD.DISTANCES){
    cat('Loading distances...\n')
    dw <- load.qiime.distance.matrix('~/drive/research/prism/data/microbiome/mine/PRISM-nov/beta-rare900/weighted_unifrac_otus-rare900-bysubj.txt')
    duw <- load.qiime.distance.matrix('~/drive/research/prism/data/microbiome/mine/PRISM-nov/beta-rare900/unweighted_unifrac_otus-rare900-bysubj.txt')

    ix <- intersect(rownames(mwas), rownames(dw))
    mwas <- mwas[ix,]
    gwas$counts <- gwas$counts[ix,]
    gwas$ped <- gwas$ped[ix,]
    metadata <- droplevels(metadata[ix,])
    dw <- dw[ix,ix]
    duw <- duw[ix,ix]
    
    if(nrow(mwas) < 4) stop(sprintf('Only %d samples after filtering.', nrow(mwas)))
    if(args$verbose) cat(sprintf('%d samples, %d MWAS features, %d GWAS features remain after dropping non-overlapping samples.\n', nrow(mwas), ncol(mwas), ncol(gwas$ped)-6))
    
    # Drop zero-variance features
    mwas <- drop.empty.mwas.features(mwas)    
    gwas <- drop.empty.gwas.features(gwas)

    if(args$verbose) cat(sprintf('%d samples, %d MWAS features, %d GWAS features remain after dropping empty features.\n', nrow(mwas), ncol(mwas), ncol(gwas$ped)-6))
    
    # save GWAS object if requested
    if(!is.null(args$save_gwas)){
        if(args$verbose) cat(sprintf('Saving gwas table to %s...\n',args$save_gwas))
        save(gwas,file=args$save_gwas)
    }
        
    # preprocess
    mwas.keepix <- which(colSums(mwas>0) >= 5 & apply(mwas,2,max) > .01)
    gwas.keepix <- which(colSums(gwas$counts > 0) >= 5)
    
    pcw <- cmdscale(dw,2)
    pcuw <- cmdscale(duw,2)
}

if(LOAD.DIVERSITY){
    # load alpha diversity
    adiv <- read.table('~/drive/research/prism/data/microbiome/mine/PRISM-nov/alpha-div.txt',sep='\t',head=T,row=1)
    adiv <- adiv[rownames(mwas),]
    shannon <- apply(mwas,1,function(xx) -sum(xx[xx>0] * log(xx[xx>0])))
}


# load study table and mapping from one to the other
# pvals <- apply(lmwas[cd.ix,],2,function(xx) cor.test(cd.scores, xx)$p.value)
# plot(cd.scores, mwas[cd.ix,'k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Lachnospiraceae;g__Roseburia'])
# cor.test(cd.scores, cmdscale(dw,1)[cd.ix])
# cor.test(cd.scores, cmdscale(duw,1)[cd.ix])
# cor.test(cd.scores, princomp(mwas[cd.ix,])$scores[cd.ix,1])
# mwas <- mwas[,taxon.subset]
lmwas <- mwas; lmwas[lmwas==0] <- min(lmwas[lmwas > 0])/2; lmwas <- log10(lmwas)

if(RUN.GX.BURDEN.ANALYSIS){
    # get mapping from Daly variants to ichip probes
    ichipmap <- read.table('~/drive/research/prism/data/genetic/2012-11-08-convert-SNPs/ichip2dbsnp_map.xls',sep='\t',head=T)
    or.table <- read.table('../../data/genetic/table2-v03052012-JCB-odds-ratios.txt',sep='\t',head=T,row=1)
    gwas.ids <- colnames(gwas$counts)
    study.ids <- as.character(ichipmap[match(colnames(gwas$counts), ichipmap[,1]),'rsID'])
    common <- match(intersect(rownames(or.table),study.ids),study.ids)
    gwas.ids <- gwas.ids[common]
    study.ids <- study.ids[common]
    gwas2study <- study.ids; names(gwas2study) <- gwas.ids
    study2gwas <- names(gwas2study); names(study2gwas) <- gwas2study
    
    
    # correct risk allele counts based on Daly risk alleles
    for(i in 1:length(gwas2study)){
#         if(is.null(args$ped)){
#             # opticall
#             first.allele <- character(nrow(gwas$counts))
#             first.allele[apply(gwas$counts,2,
#         }
        first.allele <- sapply(strsplit(gwas$ped[,study2gwas[i]],' '),'[',1)
        second.allele <- sapply(strsplit(gwas$ped[,study2gwas[i]],' '),'[',2)
        homo.ix <- which(first.allele == second.allele)
        if((gwas$counts[homo.ix[1],study2gwas[i]] == 2 && first.allele[homo.ix[1]] == or.table[gwas2study[i],'IC_nonrisk']) 
            || (gwas$counts[homo.ix[1],study2gwas[i]] == 0 && first.allele[homo.ix[1]] == or.table[gwas2study[i],'IC_risk'])) 
            gwas$counts[,study2gwas[i]] <- 2 - gwas$counts[,study2gwas[i]]
    }
    
    cd.variants <- ichipmap[match(rownames(or.table)[or.table$AssociationType != 'UC'],ichipmap[,'rsID']),1]
    uc.variants <- ichipmap[match(rownames(or.table)[or.table$AssociationType != 'CD'],ichipmap[,'rsID']),1]
    ic.variants <- ichipmap[match(rownames(or.table)[or.table$AssociationType == 'IBD'],ichipmap[,'rsID']),1]
    beta.cd <- or.table[as.character(ichipmap[match(intersect(cd.variants,colnames(gwas$counts)),ichipmap[,1]),'rsID']),'OR_CD']
    beta.cd <- or.table[as.character(ichipmap[match(intersect(cd.variants,colnames(gwas$counts)),ichipmap[,1]),'rsID']),'OR_CD']
    # cd.ix <- which(metadata$Diagnosis_0_CD_1_UC_2_IC_3_HC_4_DC == 0)
    # uc.ix <- which(metadata$Diagnosis_0_CD_1_UC_2_IC_3_HC_4_DC == 1)
    cd.ix <- which(metadata$Disease == 'CD')
    uc.ix <- which(metadata$Disease == 'UC')
    X.cd <- gwas$counts[cd.ix,intersect(cd.variants,colnames(gwas$counts))]
    cd.scores <- as.numeric(beta.cd %*% t(X.cd))
}



# COLLAPSE by gene and pathway
# load ichipmap
# ichipmap <- fast.read.table('../../data/genetic/2012-11-08-convert-SNPs/ichip2dbsnp_map-full.txt')
# ichipmap <- ichipmap[ichipmap[,'distance'] == '0',]

#gwas$counts.by.gene.all <- collapse.ichip.by.gene(gwas$counts, ichipmap, fxn.classes='all')
#gwas$counts.by.gene.missense <- collapse.ichip.by.gene(gwas$counts, ichipmap, fxn.classes='missense')
#gwas$counts.by.gene.default <- collapse.ichip.by.gene(gwas$counts, ichipmap, fxn.classes=NULL)

# collapse by pathway
# pathway.strings <- scan('~/drive/research/prism/data/genetic/pathways/c5.bp.v3.0.symbols.gmt.txt',sep='\n',what='character')
# pathway.strings <- sapply(pathway.strings,function(xx) strsplit(xx,'\t')[[1]])
# pathways <- list(); for(i in 1:length(pathway.strings)){pathways[[pathway.strings[[i]][1]]] <- pathway.strings[[i]][-(1:2)]}
# pathway.counts <- matrix(0,nrow=nrow(gwas$counts.by.gene.missense), ncol=length(pathways)); colnames(pathway.counts) <- names(pathways); rownames(pathway.counts) <- rownames(gwas$counts)
# for(i in 1:length(pathways)){common <- intersect(colnames(gwas$counts.by.gene.missense),pathways[[i]]); if(length(common) > 0) pathway.counts[,i] <- rowSums(gwas$counts.by.gene.missense[,common,drop=FALSE])}
# gwas$counts.by.pathway <- pathway.counts

# if(!is.null(args$save_gwas)){
#     if(args$verbose) cat(sprintf('Saving gwas table to %s...\n',args$save_gwas))
#     save(gwas,file=args$save_gwas)
# }


