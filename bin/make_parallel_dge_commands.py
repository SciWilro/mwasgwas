# usage
# make_parallel_dge_test_commmands.py inputdir
# input dir should have taxa.txt, md.txt, gx.txt, gx.annot.txt
import sys
from optparse import OptionParser

def make_option_parser():
    parser = OptionParser(usage="usage: %prog [options] filename",
                          version="%prog 1.0")
    parser.add_option("-i", "--microbiome_table",
                      default=None,
                      type='string',
                      help="Path to tab-delimited input table, see input_type [default %default]")
    parser.add_option("-t", "--microbiome_type",
                      default='taxa',
                      type='string',
                      help="Input table type, can be 'taxa', 'functions', or 'otutable' [default %default]")
    parser.add_option("-m", "--metadata",
                      default=None,
                      type='string',
                      help="Path to metadata file [default %default]")
    parser.add_option("-c", "--covariates",
                      default=None,
                      type='string',
                      help="Use these covariates in tests. Comma-separated list of metadata column headers. [default %default]")
    parser.add_option("-x", "--x_table",
                      default=None,
                      type='string',
                      help="Path to 'x' data file (e.g. snps counts) [default %default]")
    parser.add_option("-a", "--x_annotations",
                      default=None,
                      type='string',
                      help="Path to 'x' annotations file (gene labels for snps) [default %default]")
    parser.add_option("-X", "--x_feature",
                      default=None,
                      type='string',
                      help="Test only this feature from x table [default %default]")

    parser.add_option("-t","--trended_disp",
                      action="store_true",
                      default=False,
                      help="Estimate per-sample trended dispersion (slow) [default %default]",)
    parser.add_option("-s","--suppress_delete",
                      action="store_true",
                      default=False,
                      help="Only print commands for deletion, don't actually delete files (default %default)",)
    parser.add_option("-v","--verbose",
                      action="store_true",
                      default=False,
                      help="Verbose output (default %default)",)
    parser.add_option("-o","--out_basename",
                      type="string",
                      default=None,
                      help="Output basename (default merge-basename1-basename2)",)
    return parser

option_list <- list(
  make_option(c("-i", "--microbiome_table"), type="character", default=NULL,
              help="Path to tab-delimited input table, see input_type [default %default]",
              metavar="path"),
  make_option(c("-t", "--microbiome_type"), type="character", default='taxa',
              help="Input table type, can be 'taxa', 'functions', or 'otutable' [default %default]",
              metavar="type"),
  make_option(c("-m", "--metadata"), type="character", default=NULL,
              help="Path to metadata file [default %default]"),
  make_option(c("-c", "--covariates"), type="character", default=NULL,
              help="Use these covariates in tests. Comma-separated list of metadata column headers. [default %default]",
              metavar="covariates"),
  make_option(c("-x", "--x_table"), type="character", default=NULL,
              help="Path to 'x' data file (e.g. snps counts) [default %default]",
              metavar="path"),
  make_option(c("-a", "--x_annotations"), type="character", default=NULL,
              help="Path to 'x' annotations file (gene labels for snps) [default %default]",
              metavar="path"),
  make_option(c("-X", "--x_feature"), type="character", default=NULL,
              help="Test only this feature from x table [default %default]",
              metavar="feature"),
  make_option(c("-T", "--trended_disp"), action="store_true",
              help="Estimate per-sample trended dispersion (slow) [default %default]"),
  make_option(c("-n", "--norm_factor_method"), type="character",default='none',
              help="edgeR calcNormFactors method (none, RLE, or upperquartile) [default %default]"),
  make_option(c("-o", "--outpath"), type="character", default=NULL,
              help="path for saving output [default %default]",
              metavar="path")  
)

bsub = 'bsub -o maketable.lsf -q hour -R "rusage[mem=32]" "'
basecmd = 'Rscript $MWAS_GWAS_DIR/bin/run.dge.test.r -i input/taxa.txt -m input/md.txt -x input/gx.txt -a input/gx.annot.txt '

# read in headers from gx file
line = open('input/gx.txt','r').readline()
snps = line.strip().split('\t')[1:]

for snp in snps:
    cmd = 'bsub -o dge_%s.lsf -q hour -R "rusage[mem=4]" "' %(snp)
    cmd += basecmd + ' -o res_' + snp + '.txt -X ' + snp + '"'
    print cmd
