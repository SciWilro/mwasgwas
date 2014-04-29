# usage
# make_parallel_dge_test_commmands.py inputdir
# input dir should have taxa.txt, md.txt, gx.txt, gx.annot.txt
import sys

bsub = 'bsub -o maketable.lsf -q hour -R "rusage[mem=32]" "'
basecmd = 'Rscript $MWAS_GWAS_DIR/src/bin/run.dge.test.r -i input/taxa.txt -m input/md.txt -x input/gx.txt -a input/gx.annot.txt '

# read in headers from gx file
line = open('input/gx.txt','r').readline()
snps = line.strip().split('\t')[1:]

for snp in snps:
    cmd = 'bsub -o dge_%s.lsf -q hour -R "rusage[mem=4]" "' %(snp)
    cmd += basecmd + ' -o res_' + snp + '.txt -X ' + snp + '"'
    print cmd
