# usage:
# make-study-ID... study1.ped study2.ped > study12.pheno
import sys

# count subjects in each study
study1_ids = []
study2_ids = []
for line in open(sys.argv[1],'U'):
    if not line.startswith('#'):
        line = line[:1000]
        print '\t'.join(line.split()[:2]) + '\t' + '1'
for line in open(sys.argv[2],'U'):
    if not line.startswith('#'):
        line = line[:1000]
        print '\t'.join(line.split()[:2]) + '\t' + '2'


