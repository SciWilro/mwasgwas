# usage
# get_plink.... < plink.frq > rename_allele_list.txt
import sys

for line in sys.stdin:
    words = line.split()
    if words[2]=='-':
        print '\t'.join(words[1:4]) + '\tD\tI'
    if words[3]== '-':
        print '\t'.join(words[1:4]) + '\tI\tD'
