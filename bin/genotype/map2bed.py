# converts plink map file to liftover plain text format
#
# usage 
# python map2bed.py < infile > outfile
import sys

for line in sys.stdin:
    words = line.strip().split('\t')
    if words[0] == '23':
        words[0] = 'X'
    elif words[0] == '24':
        words[0] = 'Y'
    print 'chr%s:%s-%d' %(words[0],words[3],int(words[3]))
