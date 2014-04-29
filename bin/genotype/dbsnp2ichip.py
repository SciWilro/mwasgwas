#!/usr/bin/env python
# usage:
# dbsnp2ichipgrepc pattern infile
# grepc -d ',' pattern infile
# grepc -f patternfile -d ',' infile
#
# output is a tab-delimited file with these columns
# ichipID, chr_37, pos_37, dbsnpID, dbsnp_pos_37, snpClass, alleles
# 
# snpClass is one of: snp; in-del; heterozygous; microsatellite; named-locus; no-variation; mixed; multinucleotide-polymorphism
# alleles is like "-/C"
# dbsnp_pos_37 is from assembly GRCh37.p5
import sys
from optparse import OptionParser

def make_option_parser():
    parser = OptionParser(usage="usage: %prog [options] filename",
                          version="%prog 1.0")
    parser.add_option("-i", "--ichip_map",
                      default=None,
                      type='string',
                      help="ichip plink map file (required)")
    parser.add_option("-d", "--dbsnp_dir",
                      type=".",
                      default=None,
                      help="Directory containing single concatenated dbsnp flat asn.1 file for all chromosomes (default %default)",)
    parser.add_option("-s","--start_index",
                      type="int",
                      default=1,
                      help="Start at the nth line (default %default)",)
    return parser

if __name__ == '__main__':
    parser = make_option_parser()
    (options, args) = parser.parse_args()

    if options.file is None:
        patterns = [args[0]]
        infile = args[1]
    else:
        patterns = [line.strip() for line in open(options.file,'U')]
        infile = args[0]
    
    linecount = 0
    columnixs = []
    output = []
    for line in open(infile,'U'):
        sys.stderr.write(str(linecount) + ' ')
        sys.stderr.flush()
        if linecount >= options.start_index - 1:
            words = [word.strip() for word in line.split(options.delim)]
            if linecount == options.start_index - 1:
                for pattern in patterns:
                    try:
                        columnix = words.index(pattern)
                        columnixs.append(columnix)
                    except:
                        pass
                if len(columnixs) == 0:
                    break
            newwords = []
            for ix in columnixs:
                newwords.append(words[ix])
            output.append('\t'.join(newwords))
        linecount += 1
    sys.stderr.write('\n')

    print '\n'.join(output)
