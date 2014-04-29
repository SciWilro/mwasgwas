#!/usr/bin/env python
# usage:
# ichip2dbsnp -i ichip.map -d dbsnp_file -o outdir -x
#
# output is a tab-delimited file with these columns
# ichipID, chr_37, pos_37, dbsnpID, distance, dbsnp_pos_37, snpClass, alleles, 
# 
# one ichipID can have multiple entries. 
# snpClass is one of: snp; in-del; heterozygous; microsatellite; named-locus; no-variation; mixed; multinucleotide-polymorphism
# alleles is like "-/C"
# dbsnp_pos_37 is from assembly GRCh37.p5
import sys
import os
from optparse import OptionParser
import numpy as np
import gc

def make_option_parser():
    parser = OptionParser(usage="usage: %prog [options] filename",
                          version="%prog 1.0")
    parser.add_option("-i", "--ichip_map",
                      default=None,
                      type='string',
                      help="ichip plink map file (required)")
    parser.add_option("-x", "--maxdist",
                      default=5,
                      type='int',
                      help="maximum number of +/- basepairs to look for match (default %default)")
    parser.add_option("-d", "--dbsnp_dir",
                      type="string",
                      default='.',
                      help="Single .flat file or directory containing .flat files for all chromosomes (required)")
    return parser

if __name__ == '__main__':
    parser = make_option_parser()
    (opts, args) = parser.parse_args()

    # load all snps
    # dict of dicts
    # dict: {dbsnp_id:{'chr':chr,'pos':pos,'snpClass':snpClass,'alleles':alleles,'locus':locus,'locusID':locusID,'fxnclass':[fxnclass1,fxnclass2]}
    db = {}
    snpid = None
    ichipcount = 0
    linecount = 0
    snpcount = 0

    output = [] # will hold list of lists of columns
    
    # reverse_lookup is a dict of {pos:[snpid1, snpid2,...]}
    # for all snps within maxdist of the position
    reverse_lookup = {}
    
    if os.path.isdir(opts.dbsnp_dir):
        files = os.listdir(opts.dbsnp_dir)
    else:
        files = [opts.dbsnp_dir]

    for f in files:
        if not f.endswith(".flat"):
            continue
        sys.stderr.write(f + '\n')
        sys.stderr.flush()
        for line in open(f,'U'):
            if linecount % 1000000 == 0:
                sys.stderr.write(str(linecount) + ' ')
                sys.stderr.flush()
            linecount += 1
            if line.startswith('rs'):
                ichipcount += 1
                if snpid is not None:
                    # we just finished a snp. check for missing info
                    if not db[snpid].has_key('pos'):
                        sys.stdout.write('Warning: snp ' + snpid + ' has no position.')
                        sys.exit(1)
                    
                    # we just finished a snp. if we're past 10M lines, check for matches
                    # then start over
                    if ichipcount % 500000 == 0:
                        sys.stderr.write('ichipcount reached setpoint. checking for matches.\n')
                        sys.stderr.flush()
                        f_in = open(opts.ichip_map,'U')
                        inlinecount = 0
                        for iline in f_in:
                            # if inlinecount % 100000 == 0:
                            #     sys.stderr.write(str(inlinecount) + ' ')
                            #     sys.stderr.flush()
                            inlinecount += 1
                            words = [word.strip() for word in iline.split('\t')]
                            _id = words[1]
                            chrome = words[0]
                            if chrome == '23':
                                chrome = 'X'
                            elif chrome == '24':
                                chrome = 'Y'
                            pos = int(words[3])

                            # get list of snps in reverse lookup
                            for modpos in xrange(pos - opts.maxdist, pos + opts.maxdist + 1):
                                if not reverse_lookup.has_key(modpos):
                                    continue
                                for snp in reverse_lookup[modpos]:
                                    if chrome == db[snp]['chr']:
                                        posx = db[snp]['pos']
                                        dist = abs(pos - posx)
                                        if abs(pos - posx) <= opts.maxdist:
                                            # add a line to the output
                                            # ichipID, chr_37, pos_37, dbsnpID, distances, dbsnp_pos_37, snpClass, alleles, locus, locusID, fxnclasses
                                            fxnclasses = db[snp]['fxnclass']
                                            if len(fxnclasses) == 0:
                                                fxnclasses = ['NA']
                                            newline = [_id,chrome,str(pos),snp,str(dist),str(posx),
                                                        db[snp]['snpClass'], db[snp]['alleles'], db[snp]['locus'], db[snp]['locusID'],','.join(fxnclasses)
                                                      ]
                                            snpcount += 1
                                            output.append('\t'.join(newline))
                        f_in.close()
                        # sys.stderr.write('\n')
                        db = {}
                        reverse_lookup = {}
                        gc.collect()
                    
                sys.stderr.flush()
                words = [word.strip() for word in line.strip().split('|')]
                snpid = words[0]
                if(db.has_key(snpid)):
                    sys.stdout.write('Warning: snp ' + snpid + ' found more than once.')
                    sys.exit(1)
                db[snpid] = {}
                db[snpid]['snpClass'] = words[3]
                db[snpid]['fxnclass'] = []
                db[snpid]['locus'] = 'NA'
                db[snpid]['locusID'] = 'NA'
            else:
                if snpid is None:
                    continue
                if line.startswith('SNP'):
                    words = [word.strip() for word in line.strip().split('|')]
                    db[snpid]['alleles'] = words[1][(words[1].index('=')+2):-1]
                elif line.startswith('CTG'):
                    words = [word.strip() for word in line.strip().split('|')]
                    if words[1] == "assembly=GRCh37.p5":
                        pos = words[3][(words[3].index('=')+1):]
                        if pos != '?':
                            pos = int(pos)
                            db[snpid]['chr'] = words[2][(words[2].index('=')+1):]
                            db[snpid]['pos'] = pos
                            if not reverse_lookup.has_key(pos):
                                reverse_lookup[pos] = []
                            reverse_lookup[pos].append(snpid)
                elif line.startswith('LOC'):
                    words = [word.strip() for word in line.strip().split('|')]
                    db[snpid]['locus'] = words[1]
                    db[snpid]['locusID'] = words[2][(words[2].index('=')+1):]
                    fxnclass = words[3][(words[3].index('=')+1):]
                    if fxnclass != 'reference':
                        if fxnclass not in db[snpid]['fxnclass']:
                            db[snpid]['fxnclass'].append(fxnclass)
    
    sys.stderr.write("\n")
    sys.stderr.write(str(len(db)) + ' snps found.\n')
    sys.stderr.flush()
    
    # now iterate through illumina file, building output
    # one last time

    if len(db) > 0:
        inlinecount = 0
        for line in open(opts.ichip_map,'U'):
            # if inlinecount % 100 == 0:
            #     sys.stderr.write(str(inlinecount) + ' ')
            #     sys.stderr.flush()
            inlinecount += 1
            words = [word.strip() for word in line.split('\t')]
            _id = words[1]
            chrome = words[0]
            if chrome == '23':
                chrome = 'X'
            elif chrome == '24':
                chrome = 'Y'
            pos = int(words[3])
        
            # get list of snps in reverse lookup
            for modpos in xrange(pos - opts.maxdist, pos + opts.maxdist + 1):
                if not reverse_lookup.has_key(modpos):
                    continue
                for snp in reverse_lookup[modpos]:
                    if chrome == db[snp]['chr']:
                        posx = db[snp]['pos']
                        dist = abs(pos - posx)
                        if abs(pos - posx) <= opts.maxdist:
                            # add a line to the output
                            # ichipID, chr_37, pos_37, dbsnpID, distances, dbsnp_pos_37, snpClass, alleles, locus, locusID, fxnclasses
                            fxnclasses = db[snp]['fxnclass']
                            if len(fxnclasses) == 0:
                                fxnclasses = ['NA']
                            newline = [_id,chrome,str(pos),snp,str(dist),str(posx),
                                        db[snp]['snpClass'], db[snp]['alleles'], db[snp]['locus'], db[snp]['locusID'],','.join(fxnclasses)
                                      ]
                            snpcount += 1
                            output.append('\t'.join(newline))
        # sys.stderr.write('\n')

    sys.stderr.write(str(inlinecount) + ' variants checked\n')
    sys.stderr.write(str(snpcount) + ' SNPs matched\n')
    sys.stderr.write('Writing results to file...\n')
    sys.stderr.flush()
    print '\n'.join(output)

