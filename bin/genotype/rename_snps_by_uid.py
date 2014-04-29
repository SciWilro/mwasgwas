#!/usr/bin/env python
# renames snps to universal IDs (within a given build of genome)
#
# 
# usage:
#
# preprocess and filter:
# python rename_snps_by_uid.py -i inbasename -o outbasename

import sys
import os
import random
from subprocess import Popen, PIPE, STDOUT
from optparse import OptionParser

def make_option_parser():
    parser = OptionParser(usage="usage: %prog [options] filename",
                          version="%prog 1.0")
    parser.add_option("-i", "--inbasename",
                      default=None,
                      type='string',
                      help="plink input basename (required)")
    parser.add_option("-o", "--outbasename",
                      default=None,
                      type='string',
                      help="plink output basename (default inbasename_uid)")
    parser.add_option("-p","--print_only",
                      action="store_true",
                      default=False,
                      help="Only print commands (default %default)",)
    parser.add_option("-s","--suppress_delete",
                      action="store_true",
                      default=False,
                      help="Only print commands for deletion, don't actually delete files (default %default)",)
    parser.add_option("-v","--verbose",
                      action="store_true",
                      default=False,
                      help="Verbose output (default %default)",)
    return parser

def run_commands(commands, print_only=False, verbose=True):
    return_vals = []
    # run all commands
    for cmd in commands:
        print cmd
        if not print_only:
            proc = Popen(cmd,shell=True,universal_newlines=True,stdout=PIPE,stderr=PIPE)
            stdout, stderr = proc.communicate()
            if verbose:
                print stdout
                print stderr
            return_vals.append(proc.returncode)
    return(return_vals)
    
if __name__ == '__main__':
    parser = make_option_parser()
    (options, args) = parser.parse_args()
    
    pedfile = options.inbasename + '.ped'
    mapfile = options.inbasename + '.map'
    outfile = options.outbasename
    if outfile is None:
        outfile = options.inbasename + '_uid'

    commands = []

    # working base name will be tmp
    # recode with tabs
    tmpname = str(random.randint(1,1e9))
    # get allele freqs
    cmd = "time plink --noweb --file %s --freq --out %s" %(options.inbasename, options.inbasename)
    commands.append(cmd)

    # rename snps to universal
    cmd = 'Rscript ~/drive/prism/code/get_new_snp_ids_by_coords_and_alleles.r -m %s.map -f %s.frq -o rename-snps.txt'  %(options.inbasename, options.inbasename)
    commands.append(cmd)
    cmd = 'time plink --noweb --file %s --update-map rename-snps.txt --update-name --recode --tab --out %s'  %(options.inbasename, outfile)
    commands.append(cmd)
    
    # run all commands
    return_vals = run_commands(commands, print_only=options.print_only, verbose=options.verbose)
    
    # delete necessary files
    to_remove = ['rename-snps.txt']
    
    for f in to_remove:
        if options.suppress_delete: 
            print "rm",f
        else:
            try:
                os.remove(f)
            except OSError:
                print "Could not remove file " + f

    if not options.print_only:
        for i,cmd in enumerate(commands):
            if return_vals[i] != 0:
                print "Warning: command",cmd,"had return value",return_vals[i]

