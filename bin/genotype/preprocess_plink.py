#!/usr/bin/env python
# deduplicates subjects/snps
# flips strands to universal direction,
# renames snps to universal IDs (within a given build of genome)
#
# optional: filters snps and individuals, 
# 
# usage:
#
# preprocess and filter:
# python preprocess_plink_data.py -i inbasename -o outbasename
#
# preprocess and filter custom:
# python preprocess_plink_data.py -i inbasename -o outbasename --maf .05 --geno .2 --mind .2
#
# preprocess only, no filter:
# python preprocess_plink_data.py -i inbasename -o outbasename --nofilter

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
                      help="plink output basename (default inbasename_preprocessed)")
    parser.add_option("-p","--print_only",
                      action="store_true",
                      default=False,
                      help="Only print commands (default %default)",)
    parser.add_option("-s","--suppress_delete",
                      action="store_true",
                      default=False,
                      help="Only print commands for deletion, don't actually delete files (default %default)",)
    parser.add_option("-d","--suppress_deduplicate",
                      action="store_true",
                      default=True,
                      help="Do not deduplicate samples (default %default)",)
    parser.add_option("--maf",
                      type='float',
                      default=.001,
                      help="Filter SNPs by minor allele frequency (default %default)",)
    parser.add_option("--geno",
                      type='float',
                      default=.2,
                      help="Filter individuals by genotyping rate (default %default)",)
    parser.add_option("--mind",
                      type='float',
                      default=.2,
                      help="Filter SNPs by missingness rate (default %default)",)
    parser.add_option("--nofilter",
                      action="store_true",
                      default=False,
                      help="Suppress all filtering (default %default)",)
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
        outfile = options.inbasename + '_preprocessed'

    commands = []

    # working base name will be tmp
    # recode with tabs
    tmpname = str(random.randint(1,1e9))
    cmd = "plink --noweb --ped " + pedfile + " --map " + mapfile + " --recode --tab --out %s" %tmpname
    commands.append(cmd)
    
    # deduplicate individuals
    if not options.suppress_deduplicate:
        cmd = 'sort -u -k 2,2 %s.ped | grep -v "^#" > %s2.ped ' %(tmpname, tmpname)
        commands.append(cmd)
        commands.append('mv %s2.ped %s.ped' %(tmpname, tmpname))

	
    # if not nofilter
    if not options.nofilter:
        cmd = "time plink --noweb --file %s " %(tmpname)\
                + " --mind " + str(options.mind) \
                + " --maf " + str(options.maf) \
                + " --geno " + str(options.geno) \
                + " --recode --tab --out %s2"  %(tmpname)
        commands.append(cmd)
        commands.append('mv %s2.ped %s.ped; mv %s2.map %s.map'  %(tmpname, tmpname, tmpname, tmpname))

    # replace - N with D I
    cmd = "time plink --noweb --file %s --freq --out %s"  %(tmpname, tmpname)
    commands.append(cmd)
#     cmd = 'grep " - [ ]*[ACTG][ACTG]*\| [ACTG][ACTG]* [ ]*-" %s.frq | sed -e "s: [ ]*: :g" | sed -e "s:^ ::" | cut -f 2-4 -d " " | sed -e "s:\( - [ACGT][ACGT]*\):\1 D I:" | sed -e "s:\( [ACGT][ACGT]* -\):\1 I D:" > rename_allele_list.txt' %(tmpname)
    cmd = 'python ~/drive/prism/code/get_plink_insertion_deletion_rename_list.py < %s.frq > rename_allele_list.txt' %(tmpname)
    commands.append(cmd)

    cmd = 'time plink --noweb --file %s --update-alleles rename_allele_list.txt --out %s2 --recode'  %(tmpname, tmpname)
    commands.append(cmd)
    commands.append('mv %s2.ped %s.ped; mv %s2.map %s.map'  %(tmpname, tmpname, tmpname, tmpname))

    # use freq table to flip strands to universal orientation
    cmd = "time plink --noweb --file %s --freq --out %s"  %(tmpname, tmpname)
    commands.append(cmd)
    cmd = 'grep " C [ ]*T \|  G [ ]*T \|  T [ ]*C \|  T [ ]*G "  %s.frq | sed -e "s: [ ]*: :g" | sed -e "s:^ ::" | cut -f 2 -d " " > fliplist.txt'  %(tmpname)
    commands.append(cmd)
    cmd = 'time plink --noweb --file %s --flip fliplist.txt --recode --out %s2'  %(tmpname, tmpname)
    commands.append(cmd)
    commands.append('mv %s2.ped %s.ped; mv %s2.map %s.map'  %(tmpname, tmpname, tmpname, tmpname))
    
    # deduplicate snps
    cmd = "time plink --noweb --file %s --freq --out %s" %(tmpname, tmpname)
    commands.append(cmd)
    cmd = 'Rscript ~/drive/prism/code/deduplicate_plink_data.r -m %s.map -f %s.frq -o duplicate-snps.txt' %(tmpname, tmpname)
    commands.append(cmd)
    cmd = 'time plink --noweb --file %s --exclude duplicate-snps.txt --out %s2 --recode'  %(tmpname, tmpname)
    commands.append(cmd)
    commands.append('mv %s2.ped %s.ped; mv %s2.map %s.map'  %(tmpname, tmpname, tmpname, tmpname))
    
    # rename snps to universal
    cmd = 'Rscript ~/drive/prism/code/get_new_snp_ids_by_coords_and_alleles.r -m %s.map -f %s.frq -o rename-snps.txt'  %(tmpname, tmpname)
    commands.append(cmd)
    cmd = 'time plink --noweb --file %s --update-map rename-snps.txt --update-name --recode --tab --out %s2'  %(tmpname, tmpname)
    commands.append(cmd)
    
    # rename files
    commands.append('mv %s2.ped ' %(tmpname) + outfile + '.ped')
    commands.append('mv %s2.map ' %(tmpname) + outfile + '.map')

    # run all commands
    return_vals = run_commands(commands, print_only=options.print_only, verbose=options.verbose)
    
    # delete necessary files
    to_remove = ['%s2.log' %(tmpname),
                 '%s2.hh' %(tmpname),
                 '%s.log' %(tmpname),
                 '%s.frq' %(tmpname),
                 '%s.map' %(tmpname),
                 '%s.ped' %(tmpname),
                 '%s.nof' %(tmpname),
                 '%s2.nof' %(tmpname),
                 'fliplist.txt',
                 'rename-snps.txt',
                 'duplicate-snps.txt',
                 'rename_allele_list.txt'
        ]
    
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

