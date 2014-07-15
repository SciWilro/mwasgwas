#!/usr/bin/env python
# usage:
# python merge_studies.py -i basename1 -I basename2 -o outbasename
import sys
import os
from subprocess import Popen, PIPE, STDOUT
from optparse import OptionParser

def make_option_parser():
    parser = OptionParser(usage="usage: %prog [options] filename",
                          version="%prog 1.0")
    parser.add_option("-i", "--basename1",
                      default=None,
                      type='string',
                      help="plink basename 1 (required)")
    parser.add_option("-I", "--basename2",
                      default=None,
                      type='string',
                      help="plink basename 2 (required)")
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
    parser.add_option("-o","--out_basename",
                      type="string",
                      default=None,
                      help="Output basename (default merge-basename1-basename2)",)
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

    base1 = os.path.basename(options.basename1)
    base2 = os.path.basename(options.basename2)

    if options.out_basename is None:
        options.out_basename = '-'.join(["merge", base1, base2])

    commands = []
    # after preprocessing each data set
    # merge two, perform flipscan
    cmd = "time plink --noweb --file " + options.basename1 + \
          " --merge " + options.basename2 + ".ped " + \
          options.basename2 + ".map" + \
          " --recode --tab --out " + options.out_basename + "-pre-flipscan"
    commands.append(cmd)
    
    # run command; if failure, remove offending snps from both input files and try again
    return_vals = run_commands(commands, print_only=options.print_only, verbose=options.verbose)
    commands = []

    if return_vals[-1] != 0:
        mismatch_list_file = options.out_basename + "-pre-flipscan.missnp"
        print "Merged failed, filtering non-matching snps from",mismatch_list_file
        cmd = "time plink --noweb --file " + options.basename1 + \
            " --exclude " + mismatch_list_file + " --recode --out " + base1 + "-exclude-mismatch"
        commands.append(cmd)
        cmd = "time plink --noweb --file " + options.basename2 + \
            " --exclude " + mismatch_list_file + " --recode --out " + base2 + "-exclude-mismatch"
        commands.append(cmd)
        cmd = "time plink --noweb --file " + base1 + "-exclude-mismatch" + \
              " --merge " + base2 + "-exclude-mismatch.ped " + \
              base2 + "-exclude-mismatch.map " + \
              " --recode --tab --out " + options.out_basename + "-pre-flipscan"
        commands.append(cmd)

    # make study IDs as a phenotype file
    cmd = "python ~/drive/research/prism/src/bin/genotype/make-study-ID-as-pheno.py " + \
          options.basename1 + ".ped " + options.basename2 + ".ped > study12.pheno"
    commands.append(cmd)


    # drop snps with very low geno rate
    cmd = "time plink --noweb --file " + options.out_basename + "-pre-flipscan " + \
          "--geno 0.1 --recode --tab --out " + options.out_basename + "-pre-flipscan-genop1"
    commands.append(cmd)

    # flipscan
    cmd = "time plink --noweb --file " + options.out_basename + "-pre-flipscan-genop1 " + \
          "--flip-scan --out " + options.out_basename + "-pre-flipscan-genop1" + \
          " --pheno study12.pheno"
    commands.append(cmd)

    # munge flipscan plink output for readability
    cmd = "cat " + options.out_basename + "-pre-flipscan-genop1.flipscan | sed -e 's: [ ]*: :g' | sed -e 's:^ ::' > " + \
          options.out_basename + "-pre-flipscan-genop1.flipscan.txt"
    commands.append(cmd)


    # greedy chose snps to flip
    cmd = "time Rscript ~/drive/research/prism/src/bin/genotype/greedy_flip_strands_AT_CG.r " +\
         options.out_basename + "-pre-flipscan-genop1.flipscan.txt " +\
         options.out_basename + "-pre-flipscan-genop1.flipscan.fliplist.txt"
    commands.append(cmd)
    
    # get list of subjects to flip (from study 2)
    cmd = "grep '2$' study12.pheno | cut -f -2 > " + options.out_basename + "-pre-flipscan-genop1.flipscan.fliplist-subjects.txt"
    commands.append(cmd)

    cmd = "time plink --noweb --file " +  options.out_basename + "-pre-flipscan-genop1 " + \
        "--flip " + options.out_basename + "-pre-flipscan-genop1.flipscan.fliplist.txt " +\
        "--flip-subset " + options.out_basename + "-pre-flipscan-genop1.flipscan.fliplist-subjects.txt " +\
        "--tab --recode --out " + options.out_basename
    commands.append(cmd)
    
    # run all commands
    return_vals = run_commands(commands, print_only=options.print_only, verbose=options.verbose)
    
    # delete necessary files
    to_remove = [options.out_basename + '-pre-flipscan.log',
                 options.out_basename + '-pre-flipscan.ped',
                 options.out_basename + '-pre-flipscan.map',
                 options.out_basename + '-pre-flipscan.hh',
                 options.out_basename + '-pre-flipscan-genop1.log',
                 options.out_basename + '-pre-flipscan-genop1.ped',
                 options.out_basename + '-pre-flipscan-genop1.map',
                 options.out_basename + '-pre-flipscan-genop1.hh',
                 options.out_basename + '-pre-flipscan-genop1.flipscan.fliplist-subjects.txt',
                 options.out_basename + '-pre-flipscan-genop1.flipscan',
                 options.out_basename + '-pre-flipscan-genop1.flipscan.txt',
                 'study12.pheno'
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
