#!/usr/bin/env python
# usage:
# python prune_plink_data -i basename -ld .5 -o out_basename
import sys
import os
from subprocess import Popen, PIPE, STDOUT
from optparse import OptionParser

def make_option_parser():
    parser = OptionParser(usage="usage: %prog [options] filename",
                          version="%prog 1.0")
    parser.add_option("-i", "--basename",
                      default=None,
                      type='string',
                      help="plink basename 1 (required)")
    parser.add_option("-l", "--ld",
                      default=.8,
                      type='float',
                      help="LD threshold for collapsing (default %default)")
    parser.add_option("-p","--print_only",
                      action="store_true",
                      default=False,
                      help="Only print commands (default %default)")
    parser.add_option("-L","--log_file",
                      type='string',
                      default=None,
                      help="Log file path (default basename-log.txt)")
    parser.add_option("-v","--verbose",
                      action="store_true",
                      default=False,
                      help="Verbose output (default %default)")
    parser.add_option("-g","--suppress_genome_wide_pruning",
                  action="store_true",
                  default=False,
                  help="Do not prune based on genome-wide LD (saves time) (default %default)")
    parser.add_option("-o","--out_basename",
                      type="string",
                      default=None,
                      help="Output basename (default basename-prune)")
    return parser

def run_commands(commands, print_only=False, verbose=True):
    return_vals = []
    outputs = []
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
            outputs.append(stdout + '\n' + stderr)
    return(return_vals, outputs)
    
if __name__ == '__main__':
    parser = make_option_parser()
    (options, args) = parser.parse_args()

    if options.out_basename is None:
        base = os.path.basename(options.basename)
        options.out_basename = base + '-prune'

    if options.log_file is None:
        base = os.path.basename(options.basename)
        options.log_file = base + '-log.txt'

    commands = []

    # get scrolling window pruning inclusion list
    cmd = "time plink --noweb --file " + options.basename +\
        " --indep-pairwise 100 1 %f --out " %(options.ld) + options.out_basename + "-prune-window"
    commands.append(cmd)

    # extract
    cmd = "time plink --noweb --file " + options.basename +\
          " --extract " + options.out_basename + "-prune-window.prune.in " +\
          "--recode --out " + options.out_basename + "-prune-window"
    commands.append(cmd)

    if not options.suppress_genome_wide_pruning:
        # get whole-genome correlation table
        cmd = "time plink --noweb --file " + options.out_basename + "-prune-window " +\
            "--r2 --inter-chr --ld-window-r2 %f --out " %(options.ld) + options.out_basename + "-prune-window-prune-all"
        commands.append(cmd)

        # choose centroids (most connected within each cluster) greedily
        cmd = "time Rscript ~/drive/research/prism/code/cluster-snps.r " +\
             options.out_basename + "-prune-window-prune-all.ld > " +\
             options.out_basename + "-prune-window-prune-all.prune.out"
        commands.append(cmd)

        # exclude
        cmd = "time plink --noweb --file " + options.out_basename + "-prune-window " +\
            "--exclude " + options.out_basename + "-prune-window-prune-all.prune.out " +\
            "--recode --out " + options.out_basename
        commands.append(cmd)
    else:
        cmd = "cp " + options.out_basename + "-prune-window.ped " + options.out_basename + ".ped"
        commands.append(cmd)
        cmd = "cp " + options.out_basename + "-prune-window.map " + options.out_basename + ".map"
        commands.append(cmd)
            
    # run all commands
    return_vals, outputs = run_commands(commands, print_only=options.print_only, verbose=options.verbose)
    
    # delete necessary files
    to_remove = [options.out_basename + '-prune-window-prune-all.ld',
                 options.out_basename + '-prune-window-prune-all.log',
                 options.out_basename + '-prune-window-prune-all.prune.out',
                 options.out_basename + '-prune-window.prune.in',
                 options.out_basename + '-prune-window.prune.out',
                 options.out_basename + '-prune-window.ped',
                 options.out_basename + '-prune-window.map',
                 options.out_basename + '-prune-window.log',
        ]
    

    if not options.print_only:
        for f in to_remove:
            try:
                os.remove(f)
            except OSError:
                print "Could not remove file " + f
        # write log file
        flog = open(options.log_file,'a')
        for i,cmd in enumerate(commands):
            flog.write('\n\n' + cmd + '\n')
            flog.write(outputs[i] + '\n')
            if return_vals[i] != 0:
                print "Warning: command",cmd,"had return value",return_vals[i]
                flog.write("Warning: command" + cmd + "had return value" + str(return_vals[i]) + '\n')
        flog.close()