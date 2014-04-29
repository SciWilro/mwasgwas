# converts vcf to transposed ped
# usage:
# vcf2Tped.py tpedbasename < vcffile 
import sys
import os
pedbase = sys.argv[1]

count = 0
progress_every = 10000
pedf = open(pedbase + '.tped','w')

snpcount = 0
for line in sys.stdin:
    line = line.strip()
    if line.startswith('##'):
        pedf.write(line + '\n')
    elif line.startswith('#'):
        words = line.split('\t')
        # new header
        # chr snpid 0 pos subj1 subj2 ...
        subject_ids = words[9:]
        pedf.write('\t'.join(['#CHR','ID','DIST','POS'] + words[9:]) + '\n')
    else:
        if (count + 1) % progress_every == 0:
            sys.stderr.write(str(count+1) + '\t')
            sys.stderr.flush()
        count += 1
        words = line.split('\t')
        chr = words[0]
        pos = words[1]
        _id = words[2]
        ref = words[3] # ref allele
        alt = words[4].split(',') # alt allele
        # if _id is '.' give sequential id
        if _id == '.':
            _id = 'v_%010d' %(snpcount)
            snpcount += 1
        out = [chr,_id,'0',pos]
        alleles = [ref] + alt
        
        for i in xrange(9,len(words)):
            gt = words[i].split(':')[0]
            gt = gt.split('/')
            if(len(gt) == 1): gt = gt.split('|')
            if(len(gt) != 2): sys.stderr.write('Warning: snp ' + _id + ' subj ' + subject_ids[i-9] + ' does not have 2 genotypes: ' + '/'.join(gt) + '\n') 
            if gt[0] == '.':
                gt = '0 0'
            else:
                gt = alleles[int(gt[0])] + ' ' + alleles[int(gt[1])]
            out.append(gt)
        pedf.write('\t'.join(out) + '\n')
pedf.close()
if count > progress_every: sys.stderr.write('\n')

# write fam file
famf = open(pedbase + '.tfam','w')
for subject in subject_ids:
    famf.write('\t'.join([subject,subject,'0','0','0','0']) + '\n')
famf.close()

