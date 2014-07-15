# usage python remove_allele_from_snp_rsID.py pedfile
import sys

count = 0
for line in open(sys.argv[1],'U'):
	if count == 0:
		words = line.strip().split('\t')
		for i in xrange(1,len(words)):
			words[i] = words[i].split('_')[0]
		print '\t'.join(words)
	else:
		print line.strip()
	count += 1
