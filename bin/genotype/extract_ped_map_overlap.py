# produce ped/map file containing only those probes in common
# extract_ped_map_overlap.py pedfile mapfile outbasename
import sys
import numpy as np

pedfile = sys.argv[1]
mapfile = sys.argv[2]
outfile = sys.argv[3]

# get all map file ids
mapdict = {}
for line in open(mapfile,'U').readlines():
    words = line.strip().split('\t')
    mapdict[words[1]] = line.strip()
mapids = set(mapdict.keys())

# load ped file, indexed by column ID
pedf = open(pedfile,'U')
line = '#'
comments = []
while line.startswith('#'):
    line = pedf.readline().strip()
    comments.append(line)

pedids = comments[-2].split('\t')[6:]

# determining overlap as set
print 'determining overlap as set...'
overlap = mapids.intersection(pedids)
# determining overlap as boolean
print 'determining overlap as boolean...'
pedinmap = np.array([pedid in overlap for pedid in pedids])
outpedf = open(outfile + '.ped','w')

print len(overlap),'ids match of',len(pedids),'in ped and',len(mapids),'in map.'

# write non-header comments
outpedf.write('\n'.join(comments[:-2]) + '\n')

# write header
# metadata fields
line = comments[-2]
words = line.strip().split('\t')
outpedf.write('\t'.join(words[:6]) + '\t')
# data columns
outwords = []

for i,w in enumerate(words[6:]):
    if pedinmap[i]:
        outwords.append(w)
outpedf.write('\t'.join(outwords) + '\n')


print 'writing filtered ped file...'
line = comments[-1]
count = 0
while line != '':
    count += 1
    if count % 100 == 0:
        sys.stdout.write(str(count) + ' ')
        sys.stdout.flush()
    words = line.strip().split('\t')
    line = pedf.readline()

    # if wrong number of words, skip this line and print warning
    if len(words) != len(pedids) + 6:
        print
        print "Warning: line",count,"has only",len(words),"words, should be",len(pedids)
        continue
    # metadata fields
    outpedf.write('\t'.join(words[:6]) + '\t')
    # data columns
    outwords = []

#     outwords += words[6:][pedinmap].tolist()
        
    for i,w in enumerate(words[6:]):
        if pedinmap[i]:
            outwords.append(w)
    outpedf.write('\t'.join(outwords) + '\n')
print
outpedf.close()
pedf.close()

print 'writing filtered map file...'
outmapf = open(outfile + '.map','w')
for i,pedid in enumerate(pedids):
    if pedinmap[i]:
        outmapf.write(mapdict[pedid] + '\n')
outmapf.close()