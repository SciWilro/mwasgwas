# usage
# python convert_to_plink_map.py opticall > plink.map
import sys

fin = open(sys.argv[1],'U')

# skip the header
fin.readline()

count = 0
for line in fin:
    line = line[:1000]
    words = line.split()
    print '\t'.join([words[1], words[0], '0', words[2]])