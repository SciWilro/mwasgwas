# appends sequential numbers to seq ids to make them unique
import sys

ids = set()
dupseqs = {}
count = 0

for line in sys.stdin:
    if line.startswith('>'):
        words = line.strip().split()
        words[0] += '.%012d' %(count)
        print ' '.join(words)
        count += 1
    else:
        print line.strip()

