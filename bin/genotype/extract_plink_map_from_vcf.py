# usage
# extract_... < in.vcf > out.map
import sys

count = 0
for line in sys.stdin:
    if line.startswith('#'):
        continue
    count += 1
    if count % 1000000 == 0:
        sys.stderr.write(str(count) + ' ')
    line = line[:200]
    words = line.strip().split('\t')
    sys.stdout.write('\t'.join([words[0],words[2],'0',words[1]]) + '\n')
sys.stderr.write('\n')