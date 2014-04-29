# Extracts all rows from table that have a value from keyfile in column columnindex
# usage:
# extract_rows_by_column_value.py keyfile columnindex table 
# 
# Skips lines starting with '#'
# Column index starts at 0

import sys

keys = [line.strip() for line in open(sys.argv[1],'U')]
colix = int(sys.argv[2])

for line in open(sys.argv[3],'U'):
    if line.startswith('#'):
        print line.strip()
    else:
        words = line.split('\t')
        if words[colix] in keys:
            print line.strip()

