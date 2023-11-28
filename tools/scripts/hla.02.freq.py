import os
import sys
import re

db_file = sys.argv[1]  
hd_file = sys.argv[2]

db_dict = dict()
with open(db_file, "r") as f:
    for line in f:
        line = line.rstrip("\n")
        gene,count,freq,C = line.split()
        db_dict[gene] = freq

file_name = hd_file.replace("result.txt", "result.with.freq.txt")
with open(hd_file, "r") as f, open(file_name, "w") as output:
    for line in f:
        line = line.rstrip("\n")
        a = line.split()
        for each in a[1:]:
            if each.startswith("HLA"):
                name,b = each.split("-")
                match = re.match(r'(\w+\*\d+\:\d+)', b)
                output.write("{}\t{}\n".format(b, db_dict.get(match.group(1), "")))
