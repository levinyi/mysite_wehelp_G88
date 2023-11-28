import sys
import os

input_files = sys.argv[1:]


def addtwodimdict(thedict, key_a, key_b, val):
    ''' this is a function to add two dimetion dict '''
    if key_a in thedict:
        thedict[key_a].update({key_b: val})
    else:
        thedict.update({key_a: {key_b: val}})
    return thedict


def deal_hla(afile):
    type1 = 'NULL'
    type2 = 'NULL'

    base_name = os.path.basename(afile)
    sample, gene_name, out, txt  = base_name.split(".")
    with open(afile, "r") as f:
        for line in f:
            if line.startswith("IMGT/HLA does not exists"):
                return sample, gene_name, 'NULL', 'NULL'

            line = line.rstrip("\n")
            if line.startswith("HLA gene"):
                a,name = line.split(":")
                name = name.lstrip()
                continue
            if line.startswith("[Type 1]"):
                c = line.split("\t")
                type1 = c[1]
                continue
            if line.startswith("[Type 2]"):
                c = line.split("\t")
                type2 = c[1]
                continue
    return sample, name, type1, type2


hla_types = ["HLA-A","HLA-B","HLA-C","HLA-E","HLA-F","HLA-G","MICA","MICB","HLA-DMA","HLA-DMB",
        "HLA-DOA","HLA-DOB","HLA-DPA1","HLA-DPB1","HLA-DQA1","HLA-DQB1","HLA-DRA","HLA-DRB1",
        "HLA-DRB5","TAP1","TAP2"]


results = {}
for each_file in input_files:
    if each_file.endswith("out.txt"):
        sample, gene_name, type1, type2 = deal_hla(each_file)
        addtwodimdict(results, gene_name, 'type1', type1)
        addtwodimdict(results, gene_name, 'type2', type2)

print("HLA\t{}".format(sample))
for each in hla_types:
    if each in results:
        print("{}-1\t{}".format(each, results[each].get("type1",'NULL')))
        print("{}-2\t{}".format(each, results[each].get("type2",'NULL')))
