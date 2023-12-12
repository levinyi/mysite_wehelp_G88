import re
import sys
import os
from collections import defaultdict
import pandas as pd
from Bio import SeqIO


def deal_fasta(fasta):
    record_dict = {}
    with open(fasta,"r") as f:
        for record in SeqIO.parse(f,"fasta"):
            record_dict[record.id] = str(record.seq)
    return record_dict

def replace_base(cdna_change, homo_het):
    a_dict = {
        ('A','C'): 'M', ('A','C','G'): 'V',
        ('A','G'): 'R', ('A','C','T'): 'H',
        ('A','T'): 'W', ('A','G','T'): 'D',
        ('C','G'): 'S', ('C','G','T'): 'B',
        ('C','T'): 'Y', ('A','G','C','T'):'N',
        ('G','T'): 'K',
    }
    #################################3
    pos = 0
    ref = ''
    allele = ''

    if 'ins' in cdna_change:
        result = re.match(r"c\.(\d+)_?\d*ins([A-Z]+)", cdna_change, re.I)
        if result is not None:
            pos, ref, allele = result.groups()
    elif 'del' in cdna_change:
        '''c.1304_1319delCCACCCCAGAGCCCAC'''
        result = re.match(r"c\.(\d+)_(\d+)?del([A-Z]+)", cdna_change, re.I)
        if result is not None:
            pos, ref, allele = result.groups()
    else:
        result = re.match(r"c\.(\d+)([A-Z]+)>([A-Z]+)", cdna_change, re.I)
        if result is not None:
            pos, ref, allele = result.groups()
    pos = int(pos) - 1
    #################################
    if homo_het == 'homo':
        subt = allele
    else:
        subt = a_dict.get(tuple(sorted((ref, allele))),"Ins/Del") # type: ignore
    return pos, ref, allele, subt


def extract_cds(mutation_dict, gene_dict, cds_file, out_file):
    '''mutation_dict : {transcript_id: {cDNA_changes : het/homo }}'''
    cds_dict = deal_fasta(cds_file)
    
    raw_seq = ''
    with open(out_file, 'w') as f:
        for trans_id, cdna_dict in mutation_dict.items():
            raw_seq = cds_dict[trans_id]
            new_seq = raw_seq
            for changes, homo_het in cdna_dict.items():
                pos, ref, allele, subt = replace_base(changes, homo_het)
                new_seq = new_seq[:pos] + subt + new_seq[pos + len(subt):]
            f.write(f">{gene_dict[trans_id]}\n{new_seq}\n")


def deal_MAF_file(maf_file):
    '''从df中获取Gene，cDNA_Changes，transcript_id, homo/het'''

    df = pd.read_csv(maf_file, sep="\t")
    df = df[df['MutationDetected'] == "YES"]

    mutation_dict = defaultdict(lambda: defaultdict(str))
    gene_dict = {row['Transcript_USED']: row['Gene'] for _, row in df.iterrows()}
    
    for _, row in df.iterrows():
        mutation_dict[row['Transcript_USED']][row['cDNA_Changes']]= row['homo/het']

    return mutation_dict, gene_dict


def main():
    maf_file = sys.argv[1]  # small.add.dbsnp.MAF.xls
    cds_file = sys.argv[2]  # cds.fasta
    out_file = sys.argv[3]  # output.cds.fasta

    mutation_dict, gene_dict = deal_MAF_file(maf_file)
    extract_cds(mutation_dict, gene_dict, cds_file, out_file)


if __name__ == '__main__':
    main()
