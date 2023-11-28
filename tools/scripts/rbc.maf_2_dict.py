import sys
import pandas as pd


'''
bug here.
def check_homo(ref, allele1, allele2):
    if allele1 == ref or allele2 == ref:
        return "het"
    else:
        return "homo"
'''

def check_homo(AF):
    if float(AF) == 0.5:
        return "het"
    elif float(AF) == 1:
        return "homo"
    else:
        return AF


def addtwodimdict(thedict, key_a, key_b, val):
    ''' this is a function to add two dimetion dict '''
    if key_a in thedict:
        thedict[key_a].setdefault(key_b, []).append(val)
    else:
        thedict.update({key_a: {key_b: [val]}})
    return thedict


def maf2dict(maf_df):
    # deal with MAF file.
    # df = pd.read_csv(maf_file, sep="\t")
    df_dict = {}
    homo_dict = {}
    for index, row in maf_df.iterrows():
        gene = row['Hugo_Symbol']
        cDNA_Change = row['cDNA_Changes']  # hgvs to cDNA_Changes, Note: not cDNA_Change.
        df_dict.setdefault(gene, []).append(cDNA_Change)
        # check homo or het
        # R1 = row['Reference_Allele']
        # A1 = row['Tumor_Seq_Allele1']
        # A2 = row['Tumor_Seq_Allele2']
        AF = row['AF']
        # if check_homo(R1, A1, A2) == "het":
        if check_homo(AF) == "het":
            addtwodimdict(homo_dict, gene, "het", cDNA_Change)
        elif check_homo(AF) == "homo":
            addtwodimdict(homo_dict, gene, "homo", cDNA_Change)
        else:
            continue
            # print(f"please check AF : {AF}\t{gene}\t{cDNA_Change}")
    return df_dict, homo_dict

