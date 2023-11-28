import sys
import pandas as pd
from .multiDict import addtwodimdict


def check_homo(AF):
    if float(AF) == 0.5:
        return "het"
    elif float(AF) == 1:
        return "homo"
    else:
        return AF


def maf2dict(maf_df):
    df_dict = {}
    homo_dict = {}
    for index, row in maf_df.iterrows():
        gene = row['Hugo_Symbol']
        cDNA_Change = row['cDNA_Changes']  # hgvs to cDNA_Changes, Note: not cDNA_Change.
        df_dict.setdefault(gene, []).append(cDNA_Change)
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

