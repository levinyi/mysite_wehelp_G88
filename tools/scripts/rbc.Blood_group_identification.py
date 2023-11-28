import os
import sys
import pandas as pd
import re
import argparse


def _argparse():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('-b', dest='excel_file', action="store", help="input: Blood.Gene.metadata.xlxs")
    parser.add_argument('-p', dest='hpa_database', action="store", help="input: HPA.Gene.cDNA_Changes.xls")
    parser.add_argument('-m', dest='maf_file', action="store", help="input : maf file")
    parser.add_argument('-o', dest="output_file",action="store", help="output: final output table file")
    return parser.parse_args()


def match_allele(cDNA_Changes_list, gene_allele_df, variant_classification_dict, homo_het_set):
    '''
        cDNA_Changes_list = ['c.2561C>T','c.166A>G']
        variant_classification_dict
    '''
    matched_rows = []
    max_match_count = 0

    for index, row in gene_allele_df.iterrows():
        nucleotide_change = re.split(r';|,', row['Nucleotide_change'])
        if set(cDNA_Changes_list) == set(nucleotide_change):
            matched_rows.append(row)
            match_count = len(nucleotide_change)
            if match_count >max_match_count:
                max_match_count = match_count
        elif set(cDNA_Changes_list).issuperset(set(nucleotide_change)):
            matched_rows.append(row)
            match_count = len(nucleotide_change)
            if match_count > max_match_count:
                max_match_count = match_count
    if len(matched_rows) >0:
        max_matched_rows = [row for row in matched_rows if len(re.split(r';|,', row['Nucleotide_change'])) == max_match_count]
        if len(max_matched_rows) >=2:
            allele_name1 = "-"
            allele_name2 = "-"
            return allele_name1, allele_name2
        matched_df = pd.DataFrame(max_matched_rows)
        if 'het' in homo_het_set:
            allele_name1 = matched_df['Allele_name(Het)'].values
            allele_name2 = matched_df['Allele_name(Homo)'].values 
        else:
            allele_name1 = matched_df['Allele_name(Het)'].values
            allele_name2 = matched_df['Allele_name(Het)'].values 

        remainning_changes = set(cDNA_Changes_list) - set(re.split(r';|,', matched_df['Nucleotide_change'].iloc[0]))
        filtered_changes = [change for change in remainning_changes if variant_classification_dict.get(change) != 'Silent']
        if len(filtered_changes) >0:
            matched_rows_f = []
            max_match_count_f = 0
            for index, row in gene_allele_df.iterrows():
                nucleotide_change = re.split(r';|,', row['Nucleotide_change'])
                if set(filtered_changes) == set(nucleotide_change):
                    matched_rows_f.append(row)
                    match_count = len(nucleotide_change)
                    if match_count >max_match_count_f:
                        max_match_count_f = match_count
                    # print("this site is exactly match with ", nucleotide_change, "let's store it to matched_row_f list")
                elif set(filtered_changes).issuperset(set(nucleotide_change)):
                    matched_rows_f.append(row)
                    match_count = len(nucleotide_change)
                    if match_count > max_match_count_f:
                        max_match_count_f = match_count
            if len(matched_rows_f) >0:
                max_matched_rows = [row for row in matched_rows_f if len(re.split(r';|,', row['Nucleotide_change'])) == max_match_count_f]
                matched_df = pd.DataFrame(max_matched_rows)
                allele_name2 = matched_df['Allele_name(Het)'].values
            else:
                allele_name2 = "-"
        return allele_name1, allele_name2
    else:
        matched_df = gene_allele_df[gene_allele_df['Nucleotide_change'] == "Common"]
        allele_name1 = matched_df['Allele_name(Het)'].values
        filtered_changes = [change for change in cDNA_Changes_list if variant_classification_dict.get(change) != 'Silent']
        if len(filtered_changes) == 0:
            allele_name2 = matched_df['Allele_name(Homo)'].values
        else:
            allele_name2 = "-"
        return allele_name1, allele_name2


def remove_duplicates(group):
    '''# 将df中column_name列中同一个Gene下的重复值（除了第一个出现的值）设置为空'''
    mask = group.duplicated()
    group[mask] = ""
    return group

def main():
    parser = _argparse()

    excel_file = parser.excel_file   # input: Blood.Gene.metadata.xlxs
    hpa_db = parser.hpa_database
    maf_file   = parser.maf_file    # input: maf file
    final_output_file = parser.output_file  # output: final ouput table file

    #############################################################
    # deal with MAF file.
    maf_df = pd.read_csv(maf_file, sep="\t")
    maf_df['cDNA_Changes'] = maf_df['cDNA_Changes'].fillna('')  # 将缺失值（NaN）替换为一个空字符串（''）
    maf_df['homo/het'] = maf_df['homo/het'].fillna('')  # 将缺失值（NaN）替换为一个空字符串（''）

    new_cDNA = maf_df.groupby('Gene')['cDNA_Changes'].apply(list).reset_index(name='cDNA_Changes_list')
    new_homo = maf_df.groupby('Gene')['homo/het'].apply(set).reset_index(name='homo/het_set')
    new_df = pd.merge(new_cDNA, new_homo, on='Gene')

    variant_classification_dict = maf_df.set_index('cDNA_Changes')['Variant_Classification'].to_dict()

    ######################################################
    #### for allele name: add two columns: Allele_name1, Allele_name2
    df1 = pd.read_excel(excel_file)  # blood pdf
    new_allele_homo = []
    new_allele_het = []
    for index, row in new_df.iterrows():
        cDNA_Changes_list = row['cDNA_Changes_list']
        homo_het_set = row['homo/het_set']
        gene = row['Gene']
        gene_allele_df = df1[df1['Gene'] == gene]
        # print("Now let's start with ",gene, cDNA_Changes_list)
        if len(cDNA_Changes_list) == 1 and cDNA_Changes_list[0] == '':
            matched_rows = gene_allele_df[gene_allele_df['Nucleotide_change'] == 'Common']
            allele_name1 = matched_rows['Allele_name(Het)'].values
            allele_name2 = matched_rows['Allele_name(Homo)'].values
        else:
            allele_name1, allele_name2 = match_allele(cDNA_Changes_list, gene_allele_df, variant_classification_dict, homo_het_set)
        new_allele_homo.append(allele_name1)
        new_allele_het.append(allele_name2)
    new_df['Allele_name1'] = new_allele_homo
    new_df['Allele_name2'] = new_allele_het

    ################################################
    maf_df = pd.merge(maf_df, new_df, on="Gene")
    # print(maf_df)
    maf_df['Allele_name1'] = maf_df['Allele_name1'].apply(lambda x: "".join(x))
    maf_df['Allele_name2'] = maf_df['Allele_name2'].apply(lambda x: "".join(x))

    # maf_df.loc[maf_df.duplicated(subset='Allele_name1'), 'Allele_name1'] = ""  # bug
    maf_df['Allele_name1'] = maf_df.groupby('Gene',group_keys=False)['Allele_name1'].apply(remove_duplicates)
    maf_df['Allele_name2'] = maf_df.groupby('Gene',group_keys=False)['Allele_name2'].apply(remove_duplicates)

    #################################################
    #  for HPA
    if 'HPA' in maf_df['System'].values:
        hpa_df = pd.read_csv(hpa_db, sep="\t")
        maf_df = pd.merge(maf_df, hpa_df, how="left", on=['System', 'Gene', 'cDNA_Changes'])
        maf_df['Allele_name1'] = maf_df['homo']
        maf_df['Allele_name2'] = maf_df['het']

    #################################################

    maf_df[['Number', 'System', 'Gene', 'PanelDesigned', 'MutationDetected',
        'cDNA_Changes','Variant_Classification', 'homo/het', 'Allele_name1', 'Allele_name2']].to_excel(final_output_file, index=False)
    maf_df[['Number', 'System', 'Gene', 'PanelDesigned', 'MutationDetected',
        'cDNA_Changes','Variant_Classification', 'homo/het', 'Allele_name1', 'Allele_name2']].to_csv(final_output_file.replace("xlsx", "xls"), index=False)
    print(maf_df)

if __name__ == '__main__':
    main()

