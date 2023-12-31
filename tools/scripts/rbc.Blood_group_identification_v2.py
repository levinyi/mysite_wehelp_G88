import os
import sys
import pandas as pd
import re
import argparse
from multiDict import addtwodimdict


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


def read_maf(maf_file):
    maf_df = pd.read_csv(maf_file, sep="\t")
    # maf_df['cDNA_Changes'] = maf_df['cDNA_Changes'].fillna('')  # 将缺失值（NaN）替换为一个空字符串（''）
    # maf_df['homo/het'] = maf_df['homo/het'].fillna('')  # 将缺失值（NaN）替换为一个空字符串（''）
    return maf_df


def process_maf2dict(maf_df):
    # 将Gene，cDNA_Change，homo/het， 最终处理成如下这样的dict,每个基因中有homo，和het的列表。
    '''
        {
            "ABO": { 
                homo: ['c.826G>A','c.768C>T','c.678G>A'],
                het: [],
                }
            "ABCC1": {
                "homo": [],
                "het": []
            }
        }
    '''
    # 
    maf_dict = {}
    for index, row in maf_df.iterrows():
        gene = row['Gene']
        cDNA_Change = row['cDNA_Change']
        homo_het = row['homo/het']
        if homo_het == 'homo':
            addtwodimdict(maf_dict, gene, 'homo', cDNA_Change)
        else:
            addtwodimdict(maf_dict, gene, 'homo', cDNA_Change)
            addtwodimdict(maf_dict, gene, 'het', cDNA_Change)

    return maf_dict


def process_excel2df(excel_file):
    df = pd.read_excel(excel_file)
    # remove space in Nucleotide_change,
    return df


def core_process_of_identification_old(excel_df, gene, list1, list2):
    # Select rows where Gene column equals the given gene
    gene_df = excel_df[excel_df['Gene'] == gene]

    # Create a dictionary to store Nucleotide_change values for each allele name
    allele_nucleotide_changes = {}

    # Iterate over rows and populate the dictionary
    for _, row in gene_df.iterrows():
        allele_name = row['Allele_name(Het)']
        nucleotide_changes = str(row['Nucleotide_change']).split(';')
        
        if allele_name not in allele_nucleotide_changes:
            allele_nucleotide_changes[allele_name] = []
        
        allele_nucleotide_changes[allele_name].extend(nucleotide_changes)
    print(allele_nucleotide_changes)
    # {'ABO*A1.02': ['c.467C>T'], 'ABO*A2.01': ['c.467C>T', 'c.1061delC']}

    # Step 1
    matched_allele_homo = []
    matched_allele_het = []
    for allele_name, nucleotide_changes_list in allele_nucleotide_changes.items():
        if set(list1) == set(nucleotide_changes_list):
            matched_allele_homo.append(allele_name)

    if len(matched_allele_homo) >0:
        '''说明已经有最优的了，开始处理het'''
        for allele_name, nucleotide_changes_list in allele_nucleotide_changes.items():
            if set(list2) == set(nucleotide_changes_list):
                matched_allele_het.append(allele_name)
    
        if len(matched_allele_het) >0:
            '''说明het也有最优解了'''
            return matched_allele_homo, matched_allele_het
    else:
        '''说明没找到最优解，需要看list1的子集是不是全部包含nucleotide_changes_list'''

        for allele_name, nucleotide_changes_list in allele_nucleotide_changes.items():
            set1 = set(list1)
            set2 = set(nucleotide_changes_list)
            is_subset = set2.issubset(set1)
            if is_subset:
                # 获取交集
                intersection = set1.intersection(set2)

                # 获取非交集
                difference = set1.difference(set2)
                
                # 把非交集放到 list2 中
                list2.append(difference)
                
        # 找到最多的，返回
        


    return matched_allele_homo, matched_allele_het

def identify_matched_alleles(excel_df, gene, homo_nucleotides, het_nucleotides):
    # Filter the dataframe for the specified gene
    gene_specific_df = excel_df[excel_df['Gene'] == gene]

    # Store nucleotide changes for each allele
    allele_to_nucleotides = {}
    for _, row in gene_specific_df.iterrows():
        allele = row['Allele_name(Het)']
        changes = str(row['Nucleotide_change']).split(';')
        allele_to_nucleotides.setdefault(allele, []).extend(changes)
    # print(allele_to_nucleotides)
    
    # Identify matched alleles
    matched_alleles_homo, matched_alleles_het = [], []
    for allele, nucleotides in allele_to_nucleotides.items():
        nucleotides_set = set(nucleotides)

        if set(homo_nucleotides) == nucleotides_set:
            matched_alleles_homo.append(allele)

        if set(het_nucleotides) == nucleotides_set:
            matched_alleles_het.append(allele)
    
    print(matched_alleles_homo, matched_alleles_het)

    # Handle cases where no exact matches are found
    if not matched_alleles_homo:
        for allele, nucleotides in allele_to_nucleotides.items():
            # print(f"comparing {sorted(homo_nucleotides)} vs {sorted(nucleotides)}")
            if set(nucleotides).issubset(set(homo_nucleotides)):
                difference = set(homo_nucleotides).difference(nucleotides)
                het_nucleotides.extend(difference)

        # Re-run the matching process for het alleles after updating het_nucleotides
        for allele, nucleotides in allele_to_nucleotides.items():
            if set(het_nucleotides) == set(nucleotides):
                matched_alleles_het.append(allele)
    print(matched_alleles_homo, matched_alleles_het)

    return matched_alleles_homo, matched_alleles_het

def main(excel_file, hpa_db, maf_file, final_output_file):
    # deal with MAF file.
    maf_df = read_maf(maf_file)
    print(maf_df)

    maf_dict = process_maf2dict(maf_df)
    print(maf_dict)

    # 处理excel 数据
    excel_df = process_excel2df(excel_file)
    print(excel_df)
    # 按homo，het去处理。
    gene = 'ABO'
    list1 = maf_dict['ABO'].get('homo', [])
    list2 = maf_dict['ABO'].get('het',  [])
    print(gene, list1, list2)
    identify_matched_alleles(excel_df, gene, list1, list2)

    # for gene, cdna_dict in maf_dict.items():
    #     list1 = cdna_dict.get('homo',[])
    #     list2 = cdna_dict.get('het', [])
    #     matched_allele_homo, matched_allele_het = identify_matched_alleles(excel_df, gene, list1, list2)
    #     print(f'Gene: {gene}, Matched Allele (Homo): {matched_allele_homo}, Matched Allele (Het): {matched_allele_het}')


    '''
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
        
        # print("Now let's start with ", gene, cDNA_Changes_list)
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
    # output
    maf_df[['Number', 'System', 'Gene', 'PanelDesigned', 'MutationDetected',
        'cDNA_Changes','Variant_Classification', 'homo/het', 'Allele_name1', 'Allele_name2']].to_excel(final_output_file, index=False)
    maf_df[['Number', 'System', 'Gene', 'PanelDesigned', 'MutationDetected',
        'cDNA_Changes','Variant_Classification', 'homo/het', 'Allele_name1', 'Allele_name2']].to_csv(final_output_file.replace("xlsx", "csv"), index=False)

    print(maf_df)
    '''

if __name__ == '__main__':
    parser = _argparse()

    excel_file = parser.excel_file   # input: Blood.Gene.metadata.xlxs6
    hpa_db = parser.hpa_database
    maf_file   = parser.maf_file    # input: maf file
    final_output_file = parser.output_file  # output: final ouput table file

    if len(sys.argv) == 1:
        print(parser.help())
        sys.exit(1)

    main(excel_file, hpa_db, maf_file, final_output_file)

