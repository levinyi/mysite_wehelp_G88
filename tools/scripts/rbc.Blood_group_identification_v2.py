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
                "homo": ['c.826G>A', 'c.768C>T', 'c.678G>A'],
                "het": ['c.106G>T', 'c.188G>A', 'c.189C>T', 'c.220C>T', 'c.261delG'],
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
        elif homo_het == 'het':
            # addtwodimdict(maf_dict, gene, 'homo', cDNA_Change)
            addtwodimdict(maf_dict, gene, 'het', cDNA_Change)

    return maf_dict


def process_excel2df(excel_file):
    df = pd.read_excel(excel_file)
    # remove space in Nucleotide_change,
    return df


def specified_gene(df, gene):
    # Filter the dataframe for the specified gene
    gene_specific_df = df[df['Gene'] == gene]

    # Store nucleotide changes for each allele
    allele_to_nucleotides = {}
    common_allele_dict = {}
    for _, row in gene_specific_df.iterrows():
        allele = row['Allele_name(Het)']
        changes = str(row['Nucleotide_change']).split(';')
        common = row['Allele_name(Homo)']
        allele_to_nucleotides.setdefault(allele, []).extend(changes)
        common_allele_dict[allele] = common # bug

    return allele_to_nucleotides, common_allele_dict


def identify_matched_alleles(excel_df, gene, homo_nucleotides, het_nucleotides):
    homo_nucleotides_set = set(homo_nucleotides)
    # print(f"I'm processing {gene} {homo_nucleotides}, {het_nucleotides}")
    allele_to_nucleotides, common_allele_dict = specified_gene(excel_df, gene)

    # Identify matched alleles in list 1:
    allele_paired = []
    for allele1, nucleotides in allele_to_nucleotides.items():
        nucleotides_set = set(nucleotides)

        if nucleotides_set.issubset(homo_nucleotides_set):
            # 把剩下的归到list2中，但不能把homo的再归到list2，就重复了，所以再加一步判断。
            # 也就是说这个homo的一定要在list1中被使用才行。
            rest_site = []
            for each in homo_nucleotides:
                if each not in nucleotides:
                    rest_site.append(each)

            if len(rest_site) != 0:
                het_nucleotides_set = het_nucleotides + rest_site
                if len(het_nucleotides) == len(set(het_nucleotides_set)):
                    # 表示合并后和合并前是相同的表，那么有可能是把homo的位点合并了，这是不允许的，返回error。
                    allele_paired.append((allele1, 'error'))
                else:
                    for allele2, nucleotides2 in allele_to_nucleotides.items():
                        nucleotides_set = set(nucleotides2)
                        if nucleotides_set == set(het_nucleotides_set):
                            # prefect match !
                            allele_paired.append((allele1, allele2))
                        else:
                            allele_paired.append((allele1, '-'))
            else:
                if len(het_nucleotides) == 0:
                    common_allele = common_allele_dict[allele1]
                    allele_paired.append((allele1, common_allele))
                else:
                    for allele2, nucleotides2 in allele_to_nucleotides.items():
                        nucleotides_set = set(nucleotides2)
                        if nucleotides_set == set(het_nucleotides):
                            # prefect match !
                            allele_paired.append((allele1, allele2))
                        else:
                            allele_paired.append((allele1, '-'))

    # 使用集合去除冗余
    unique_tuples = set(frozenset(item) for item in allele_paired)
    
    # 剔除掉包含有error的元组
    unique_tuples = [item for item in unique_tuples if 'error' not in item]

    # 将结果转换回列表
    result_list = [tuple(item) for item in unique_tuples]
    
    if len(result_list) == 0:
        result_list = [("-","-")]
    
    return result_list


def main(excel_file, hpa_db, maf_file, final_output_file):
    # deal with MAF file.
    maf_df = read_maf(maf_file)
    # print(maf_df)

    maf_dict = process_maf2dict(maf_df)
    # print(maf_dict)

    # 处理excel 数据
    excel_df = process_excel2df(excel_file)
    # print(excel_df)

    new_df = pd.DataFrame()
    for gene, cdna_dict in maf_dict.items():
        list1 = cdna_dict.get('homo',[]) + cdna_dict.get('het',[])
        list2 = cdna_dict.get('homo', [])
        print(gene, list1, list2)
        name_list = identify_matched_alleles(excel_df, gene, list1, list2)
        for allele1, allele2 in name_list:
            new_df = new_df._append({'Gene':gene, 'Allele_name1': allele1, 'Allele_name2': allele2}, ignore_index=True)
    
    # print("new_df")
    # print(new_df)
    ###############################################
    maf_df = pd.merge(maf_df, new_df, how='left', on='Gene')
    # print("maf_df")
    # print( maf_df)
    # maf_df['Allele_name1'] = maf_df['Allele_name1'].apply(lambda x: "".join(x))
    # maf_df['Allele_name2'] = maf_df['Allele_name2'].apply(lambda x: "".join(x))
    maf_df['Allele_name1'] = maf_df['Allele_name1'].apply(lambda x: "".join(x) if isinstance(x, list) else x)
    maf_df['Allele_name2'] = maf_df['Allele_name2'].apply(lambda x: "".join(x) if isinstance(x, list) else x)

    # 
    maf_df['Allele_name1'] = maf_df.groupby('Gene',group_keys=False)['Allele_name1'].apply(remove_duplicates)
    maf_df['Allele_name2'] = maf_df.groupby('Gene',group_keys=False)['Allele_name2'].apply(remove_duplicates)

    #################################################
    #  for HPA
    if 'HPA' in maf_df['System'].values:
        hpa_df = pd.read_csv(hpa_db, sep="\t")
        maf_df = pd.merge(maf_df,hpa_df, how="left", on=['Gene','System','cDNA_Changes'])
        maf_df['Allele_name1'] = maf_df['homo']
        maf_df['Allele_name2'] = maf_df['het']

    #################################################
    # output
    maf_df[['Number', 'System', 'Gene', 'PanelDesigned', 'MutationDetected',
        'cDNA_Changes','Variant_Classification', 'homo/het', 'Allele_name1', 'Allele_name2']].to_excel(final_output_file, index=False)
    maf_df[['Number', 'System', 'Gene', 'PanelDesigned', 'MutationDetected',
        'cDNA_Changes','Variant_Classification', 'homo/het', 'Allele_name1', 'Allele_name2']].to_csv(final_output_file.replace("xlsx", "csv"), index=False)

    print(maf_df)


if __name__ == '__main__':
    parser = _argparse()

    excel_file = parser.excel_file   # input: Blood.Gene.metadata.v2.xlsx
    hpa_db = parser.hpa_database
    maf_file   = parser.maf_file    # input: maf file
    final_output_file = parser.output_file  # output: final ouput table file

    if len(sys.argv) == 1:
        print(parser.help())
        sys.exit(1)

    main(excel_file, hpa_db, maf_file, final_output_file)

