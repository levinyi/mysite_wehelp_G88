import pandas as pd
import argparse
import re


def _argparse():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('-m', dest='maf_file', action="store", help="input : maf file")
    parser.add_argument('-s', dest='dbSNP_db', action="store", help="input: dbsnp_file.xls")
    parser.add_argument('-t', dest='template_file', action="store", help="input: template_file.xls")
    parser.add_argument('-o', dest="output_file",action="store", help="output: final output table file")
    return parser.parse_args()

def fill_hgvs(row):
    ''' 根据条件填充 hgvs 列的值 '''
    special_site = {
        'ABO': 3,
        'BSG': 627,
        'CD44': 0,
        'CD55': 6,
        'CR1': -1350,
        # 'FUT2': -33,
        'PIGG': 24,
        'SLC14A1': 396,
        'SLC44A2': 4
    }
    # 如果 hgvs 列不为空，则直接返回原值
    if pd.notna(row['hgvs']):
        return row['hgvs']
    # 如果 MutationDetected 列的值为 "NO"，则返回原值（也就是 np.nan）
    elif row['MutationDetected'] == 'NO':
        return row['hgvs']
    # 如果 MutationDetected 列的值为 "YES"，则按照规则填充 hgvs 列的值
    else:
        # 如果 Gene 在 special_site 字典中
        if row['Gene'] in special_site:
            # 返回 cDNA_Change 值 加上 special_site 字典中的值
            # 例如：c.6178A>T + 3 = c.6181A>T
            # 先使用re将数字提取出来，然后再加上 special_site 字典中的值
            num = int(re.findall(r'\d+', row['cDNA_Change'])[0]) + special_site[row['Gene']]
            return re.sub(r'\d+', str(num), row['cDNA_Change'])
        else:
            # 返回 cDNA_Change 的值
            return row['cDNA_Change']

def deal_dbsnp(dbsnp_df):
    ''' 处理dbSNP database，返回一个小的dataframe '''
    # 第28列是转录本NM号和hgvs信息：
    dbsnp_df[['transcript', 'hgvs']] = dbsnp_df.iloc[:, 28].str.split(':', n=1, expand=True)
    # 去掉版本号：.2/.3等
    dbsnp_df['transcript'] = dbsnp_df['transcript'].str[:-2]

    # 使用rename函数同时重命名多个列，50列是rs_number， 46列是allele_site。
    dbsnp_df.rename(columns={50: 'rs_number', 46: 'allele_site'}, inplace=True)
    dbsnp_df['rs_number'] = 'rs' + dbsnp_df['rs_number'].astype(str) 
    dbsnp_df['last3bp'] = dbsnp_df['hgvs'].str[-3:]  # 可能有bug，尤其是多个碱基的插入或缺失时，需要注意。
    
    # select columns
    dbsnp_df = dbsnp_df[['transcript', 'rs_number', 'hgvs', 'last3bp']]
    # drop duplicates
    dbsnp_df = dbsnp_df.drop_duplicates(subset=['transcript', 'rs_number', 'last3bp'])
    return dbsnp_df


def main():
    parser = _argparse()

    maf_file   = parser.maf_file    # input: maf file
    dbSNP_file = parser.dbSNP_db    # dbsnp file
    template_file = parser.template_file 
    out_file = parser.output_file  # output: final ouput table file

    ########################### dbSNP file ###########################################
    # 读取dbsnp_df时，直接指定需要转换为字符串的列
    dbsnp_df = pd.read_table(dbSNP_file, sep="\t", header=None, dtype={9: str, 18: str, 27: str, 50: str, 63: str, 65: str})
    dbsnp_df = deal_dbsnp(dbsnp_df)
    

    # print(dbsnp_df)
    ########################## template file ###################################
    template_df = pd.read_table(template_file, sep="\t")
    template_df['Number'] = template_df['Number'].apply(lambda x: str(x).zfill(3) if x !='/' else x)
    # print("template_df", template_df)
    ########################## MAF file ##################################
    maf_df = pd.read_table(maf_file, sep="\t", dtype={'Start_Position': str, 'End_Position': str})
    # for special site,这里只是将纯合的位点删掉。
    maf_df = maf_df[~((maf_df['Hugo_Symbol'] == 'CR1') & (maf_df['cDNA_Change'] == 'c.6178A>T') & (maf_df['AF'] == 1))]
    maf_df = maf_df[~((maf_df['Hugo_Symbol'] == 'GCNT2') & (maf_df['cDNA_Change'] == 'c.816C>G') & (maf_df['AF'] == 1))]
    merged_df = pd.merge(template_df, maf_df, left_on="Gene", right_on='Hugo_Symbol', how='left')
    # only print rows that Gene is 'KEL' 
    merged_df['MutationDetected'] = merged_df['Hugo_Symbol'].apply(lambda x : 'NO' if pd.isna(x) else 'YES')
    merged_df['homo/het'] = merged_df['AF'].apply(lambda x: "het" if x == 0.5 else 'homo' if x == 1 else float('nan'))
    merged_df['Transcript_stripped'] = merged_df['Transcript_USED'].str[:-2]
    merged_df['last3bp'] = merged_df['cDNA_Change'].str[-3:]
    # print(merged_df)
    # print(dbsnp_df)
    # select transcript== 'NM_000420' and rs_number== 'rs2293266' and  last3bp == 'T>C' in dbsnp_df
    # print(dbsnp_df[(dbsnp_df['transcript'] == 'NM_000420') & (dbsnp_df['rs_number'] == 'rs2293266')])

    # bug here, 合并时，如果组合中没有同时满足这三个，就会出现重复行？
    merged_df = pd.merge(merged_df, dbsnp_df, left_on=['Transcript_stripped','dbSNP_ID','last3bp'], right_on=['transcript', 'rs_number','last3bp'], how='left')
    # print(merged_df[merged_df['Gene'] == 'KEL'])
    #########################################################
    # for empty rs or empty hgvs
    ## for debug.
    merged_df['hgvs'] = merged_df.apply(fill_hgvs, axis=1)
    # print(merged_df)
    
    #########################################################
    ############ for special sites:
    ####### 1.FUT2 all site should offset 33bp.
    if 'FUT2' in merged_df['Gene'].values:
        merged_df.loc[merged_df['Gene'] == 'FUT2', 'hgvs'] = merged_df.loc[merged_df['Gene'] == 'FUT2', 'hgvs'].str.replace(r'c\.(\d+)(.*)', lambda m: f"c.{int(m.group(1)) - 33}{m.group(2)}", regex=True)
    
    ####### 2.ABO c.260_262insG mute. otherwise, c.261delG
    # 根据条件进行操作
    ABO_df = merged_df[(merged_df['Hugo_Symbol'] == 'ABO') & (merged_df['MutationDetected']=='YES')] # bug when 'No'?
    # print(ABO_df)
    if not ABO_df.empty:
        # 1) 当 hgvs 列没有 c.260_262insG 时，新增一行
        if 'c.260_262insG' not in ABO_df['hgvs'].values:
            new_row = pd.DataFrame({
                'Number':['001'], 'System':['ABO'], 'Gene':['ABO'], 'PanelDesigned':['YES'],'Transcript_USED': ['NM_020469.2'],'Hugo_Symbol': ['ABO'], 
                'Variant_Classification':'Splice_Site','cDNA_Change':['c.260_262insG'], 'MutationDetected' : ['YES'],
                'dbSNP_ID': ['rs8176719'], 'AF': [1], 'homo/het':['homo'],
                'rs_number': ['rs8176719'],'hgvs': ['c.261delG'],
            })
            merged_df = pd.concat([merged_df, new_row], ignore_index=True)
            # print(merged_df)
        # 2) 当 hgvs 列为 c.260_262insG 且 AF 等于 0.5 时，将 hgvs 的值改为 'c.261delG'
        merged_df.loc[(merged_df['Hugo_Symbol'] == 'ABO') & (merged_df['hgvs'] == 'c.260_262insG') & (merged_df['AF'] == 0.5), 'hgvs'] = 'c.261delG'
        # 3) 当 hgvs 列为 c.260_262insG 且 AF 等于 1 时，删除该行,因为纯合的插入G就跟ref一致，相当于没有突变。
        merged_df = merged_df[~((merged_df['Hugo_Symbol'] == 'ABO') & (merged_df['hgvs'] == 'c.260_262insG') & (merged_df['AF'] == 1))]
    
    ####### 3. CR  c.4828A>T, 如果发现该基因，则删除该行信息。以后补充纯合。
    merged_df.loc[(merged_df['Hugo_Symbol'] == 'CR1') & (merged_df['hgvs'] == 'c.4828A>T') & (merged_df['AF'] == 0.5), 'hgvs'] = 'c.4828T>A'

    ####### 4. GCNT2 c.816C>G，如果发现该基因，则删除该行信息。
    merged_df.loc[(merged_df['Hugo_Symbol'] == 'GCNT2') & (merged_df['hgvs'] == 'c.816C>G') & (merged_df['AF'] == 0.5), 'hgvs'] = 'c.816G>C'
    ###################################### DONE ###############################
    
    merged_df = merged_df.rename(columns={'hgvs':'cDNA_Changes'})
    
    ####################### output 
    merged_df.to_csv(out_file, sep="\t", index = False, header=True)
    merged_df[['Number','System', 'Gene', 'PanelDesigned', 'MutationDetected',
        'cDNA_Changes','Variant_Classification','homo/het']].to_csv(out_file+'.site.only.xls',sep="\t", index=False)
    # print(merged_df)

if __name__ == '__main__':
    main()