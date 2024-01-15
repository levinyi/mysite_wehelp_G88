import subprocess
import pandas as pd
import os
import sys
sys.path.append(os.path.realpath(__file__))
import multiDict


def func(listTemp, n):
    '''Group the elements in the list in pairs'''
    for i in range(0, len(listTemp), n):
        yield listTemp[i:i + n]

result_path = sys.argv[1]
output_sites = sys.argv[2]

command = ["find", result_path, "-name", "HLA-HD_Result_final.result.txt"]
proc = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
out, err = proc.communicate()
file_list = [i.rstrip("\n") for i in out.decode().split('\n') if i]

print("find {} samples in current directory".format(len(file_list)))

# 初始化空 DataFrame
df = pd.DataFrame()

for each in file_list:
    sample_name = os.path.dirname(each).split("/")[-3] # /xxx/test-autoupload/S25002365_L01_P135/HLA-HD_Resule/result/HLA-HD_Result_final.result.txt
    sample_dict = {}
    with open(each,"r") as f:
        for line in f:
            line = line.rstrip("\n")
            a = line.split("\t")
            gene = a[0]
            for each in func(a[1:], 2):
                type1 = each[0]
                type2 = each[1]
                multiDict.addtwodimdict(sample_dict,gene,"type1",type1)
                multiDict.addtwodimdict(sample_dict,gene,"type2",type2)
    for gene, values in sample_dict.items():
        new_row = pd.DataFrame([{"sample_name": sample_name, "gene": gene, 
                                "type1": ";".join(values['type1']), 
                                "type2": ";".join(values['type2'])}])
        df = pd.concat([df, new_row], ignore_index=True)

# 创建一个空字典来存储每个样本的所有基因信息
sample_dict = {}

for index, row in df.iterrows():
    sample_name = row['sample_name']
    gene = row['gene']
    type1 = row['type1']
    type2 = row['type2']

    if sample_name not in sample_dict:
        sample_dict[sample_name] = {}

    # 将 'type1' 和 'type2' 信息添加为独立的列
    sample_dict[sample_name][f'{gene}-type1'] = type1
    sample_dict[sample_name][f'{gene}-type2'] = type2

# 现在将 sample_dict 转换为 DataFrame
new_df = pd.DataFrame.from_dict(sample_dict, orient='index')

# 重新命名列，如果需要
new_df.columns = ['HLA-' + col for col in new_df.columns]

if output_sites != "all_sites":
    sites = ['A', 'B','C','DRB1','DQB1','DPB1','DQA1','DPA1','DRB3','DRB4','DRB5']

    # 构建一个新的列名列表，包含每个基因的 type1 和 type2
    columns_to_keep = []
    for site in sites:
        columns_to_keep.append(f'HLA-{site}-type1')
        columns_to_keep.append(f'HLA-{site}-type2')

    # 筛选出新 DataFrame 中的相关列
    new_df = new_df[columns_to_keep ]

# 将index转换为列，叫做 'sample_name'
new_df.reset_index(inplace=True)
new_df.rename(columns={'index': 'Sample'}, inplace=True)

output_file = os.path.join(result_path, "HLA-HD_Result_final.summary.xls")
new_df.to_csv(output_file, sep="\t", index=False)

print("Output file: HLA-HD_Result_final.summary.xls")
print("Summary HLA-HD result done")

