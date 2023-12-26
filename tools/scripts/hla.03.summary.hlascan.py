import os
import sys
import pandas as pd
from commonFunction import find_files_by_suffix


result_path = sys.argv[1]

file_list = find_files_by_suffix(".HLAscan.results.txt", result_path)
print("find {} samples in current directory".format(len(file_list)))

big_df = pd.DataFrame()
for each in file_list:
    print(each)
    if os.path.getsize(each) == 0:
        continue
    df = pd.read_csv(each, sep="\t")
    big_df = big_df._append(df)

# 转置big_df
big_df = big_df.T
print(big_df)
# 重置索引
big_df = big_df.reset_index()
# 重命名列名
print(big_df)
big_df.to_csv(os.path.join(result_path, "HLAscan_Result_final.summary.xls"), index=False, header=True)

print("Output file: HLAscan_Result_final.summary.xls")
print("Summary HLAscan result done")

