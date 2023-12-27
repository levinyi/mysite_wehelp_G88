import os
import sys
import pandas as pd
from commonFunction import find_files_by_suffix
# import ast  # Import the ast module for literal_eval


result_path = sys.argv[1]

file_list = find_files_by_suffix(".HLAscan.results.txt", result_path)
print("find {} samples in current directory".format(len(file_list)))

big_df = pd.DataFrame()
for each in file_list:
    if os.path.getsize(each) == 0:
        continue
    df = pd.read_csv(each, sep="\t")

    # Set 'HLA' column as the index
    df.set_index('HLA', inplace=True)

    # Transpose the DataFrame
    df = df.T

    # Convert each cell value to a semicolon-separated string
    df = df.map(lambda x: ';'.join(eval(x) if isinstance(x, str) else x) if pd.notna(x) else None)

    big_df = big_df._append(df)

# print(big_df)
# Rename the index header to 'Sample'
# set index to first column
big_df = big_df.reset_index()
big_df = big_df.rename(columns={'index': 'Sample'})
print(big_df)
big_df.to_csv(os.path.join(result_path, "HLAscan_Result_final.summary.csv"), index=False, header=True)

print("Output file: HLAscan_Result_final.summary.xls")
print("Summary HLAscan result done")

