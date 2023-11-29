import pandas as pd
import os
import sys
from collections import defaultdict
import argparse

def deal_hpa_db(hpa_db):
    df = pd.read_csv(hpa_db, sep="\t", skiprows=1, names=['System', 'Gene', 'cDNA_Changes', 'homo', 'het'])
    db_dict = df.set_index('cDNA_Changes')[['homo', 'het']].T.to_dict()
    return db_dict

def main(project_dir, hpa_db, outputfile):
    # script_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    # hpa_db = os.path.join(script_dir, 'database',"HPA.Gene.cDNA_Changes.xls")

    hpa_db_dict = deal_hpa_db(hpa_db)
    file_list = [i.rstrip("\n") for i in os.popen(f"find {project_dir} -type f -name *brief.table.xls").readlines()]
    print(f"find {len(file_list)} samples in {project_dir} directory")

    header= ["Sample"] + [f"HPA-{i}"  for i in range(1, 36)]
    big_list = [header]

    for each in file_list:
        content = []
        sample_name = os.path.basename(os.path.dirname(os.path.dirname(each)))
        content.append(sample_name)

        df = pd.read_csv(each)
        df_hpa = df[df["System"] == "HPA"]
        # what if df_hpa is empty?
        
        cdna_changes = df_hpa['cDNA_Changes'].values.tolist()

        for cdna_change, homo_het in hpa_db_dict.items():
            if cdna_change in cdna_changes:
                df_filtered = df_hpa[(df_hpa['cDNA_Changes'] == cdna_change)]
                if (df_filtered['homo/het'] == 'homo').any():
                    content.append(homo_het['homo'])
                elif (df_filtered['homo/het'] == 'het').any():
                    content.append(homo_het['het'])
            else:
                content.append(f"{homo_het['homo']},{homo_het['het']}")

        big_list.append(content)
    big_df = pd.DataFrame(big_list[1:], columns=big_list[0])
    # print(big_df)
    big_df.to_csv(outputfile, index=False, header=True)
    print("all done!")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    # 定义参数 -i -d -o 
    parser.add_argument('-i', dest='project_dir', action='store', help="project dir")
    parser.add_argument('-d', dest='hpa_db',     action='store', help="hpa database file")
    parser.add_argument('-o', dest='outputfile', action="store", help="output file full path")

    # 解析命令行参数
    args = parser.parse_args()

    # 使用参数
    project_dir = args.project_dir
    hpa_db = args.hpa_db   ## "/path/to/database/HPA.Gene.cDNA_Changes.xls"
    outputfile = args.outputfile  #  "/path/to/HPA.summary.report.xls")

    main(project_dir, hpa_db, outputfile)
