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
command = ["find", result_path, "-name", "HLA-HD_Result_final.result.txt"]
proc = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
out, err = proc.communicate()
file_list = [i.rstrip("\n") for i in out.decode().split('\n') if i]

print("find {} samples in current directory".format(len(file_list)))
with open(os.path.join(result_path, "HLA-HD_Result_final.summary.xls"),"w") as output:
    for each in file_list:
        sample_name = os.path.dirname(each).split("/")[-3] # /xxx/test-autoupload/S25002365_L01_P135/HLA-HD_Resule/result/HLA-HD_Result_final.result.txt
             # /data/storeData/ztron/rawdata/HLA/analysis/test_hla_4/FT100020530_L01_UDB-P045/mapping/HLA-HD_Result/result/HLA-HD_Result_final.result.txt
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
        results = []
        for gene, values in sample_dict.items():
            results.append(";".join(values['type1']))
            results.append(";".join(values['type2']))
        output.write("{}\t{}\n".format(sample_name, "\t".join(results)))

print("output: HLA-HD_Result_final.summary.xls")
print("done")

