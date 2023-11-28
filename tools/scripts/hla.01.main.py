import subprocess
import sys
import os
import fnmatch
from concurrent.futures import ProcessPoolExecutor

os.environ["PATH"] += os.pathsep + "/data/webapp/mysite/pipeline/software/hlahd/hlahd.1.7.0/bin"
print(os.environ["PATH"])

def find_fastq_files(directory):
    fastq_files = []
    for root, dirs, files in os.walk(directory):
        for filename in files:
            if fnmatch.fnmatch(filename, '*.fq.gz'):
                fastq_files.append(os.path.join(root, filename))
    return fastq_files

def find_files_by_suffix(suffix, path):
    result = []
    for root, dirs, files in os.walk(path):
        for file in files:
            if file.endswith(suffix):
                result.append(os.path.join(root, file))
    return result

def process_fastq_files(fastq_files_list):
    sample_dict = {}
    
    for filepath in fastq_files_list:
        filename = os.path.basename(filepath)  # 获取文件名，如'CX1680_Raw_1.fq.gz'
        sample_name = filename.split('.')[0][:-2]   # 提取样本名称，如'CX1680_Raw'
        file_number = filename.split('_')[-1]  # 提取文件编号，如'1.fq.gz'
        file_key = 'fq1' if file_number == '1.fq.gz' else 'fq2'  # 判断文件类型
        
        if sample_name in sample_dict:
            sample_dict[sample_name][file_key] = filepath
        else:
            sample_dict[sample_name] = {file_key: filepath}
    
    return sample_dict

def analyze_sample(sample_name, sample_files, project_dir, software_path, database_path, script_path, ref_fa):
    fq1 = sample_files['fq1']
    fq2 = sample_files['fq2']
    sample_dir = os.path.join(project_dir, sample_name)
    os.makedirs(sample_dir, exist_ok=True)
    print("sample_name: ", sample_name)

    # fastqc
    subprocess.run(["fastqc", "--threads", "8", "--format", "fastq", "--outdir", sample_dir, "--quiet", fq1, fq2])
    # bwa
    bwa_command = (
        f"bwa mem -t 8 -Y -H '@HD\\tVN:1.5\\tGO:none\\tSO:coordinate' "
        f"-R '@RG\\tID:TEST\\tSM:TEST\\tLB:Target\\tPL:BGI\\tPU:HVW2MCCXX:6:none' "
        f"{ref_fa} {fq1} {fq2} | "
        f"samtools view -buhS -t {ref_fa}.fai - | "
        f"samtools sort -o {sample_dir}/{sample_name}.hla.sortedByCoord.bam -"
    )
    subprocess.run(bwa_command, shell=True)

    subprocess.run(["samtools", "index", f"{sample_dir}/{sample_name}.hla.sortedByCoord.bam"])
    
    # extract fastq
    subprocess.run(["samtools", "view", "-@", "8", "-bF", "4", f"{sample_dir}/{sample_name}.hla.sortedByCoord.bam", ">", f"{sample_dir}/{sample_name}.mapped.hla.bam"])
    subprocess.run(["samtools", "view", "-@", "8", "-bf", "4", f"{sample_dir}/{sample_name}.hla.sortedByCoord.bam", ">", f"{sample_dir}/{sample_name}.Unmapped.hla.bam"])
    subprocess.run(["samtools", "fastq", "-@", "8", "-1", f"{sample_dir}/{sample_name}.mapped.hla.1.fastq",
                    "-2", f"{sample_dir}/{sample_name}.mapped.hla.2.fastq",
                    "-n", f"{sample_dir}/{sample_name}.mapped.hla.bam"])

    # hlahd
    hlahd_command = [
        "hlahd.sh", "-t", "8", "-m", "100", "-c", "0.95",
        "-f", f"{software_path}/hlahd/hlahd.1.7.0/freq_data",
        f"{sample_dir}/{sample_name}.mapped.hla.1.fastq",
        f"{sample_dir}/{sample_name}.mapped.hla.2.fastq",
        f"{software_path}/hlahd/hlahd.1.7.0/HLA_gene.split.txt",
        f"{software_path}/hlahd/hlahd.1.7.0/dictionary/",
        "HLA-HD_Result", f"{sample_dir}/"
    ]
    subprocess.run(hlahd_command)
    
    subprocess.run(["python3", f"{script_path}/hla.freq.py", f"{database_path}/hla.cwd.xls", f"{sample_dir}/HLA-HD_Result/result/HLA-HD_Result_final.result.txt"])

    return f"{sample_dir}/{sample_name}_1_fastqc.html", f"{sample_dir}/{sample_name}_2_fastqc.html"

def main(upload_dir, project_dir, software_path, database_path, script_path, ref_fa):
    # ## start analysis !!!!!
    fastq_list = find_files_by_suffix(".fq.gz", upload_dir)
    sample_dict = process_fastq_files(fastq_list)

    result_list = []
    with ProcessPoolExecutor(max_workers=8) as executor:
        futures = [executor.submit(analyze_sample, sample_name, sample_files, project_dir, software_path, database_path, script_path, ref_fa) 
                   for sample_name, sample_files in sample_dict.items()]

        for future in futures:
            result_list.extend(future.result())

    subprocess.run(["python", f"{script_path}/hla.summary.hla-hd.py", f"{project_dir}"])

    ###########################################################################################
    # write file list that need to be packaged into a file.
    with open(os.path.join(project_dir, "need_to_be_packaged.txt"), "w") as f:
        f.write(f"{project_dir}/HLA-HD_Result_final.summary.xls\n")
        for result in result_list:
            f.write(f"{result}\n")


if __name__ == "__main__":
    upload_dir = sys.argv[1] # /data/webapp/mysite/pipeline/project_hla/8833a5c2-0fb1-4c13-b734-6acace996f08
    project_name = sys.argv[2]  # test-autoupload

    project_dir = os.path.join(upload_dir , project_name)
    os.makedirs(project_dir, exist_ok=True)

    BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    software_path = os.path.join(BASE_DIR, "software")
    database_path = os.path.join(BASE_DIR, "database")
    script_path = os.path.join(BASE_DIR, "scripts")
    ref_fa = os.path.join(BASE_DIR, "ref/hla/hla_gen.fasta")

    main(upload_dir, project_dir, software_path, database_path, script_path, ref_fa)
