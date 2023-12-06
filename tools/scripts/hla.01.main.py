import subprocess
import sys
import os
import fnmatch
import concurrent.futures
from commonFunction import deal_fastqc, find_files_by_suffix, process_fastq_files

BASE_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
software_path = os.path.join(BASE_DIR, "pipeline/software")
os.environ["PATH"] += os.pathsep + f"{software_path}/hlahd/hlahd.1.7.0/bin"
print(os.environ["PATH"])


def trim_fastq(sample_name, sample_file, project_dir, software_path):
    print("I'm trim_fastq function!")
    seqkit = os.path.join(software_path, "seqkit/seqkit")
    trim_threads = os.cpu_count()
    trim_portion = 0.1
    seqkit_seed = 120
    trim_string = "1M"
    trim_number = int(float(trim_string.rstrip("M"))*1000000)  # 1000000
    
    fq1_path = sample_file['fq1']
    fq2_path = sample_file['fq2']
    
    trim_data_dir = os.path.join(project_dir, sample_name, f"downsample_{trim_string}")

    os.makedirs(trim_data_dir, exist_ok=True)
    subprocess.run(f"{seqkit} sample -j {trim_threads} -p {trim_portion} -s {seqkit_seed} {fq1_path} | {seqkit} head -n {trim_number} -o {trim_data_dir}/{sample_name}_1.fq.gz --quiet", shell=True)
    subprocess.run(f"{seqkit} sample -j {trim_threads} -p {trim_portion} -s {seqkit_seed} {fq2_path} | {seqkit} head -n {trim_number} -o {trim_data_dir}/{sample_name}_2.fq.gz --quiet", shell=True)


def analyze_sample(sample_name, sample_files, project_dir, software_path, database_path, script_path, ref_fa):
    print("I'm analyze_sample function!")
    fq1 = sample_files['fq1']
    fq2 = sample_files['fq2']
    threads = os.cpu_count()

    sample_dir = os.path.join(project_dir, sample_name, 'mapping')
    os.makedirs(sample_dir, exist_ok=True)
    # print("sample_name: ", sample_name)
    
    # bwa
    bwa_command = (
        f"bwa mem -t {threads} -Y -H '@HD\\tVN:1.5\\tGO:none\\tSO:coordinate' "
        f"-R '@RG\\tID:TEST\\tSM:TEST\\tLB:Target\\tPL:BGI\\tPU:HVW2MCCXX:6:none' "
        f"{ref_fa} {fq1} {fq2} | "
        f"samtools view -buhS -t {ref_fa}.fai - | "
        f"samtools sort -o {sample_dir}/{sample_name}.hla.sortedByCoord.bam -"
    )
    subprocess.run(bwa_command, shell=True)

    subprocess.run(["samtools", "index", f"{sample_dir}/{sample_name}.hla.sortedByCoord.bam"])
    
    subprocess.run(f'samtools mpileup -f {ref_fa} {sample_dir}/{sample_name}.hla.sortedByCoord.bam -aa -A -Q 13 -o {sample_dir}/{sample_name}.pileup.Q13.out', shell=True)
    subprocess.run(f'python {script_path}/03_mileup_stats.py {sample_dir}/{sample_name}.pileup.Q13.out > {sample_dir}/{sample_name}.pileup.Q13.out.stats', shell=True)
    
    # extract fastq
    subprocess.run(f"samtools view -@ {threads} -bF 4 {sample_dir}/{sample_name}.hla.sortedByCoord.bam > {sample_dir}/{sample_name}.mapped.hla.bam", shell=True)
    subprocess.run(f"samtools view -@ 8 -bf 4 {sample_dir}/{sample_name}.hla.sortedByCoord.bam > {sample_dir}/{sample_name}.Unmapped.hla.bam", shell=True)
    subprocess.run(f"samtools fastq -@ {threads} -1 {sample_dir}/{sample_name}.mapped.hla.1.fastq -2 {sample_dir}/{sample_name}.mapped.hla.2.fastq -n {sample_dir}/{sample_name}.mapped.hla.bam", shell=True)
    
    # hlahd
    hlahd_command = (
        f"hlahd.sh -t {threads} -m 100 -c 0.95 "
        f"-f {software_path}/hlahd/hlahd.1.7.0/freq_data "
        f"{sample_dir}/{sample_name}.mapped.hla.1.fastq "
        f"{sample_dir}/{sample_name}.mapped.hla.2.fastq "
        f"{software_path}/hlahd/hlahd.1.7.0/HLA_gene.split.txt "
        f"{software_path}/hlahd/hlahd.1.7.0/dictionary/ "
        "HLA-HD_Result "
        f"{project_dir}/{sample_name}/"
    )
    subprocess.run(hlahd_command, shell=True)
    
    subprocess.run(f"python {script_path}/hla.02.freq.py {database_path}/hla.cwd.xls {project_dir}/{sample_name}/HLA-HD_Result/result/HLA-HD_Result_final.result.txt", shell=True)
    
    return  f"{project_dir}/{sample_name}/HLA-HD_Result/result/HLA-HD_Result_final.result.txt"


def main(data_dir, project_dir, software_path, database_path, script_path, ref_fa):
    fastq_list  = find_files_by_suffix(".fq.gz", data_dir)
    sample_dict = process_fastq_files(fastq_list)
    print("sample_dict: ", sample_dict)
    
    # Step 1 对原始数据进行 fastqc， 不需要等待执行结果。下面的fastqc会覆盖掉这个，想个办法解决，换个名字？
    with concurrent.futures.ProcessPoolExecutor() as executor:
        for sample_name, sample_files in sample_dict.items():
            executor.submit(deal_fastqc, sample_name, sample_files, project_dir, software_path)

    ### Step 2 对原始数据进行downsample. 要等待执行结果。
    ##################### 处理downsample文件夹下的fastq文件 #####################################
    with concurrent.futures.ProcessPoolExecutor() as executor:
        futures = [executor.submit(trim_fastq, sample_name, sample_files, project_dir, software_path) 
                   for sample_name, sample_files in sample_dict.items()]

        for future in concurrent.futures.as_completed(futures):
            future.result()

    ### Step 3 对downsample的数据进行fastqc， 不需要等待执行结果。
    fastq_list  = find_files_by_suffix(".fq.gz", project_dir)
    sample_dict = process_fastq_files(fastq_list)
    print("downsample_dict: ", sample_dict)

    # fastqc 
    with concurrent.futures.ProcessPoolExecutor() as executor:
        for sample_name, sample_files in sample_dict.items():
            executor.submit(deal_fastqc, sample_name, sample_files, project_dir, software_path, redirct=True)

    ### 对downsample数据进行 bwa + hlahd， 要等待执行结果。
    result_list = []
    with concurrent.futures.ProcessPoolExecutor() as executor:
        futures = [executor.submit(analyze_sample, sample_name, sample_files, project_dir, software_path, database_path, script_path, ref_fa) 
                   for sample_name, sample_files in sample_dict.items()]

        for future in concurrent.futures.as_completed(futures):
            result_list.extend(future.result())

    subprocess.run(f"python {script_path}/hla.03.summary.hla-hd.py {project_dir}",shell=True)

    ###########################################################################################
    # multiqc
    # multiqc was installed through pip install.
    subprocess.run(f"multiqc  . --outdir {project_dir} --quiet", shell=True)

    # write file list that need to be packaged into a file.
    with open(os.path.join(project_dir, "need_to_be_packaged.txt"), "w") as f:
        f.write(f"{project_dir}/HLA-HD_Result_final.summary.xls\n")

        # 添加其他要打包的文件， 需要确认！
        # for result in result_list:
        #     f.write(f"{result}\n")
    print("finished!")


if __name__ == "__main__":
    data_dir     = sys.argv[1]
    project_name = sys.argv[2]  # test1
    project_dir  = sys.argv[3]  # /data/storeData/ztron/rawdata/HLA/analysis/test3
    os.makedirs(project_dir, exist_ok=True)

    BASE_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    software_path = os.path.join(BASE_DIR, "pipeline/software")
    database_path = os.path.join(BASE_DIR, "pipeline/database")
    ref_fa = os.path.join(BASE_DIR, "pipeline/ref/hla/hla_gen.fasta")
    script_path = os.path.join(BASE_DIR, "tools/scripts")
    print(f"BASE_DIR: {BASE_DIR}")
    print(f"software_path: {software_path}")
    print(f"database_path: {database_path}")
    print(f"ref_fa: {ref_fa}")
    print(f"script_path: {script_path}")
    print(f"project_dir: {project_dir}")

    main(data_dir, project_dir, software_path, database_path, script_path, ref_fa)
