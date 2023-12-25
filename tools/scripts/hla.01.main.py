import subprocess
import sys
import os
import fnmatch
import concurrent.futures
from commonFunction import deal_fastqc, find_files_by_suffix, process_fastq_files

BASE_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
software_path = os.path.join(BASE_DIR, "pipeline/software")
os.environ["PATH"] += os.pathsep + f"{software_path}/hlahd/hlahd.1.7.0/bin"
# print(os.environ["PATH"])


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


def analyze_sample(sample_name, sample_files, project_dir, software_path, database_path, script_path, ref_fa, software_list):
    print("I'm analyze_sample function!")
    fq1 = sample_files['fq1']
    fq2 = sample_files['fq2']
    threads = os.cpu_count()

    sample_dir = os.path.join(project_dir, sample_name, 'mapping')
    os.makedirs(sample_dir, exist_ok=True)
    
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
    subprocess.run(f'python {script_path}/bwa_mileup_stats.py {sample_dir}/{sample_name}.pileup.Q13.out > {sample_dir}/{sample_name}.pileup.Q13.out.stats', shell=True)
    
    # extract fastq
    subprocess.run(f"samtools view  -@ {threads} -bF 4 {sample_dir}/{sample_name}.hla.sortedByCoord.bam > {sample_dir}/{sample_name}.mapped.hla.bam", shell=True)
    subprocess.run(f"samtools view  -@ {threads} -bf 4 {sample_dir}/{sample_name}.hla.sortedByCoord.bam > {sample_dir}/{sample_name}.Unmapped.hla.bam", shell=True)
    subprocess.run(f"samtools fastq -@ {threads} -1    {sample_dir}/{sample_name}.mapped.hla.1.fastq -2 {sample_dir}/{sample_name}.mapped.hla.2.fastq -n {sample_dir}/{sample_name}.mapped.hla.bam", shell=True)
    
    return_list = []
    if 'hlahd' in software_list:
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
        subprocess.run(f"python {script_path}/hla.02.freq.py {database_path}/hla.cwd.xls {project_dir}/{sample_name}/HLA-HD_Result/result/{sample_name}.HLA-HD_Result_final.result.txt", shell=True)
        return_list.append(f"{project_dir}/{sample_name}/HLA-HD_Result/result/{sample_name}.HLA-HD_Result_final.result.txt")
    if 'hla_scan' in software_list:
        # hla_scan
        print("Start hla_scan!")
        hla_scan_dir = os.path.join(project_dir, sample_name, "HLAscan_Result")
        os.makedirs(hla_scan_dir, exist_ok=True)
        hla_genes = ["HLA-A","HLA-B","HLA-C","HLA-DMA","HLA-DMB","HLA-DOA","HLA-DOB","HLA-DPA1","HLA-DPB1","HLA-DQA1","HLA-DQB1","HLA-DRA","HLA-DRB1","HLA-DRB5","HLA-E","HLA-F","HLA-G","MICA","MICB","TAP1","TAP2"]
        with open(f"{hla_scan_dir}/{sample_name}.HLAscan.shell.sh", "w") as f:
            for gene in hla_genes:
                f.write(f"{software_path}/hla_scan/hla_scan -t {threads} -l {sample_dir}/{sample_name}.mapped.hla.1.fastq -r {sample_dir}/{sample_name}.mapped.hla.2.fastq -d {software_path}/hla_scan/db/HLA-ALL.IMGT -g {gene} > {hla_scan_dir}/{sample_name}.{gene}.out.txt\n")
                subprocess.run(f"{software_path}/hla_scan/hla_scan -t {threads} -l {sample_dir}/{sample_name}.mapped.hla.1.fastq -r {sample_dir}/{sample_name}.mapped.hla.2.fastq -d {software_path}/hla_scan/db/HLA-ALL.IMGT -g {gene} > {hla_scan_dir}/{sample_name}.{gene}.out.txt", shell=True)
        subprocess.run(f"python3 {script_path}/hla.04.merge_HLAscan_result.py {hla_scan_dir}/{sample_name}*.out.txt > {hla_scan_dir}/{sample_name}.HLAscan.results.txt\n", shell=True)
        return_list.append(f"{hla_scan_dir}/{sample_name}.HLAscan.results.txt")
    
    if 'OptiType' in software_list:
        # optitype
        print("Start OptiType!")
        opt_dir = os.path.join(project_dir, sample_name ,"OptiType_Result")
        os.makedirs(opt_dir)
        OptiType = os.path.join(software_path, "OptiType-1.3.5", "OptiTypePipeline.py")
        razers3 = os.path.join(software_path, "razers", "razers3")
        OptiType_ref = os.path.join(ref_path, "hla_reference_dna.fasta")
        subprocess.run(f"{razers3} -tc {threads} -i 95 -m 1 -dr 0 -o {opt_dir}/sample.fished_1.bam {OptiType_ref} {fq1}\n", shell=True)
        subprocess.run(f"{razers3} -tc {threads} -i 95 -m 1 -dr 0 -o {opt_dir}/sample.fished_2.bam {OptiType_ref} {fq2}\n", shell=True)
        subprocess.run(f"samtools bam2fq {opt_dir}/sample.fished_1.bam > {opt_dir}/sample.fished_1.fastq\n", shell=True)
        subprocess.run(f"samtools bam2fq {opt_dir}/sample.fished_2.bam > {opt_dir}/sample.fished_2.fastq\n", shell=True)
        subprocess.run(f"python2.7 {OptiType} -i {opt_dir}/sample.fished_1.fastq {opt_dir}/sample.fished_2.fastq --dna -v -o {opt_dir} -p {sample_name}\n", shell=True)
        return_list.append(f"{opt_dir}/{sample_name}.result.tsv")
    
    return return_list


def main(data_dir, project_dir, software_path, database_path, script_path, ref_fa, software_list):
    fastq_list  = find_files_by_suffix(".fq.gz", data_dir)
    sample_dict = process_fastq_files(fastq_list)
    
    # Step 1 对原始数据进行 fastqc， 不需要等待执行结果。
    with concurrent.futures.ProcessPoolExecutor() as executor:
        for sample_name, sample_files in sample_dict.items():
            executor.submit(deal_fastqc, sample_name, sample_files, project_dir, software_path, redirct=True)

    '''
    if 'hlahd' in software_list:
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
        with concurrent.futures.ProcessPoolExecutor() as executor:
            for sample_name, sample_files in sample_dict.items():
                executor.submit(deal_fastqc, sample_name, sample_files, project_dir, software_path)
    else:
        print("No need to downsample!")
    
    '''
    ### Step 4 multiqc  # multiqc was installed through pip install.
    subprocess.run(f"multiqc {project_dir}  --outdir {project_dir}", shell=True)

    ### Step 5 对downsample数据进行 bwa + hlahd， 要等待执行结果。
    result_list = []
    with concurrent.futures.ProcessPoolExecutor() as executor:
        futures = [executor.submit(analyze_sample, sample_name, sample_files, project_dir, software_path, database_path, script_path, ref_fa, software_list) 
                   for sample_name, sample_files in sample_dict.items()]

        for future in concurrent.futures.as_completed(futures):
            result_list.extend(future.result())
    
    print(f"result_list: {result_list}")
    ###########################################################################################
    # write file list that need to be packaged into a file.
    with open(os.path.join(project_dir, "need_to_be_packaged.txt"), "w") as f:
        if 'hlahd' in software_list:
            # 额外统计HLA-HD的结果，并写入文件
            subprocess.run(f"python {script_path}/hla.03.summary.hla-hd.py {project_dir}",shell=True)
            f.write(f"{project_dir}/HLA-HD_Result_final.summary.xls\n")

        # 添加其他要打包的文件， 需要确认！
        for result in result_list:
            print(f"result: {result}")
            f.write(f"{result}\n")
    print("Package result finished!")


if __name__ == "__main__":
    data_dir     = sys.argv[1]
    project_name = sys.argv[2]  # test1
    project_dir  = sys.argv[3]  # /data/storeData/ztron/rawdata/HLA/analysis/test3
    software_list = sys.argv[4].split(",")  # ['hla-hd', 'hla_scan', 'optitype']
    print(f"software_list in hla.main.py 1st : {software_list}")
    os.makedirs(project_dir, exist_ok=True)

    BASE_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    software_path = os.path.join(BASE_DIR, "pipeline/software")
    database_path = os.path.join(BASE_DIR, "pipeline/database")
    ref_path = os.path.join(BASE_DIR, "pipeline/ref/hla")
    ref_fa = os.path.join(ref_path, "hla_gen.fasta")
    script_path = os.path.join(BASE_DIR, "tools/scripts")

    main(data_dir, project_dir, software_path, database_path, script_path, ref_fa, software_list)
