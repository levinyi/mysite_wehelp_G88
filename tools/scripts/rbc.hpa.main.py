import os
import sys
import argparse
import subprocess
import concurrent.futures
from commonFunction import deal_fastqc, find_files_by_suffix, process_fastq_files


def run_bwa_mem(fq1, fq2, ref_fa, sample_name, project_dir):
    # BWA MEM command
    print(f"BWA MEM Start: {sample_name}")

    sample_folder = os.path.join(project_dir, sample_name)
    os.makedirs(sample_folder, exist_ok=True)

    threads = 8

    bwa_command = (
        f"bwa mem -t {threads} -Y -H '@HD\\tVN:1.5\\tGO:none\\tSO:coordinate' "
        f"-R '@RG\\tID:TEST\\tSM:TEST\\tLB:Target\\tPL:BGI\\tPU:HVW2MCCXX:6:none' "
        f"{ref_fa} {fq1} {fq2} | "
        f"samtools view -buhS -t {ref_fa}.fai - | "
        f"samtools sort -o {project_dir}/{sample_name}/{sample_name}.Rawsample.hg38.sortedByCoord.bam -"
    )
    subprocess.run(bwa_command, shell=True)
    subprocess.run([
        f"samtools index {project_dir}/{sample_name}/{sample_name}.Rawsample.hg38.sortedByCoord.bam"],
        shell=True)

def run_gatk(ref_fa, sample_name, project_dir):
    # GATK HaplotypeCaller command
    subprocess.run([
        f"{software_path}/gatk/gatk --java-options \"-Xmx50G\" HaplotypeCaller -R {ref_fa} -I {project_dir}/{sample_name}/{sample_name}.Rawsample.hg38.sortedByCoord.bam -L {database_path}/Blood.gene.bed -O {project_dir}/{sample_name}/{sample_name}.Rawsample.output.raw.vcf --QUIET true --verbosity ERROR --create-output-variant-index false"],
        shell=True)

    # GATK VariantFiltration command
    subprocess.run([
        f"{software_path}/gatk/gatk --java-options \"-Xmx50G\" VariantFiltration -R {ref_fa} -V {project_dir}/{sample_name}/{sample_name}.Rawsample.output.raw.vcf -filter \"QD < 2.0\"  --filter-name \"QDFilter\" -filter \"MQ < 40.0\" --filter-name \"MQFilter\" -filter \"FS > 60.0\" --filter-name \"FSFilter\" --create-output-variant-index false --output {project_dir}/{sample_name}/{sample_name}.Rawsample.output.filtered.vcf"],
        shell=True)

    # GATK Funcotator command
    subprocess.run([
        f"{software_path}/gatk/gatk --java-options \"-Xmx50G\" Funcotator --data-sources-path {database_path}/funcotator_dataSources.v1.7.20200521s --ref-version hg38 --output-file-format MAF --reference {ref_fa} --exclude-field Center --exclude-field Tumor_Sample_Barcode --exclude-field Matched_Norm_Sample_Barcode --exclude-field Match_Norm_Seq_Allele1 --exclude-field Match_Norm_Seq_Allele2 --exclude-field Tumor_Validation_Allele1 --exclude-field Tumor_Validation_Allele2 --exclude-field Match_Norm_Validation_Allele1 --exclude-field Match_Norm_Validation_Allele2 --exclude-field Verification_Status --exclude-field Validation_Status --exclude-field Mutation_Status --exclude-field Sequencing_Phase --exclude-field Sequence_Source --exclude-field Validation_Method --exclude-field Score --exclude-field BAM_File --exclude-field Sequencer --exclude-field Tumor_Sample_UUID --exclude-field Matched_Norm_Sample_UUID --variant {project_dir}/{sample_name}/{sample_name}.Rawsample.output.filtered.vcf --output {project_dir}/{sample_name}/{sample_name}.Rawsample.funcotated.raw.MAF --remove-filtered-variants true"],
        shell=True)

    # GREP remove MAF header
    subprocess.run([
        f"grep -v \"^#\" {project_dir}/{sample_name}/{sample_name}.Rawsample.funcotated.raw.MAF > {project_dir}/{sample_name}/{sample_name}.Rawsample.funcotated.MAF.xls"],
        shell=True)
    
def deal_maf(sample_name, project_dir, output_templates_file, script_path, database_path, hpa=False):
    # Python scripts filter MAF table
    if hpa:
        run_command1 = (
            f"python3 {script_path}/rbc.deal_maf_table.py  "
            f"-m {project_dir}/{sample_name}/{sample_name}.Rawsample.funcotated.MAF.xls "
            f"-o {project_dir}/{sample_name}/{sample_name}.Rawsample.funcotated.small.MAF.xls "
            f"--hpa "
        )
    else:
        run_command1 = (
            f"python3 {script_path}/rbc.deal_maf_table.py  "
            f"-m {project_dir}/{sample_name}/{sample_name}.Rawsample.funcotated.MAF.xls "
            f"-o {project_dir}/{sample_name}/{sample_name}.Rawsample.funcotated.small.MAF.xls "
        )
    subprocess.run(run_command1, shell=True)

    # add dbsnp HVGS info
    run_command2 = (
        f"python3 {script_path}/rbc.maf_add_dbSNP_HVGS.py "
        f"-m {project_dir}/{sample_name}/{sample_name}.Rawsample.funcotated.small.MAF.xls "
        f"-s {database_path}/dbSNP.db.txt -t {output_templates_file} "
        f"-o {project_dir}/{sample_name}/{sample_name}.Rawsample.funcotated.small.add.dbsnp.MAF.xls"
    )
    subprocess.run(run_command2, shell=True)


def blood_identy(sample_name, project_dir, database_path, script_path):
    # identify blood group info
    run_command = (
        f"python {script_path}/rbc.Blood_group_identification_v2.py "
        f"-b {database_path}/Blood.Gene.metadata.manually.checked.v2.xlsx "
        f"-p {database_path}/HPA.Gene.cDNA_Changes.xls "
        f"-m {project_dir}/{sample_name}/{sample_name}.Rawsample.funcotated.small.add.dbsnp.MAF.xls "
        f"-o {project_dir}/{sample_name}/{sample_name}.Rawsample.funcotated.brief.table.xlsx"
    )
    subprocess.run(run_command, shell=True)

    result_list = [f"{project_dir}/{sample_name}/{sample_name}.Rawsample.funcotated.brief.table.xlsx"]

    return result_list


def extract_cds(sample_name, project_dir, script_path, database_path):
    # extract CDS sequence
    run_command = (
        f"python {script_path}/rbc.extractCDS.py "
        f"{project_dir}/{sample_name}/{sample_name}.Rawsample.funcotated.small.add.dbsnp.MAF.xls "
        f"{database_path}/RBC.HPA.CDS.60.fasta "
        f"{project_dir}/{sample_name}/{sample_name}.Rawsample.CDS.fasta"
    ) 
    subprocess.run(run_command, shell=True)

    result_list = [f"{project_dir}/{sample_name}/{sample_name}.Rawsample.CDS.fasta"]
    
    return result_list


def main(data_dir, project_dir, software_path, database_path, script_path, ref_fa, output_templates_file, hpa):
    fastq_list  = find_files_by_suffix(".fq.gz", data_dir)
    # 如果没有fq.gz文件，就退出
    if len(fastq_list) == 0:
        print("No fq.gz files found!")
        sys.exit(1)

    sample_dict = process_fastq_files(fastq_list)
    # print("sample_dict: ", sample_dict)

    print("Start Step1 Fastqc!")
    # Step 1 对原始数据进行 fastqc， 不需要等待执行结果。
    fastqc_files = find_files_by_suffix("fastqc.html", project_dir)
    if len(fastqc_files) == 0:
        with concurrent.futures.ProcessPoolExecutor() as executor:
            for sample_name, sample_files in sample_dict.items():
                executor.submit(deal_fastqc, sample_name, sample_files, project_dir, software_path, redirct=True)
    else:
        print("\tfastqc files exist, skip fastqc!")

    multiqc = find_files_by_suffix("multiqc_report.html", project_dir)
    if len(multiqc) == 0:
        subprocess.run(f"multiqc {project_dir}  --outdir {project_dir}", shell=True)
    else:
        print("\tmultiqc files exist, skip multiqc!")
    
    ######################################################################
    print("Start Step2 BWA MEM!")
    result_list = []
    cpu_count = int(os.cpu_count()/2) - 6  # 为了不影响其他任务，这里减去6个核心，30个核心
    
    bam_files = find_files_by_suffix(".bam", project_dir)

    if len(bam_files) == 0:
        with concurrent.futures.ProcessPoolExecutor(max_workers=cpu_count) as executor:
            futures = [executor.submit(run_bwa_mem, sample_files['fq1'], sample_files['fq2'], ref_fa, sample_name, project_dir)
                    for sample_name, sample_files in sample_dict.items()]
    else:
        print("\tBAM files exist, skip analysis!")
    
    ######################################################################
    print("Start Step3 GATK!")
    maf_files = find_files_by_suffix(".MAF", project_dir)
    print(f"find maf_files: {len(maf_files)}")
    if len(maf_files) == 0:
        with concurrent.futures.ProcessPoolExecutor(max_workers=cpu_count) as executor:
            futures = [executor.submit(run_gatk, ref_fa, sample_name, project_dir)
                    for sample_name, sample_files in sample_dict.items()]
    else:
        print("\tMAF files exist, skip analysis!")

    ######################################################################
    print("Start Step4 deal MAF!")
    small_maf_files = find_files_by_suffix(".add.dbsnp.MAF.xls", project_dir)
    print("find small_maf_files: ", len(small_maf_files))
    if len(small_maf_files) == 0:
        with concurrent.futures.ProcessPoolExecutor(max_workers=cpu_count) as executor:
            futures = [executor.submit(deal_maf, sample_name, project_dir, output_templates_file, script_path, database_path, hpa)
                    for sample_name, sample_files in sample_dict.items()]
    else:
        print("\tsmall MAF files exist, skip analysis!")

    ######################################################################
    print("Start Step5 blood identity!")
    brief_table_files = find_files_by_suffix("brief.table.xlsx", project_dir)
    print("find brief_table_files: ", len(brief_table_files))
    if len(brief_table_files) == 0:
        with concurrent.futures.ProcessPoolExecutor(max_workers=cpu_count) as executor:
            futures = [executor.submit(blood_identy, sample_name, project_dir, database_path, script_path)
                    for sample_name, sample_files in sample_dict.items()]
            
            for future in concurrent.futures.as_completed(futures):
                result_list.extend(future.result())
    else:
        print("\tbrief table files exist, skip analysis!")
    
    
    ######################################################################
    print("Start Step6 extract CDS!")
    cds_files = find_files_by_suffix(".CDS.fasta", project_dir)
    print("find cds_files: ", len(cds_files))
    if len(cds_files) == 0:
        with concurrent.futures.ProcessPoolExecutor(max_workers=cpu_count) as executor:
            futures = [executor.submit(extract_cds, sample_name, project_dir, script_path, database_path)
                    for sample_name, sample_files in sample_dict.items()]
            
            for future in concurrent.futures.as_completed(futures):
                result_list.extend(future.result())
        cds_files = find_files_by_suffix(".CDS.fasta", project_dir)
        print("After: find cds_files: ", len(cds_files))
    else:
        print("\tCDS files exist, skip analysis!")

    
    ######################################################################
    if hpa:
        print("Start Step7 HPA report!")
        print("Start HPA report!")
        subprocess.run(f"python3 {script_path}/hpa.summary.report.py -i {project_dir} -d {database_path}/HPA.Gene.cDNA_Changes.xls -o {project_dir}/HPA.summary.report.xls", shell=True)
        
        # write file list that need to be packaged into a file.
        with open(os.path.join(project_dir, "need_to_be_packaged.txt"), "w") as f:
            f.write(f"{project_dir}/HPA.summary.report.xls\n")
            
            # 添加其他要打包的文件， 需要确认！
            for each in result_list:
                f.write(each + "\n")
    else:
        print("Start Step7 RBC report!")
        with open(os.path.join(project_dir, "need_to_be_packaged.txt"), "w") as f:
            for each in result_list:
                f.write(each + "\n")
    
    print("finished!")


# 在主脚本中
if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    # 定义位置参数
    parser.add_argument('data_dir', type=str, help='The directory where data is stored')
    parser.add_argument('project_name', type=str, help='The name of the project')
    parser.add_argument('output_templates_file', type=str, help='The output templates file path')
    parser.add_argument('result_path', type=str, help='The result path')

    # 定义可选参数
    parser.add_argument('--hpa', action='store_true', help='If this is a hpa project', default=False)

    # 解析命令行参数
    args = parser.parse_args()

    # 使用参数
    data_dir = args.data_dir
    project_dir = args.result_path

    BASE_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    software_path = os.path.join(BASE_DIR, "pipeline/software")
    database_path = os.path.join(BASE_DIR, "pipeline/database")

    script_path = os.path.join(BASE_DIR, "tools/scripts")
    ref_fa = os.path.join(BASE_DIR, "pipeline/ref/hg38/Homo_sapiens_assembly38.fasta")

    main(data_dir, project_dir, software_path, database_path, script_path, ref_fa, 
        args.output_templates_file, args.hpa)
