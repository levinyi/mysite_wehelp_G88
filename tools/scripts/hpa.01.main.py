import subprocess
import sys
import os
import fnmatch
import concurrent.futures
from commonFunction import find_files_by_suffix, process_fastq_files


def analyze_sample(sample_name, sample_files, project_dir, software_path, database_path, script_path, ref_fa, output_templates_file):
    fq1 = sample_files['fq1']
    fq2 = sample_files['fq2']
    sample_folder = os.path.join(project_dir, sample_name)
    os.makedirs(sample_folder, exist_ok=True)

    threads = os.cpu_count()

    # BWA MEM command
    print("BWA MEM Start!")
    subprocess.run([
        f"bwa mem -t {threads} -v 1 -Y -H \"@HD\\tVN:1.5\\tGO:none\\tSO:coordinate\" -R \"@RG\\tID:TEST\\tSM:TEST\\tLB:Target\\tPL:BGI\\tPU:HVW2MCCXX:6:none\"  {ref_fa} {fq1} {fq2} | samtools view -buhS -t {ref_fa}.fai - | samtools sort -o {project_dir}/{sample_name}/{sample_name}.Rawsample.hg38.sortedByCoord.bam -"],
        shell=True)
    
    # Samtools index command
    subprocess.run([
        f"samtools index {project_dir}/{sample_name}/{sample_name}.Rawsample.hg38.sortedByCoord.bam"],
        shell=True)

    # GATK HaplotypeCaller command
    subprocess.run([
        f"{software_path}/gatk/gatk --java-options \"-Xmx30G\" HaplotypeCaller -R {ref_fa} -I {project_dir}/{sample_name}/{sample_name}.Rawsample.hg38.sortedByCoord.bam -L {database_path}/Blood.gene.bed -O {project_dir}/{sample_name}/{sample_name}.Rawsample.output.raw.vcf --QUIET true --verbosity ERROR --create-output-variant-index false"],
        shell=True)

    # GATK VariantFiltration command
    subprocess.run([
        f"{software_path}/gatk/gatk VariantFiltration -R {ref_fa} -V {project_dir}/{sample_name}/{sample_name}.Rawsample.output.raw.vcf -filter \"QD < 2.0\"  --filter-name \"QDFilter\" -filter \"MQ < 40.0\" --filter-name \"MQFilter\" -filter \"FS > 60.0\" --filter-name \"FSFilter\" --create-output-variant-index false --output {project_dir}/{sample_name}/{sample_name}.Rawsample.output.filtered.vcf"],
        shell=True)

    # GATK Funcotator command
    subprocess.run([
        f"{software_path}/gatk/gatk --java-options \"-Xmx30G\" Funcotator --data-sources-path {database_path}/funcotator_dataSources.v1.7.20200521s --ref-version hg38 --output-file-format MAF --reference {ref_fa} --exclude-field Center --exclude-field Tumor_Sample_Barcode --exclude-field Matched_Norm_Sample_Barcode --exclude-field Match_Norm_Seq_Allele1 --exclude-field Match_Norm_Seq_Allele2 --exclude-field Tumor_Validation_Allele1 --exclude-field Tumor_Validation_Allele2 --exclude-field Match_Norm_Validation_Allele1 --exclude-field Match_Norm_Validation_Allele2 --exclude-field Verification_Status --exclude-field Validation_Status --exclude-field Mutation_Status --exclude-field Sequencing_Phase --exclude-field Sequence_Source --exclude-field Validation_Method --exclude-field Score --exclude-field BAM_File --exclude-field Sequencer --exclude-field Tumor_Sample_UUID --exclude-field Matched_Norm_Sample_UUID --variant {project_dir}/{sample_name}/{sample_name}.Rawsample.output.filtered.vcf --output {project_dir}/{sample_name}/{sample_name}.Rawsample.funcotated.raw.MAF --remove-filtered-variants true"],
        shell=True)

    # GREP command
    subprocess.run([
        f"grep -v \"^#\" {project_dir}/{sample_name}/{sample_name}.Rawsample.funcotated.raw.MAF > {project_dir}/{sample_name}/{sample_name}.Rawsample.funcotated.MAF.xls"],
        shell=True)

    # Python scripts execution
    subprocess.run([
        f"python3 {script_path}/deal_maf_table.py  -m {project_dir}/{sample_name}/{sample_name}.Rawsample.funcotated.MAF.xls -o {project_dir}/{sample_name}/{sample_name}.Rawsample.funcotated.small.MAF.xls"],
        shell=True)

    subprocess.run([
        f"python3 {script_path}/maf_add_dbSNP_HVGS.py -m {project_dir}/{sample_name}/{sample_name}.Rawsample.funcotated.small.MAF.xls -s {database_path}/dbSNP.db.txt -t {output_templates_file} -o {project_dir}/{sample_name}/{sample_name}.Rawsample.funcotated.small.add.dbsnp.MAF.xls"],
        shell=True)

    subprocess.run([
        f"python3 {script_path}/Blood_group_identification.py -b {database_path}/Blood.Gene.metadata.manually.checked.v2.xlsx -m {project_dir}/{sample_name}/{sample_name}.Rawsample.funcotated.small.add.dbsnp.MAF.xls -o {project_dir}/{sample_name}/{sample_name}.Rawsample.funcotated.brief.table.xlsx"],
        shell=True)

    result_list = [f"{project_dir}/{sample_name}/{sample_name}.Rawsample.funcotated.brief.table.xlsx"]
    return result_list


def main(data_dir, project_dir, software_path, database_path, script_path, ref_fa, output_templates_file):
    fastq_list  = find_files_by_suffix(".fq.gz", data_dir)
    sample_dict = process_fastq_files(fastq_list)
    print("sample_dict: ", sample_dict)

    ### 主要分析步骤
    print("Start multiprocess analysis!")
    result_list = []
    with concurrent.futures.ProcessPoolExecutor() as executor:
        futures = [executor.submit(analyze_sample, sample_name, sample_files, project_dir, software_path, database_path, script_path, ref_fa, output_templates_file) 
                   for sample_name, sample_files in sample_dict.items()]

        for future in concurrent.futures.as_completed(futures):
            result_list.extend(future.result())

    subprocess.run(f"python3 {script_path}/hpa.summary.report.py -i {project_dir} -o {project_dir}/{project_name}.Rawsample.funcotated.brief.table.xlsx", shell=True)
    ##########################################################################
    # write file list that need to be packaged into a file.
    with open(os.path.join(project_dir, "need_to_be_packaged.txt"), "w") as f:
        f.write(f"{project_dir}/{project_name}.Rawsample.funcotated.brief.table.xlsx\n")
        
        # 添加其他要打包的文件， 需要确认！
        # for each in result_list:
        #     f.write(each + "\n")
    print("finished!")


if __name__ == "__main__":
    data_dir     = sys.argv[1]
    project_name = sys.argv[2]
    output_templates_file = sys.argv[3]
    project_dir = sys.argv[4]
    os.makedirs(project_dir, exist_ok=True)

    BASE_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    software_path = os.path.join(BASE_DIR, "pipeline/software")
    database_path = os.path.join(BASE_DIR, "pipeline/database")
    ref_fa = os.path.join(BASE_DIR, "pipeline/ref/hg38/Homo_sapiens_assembly38.fasta")
    script_path = os.path.join(BASE_DIR, "tools/scripts")
    print(f"data_dir: {data_dir}")
    print(f"project_name: {project_name}")
    print(f"BASE_DIR: {BASE_DIR}")
    print(f"softwarepath: {software_path}")

    main(data_dir, project_dir, software_path, database_path, script_path, ref_fa, output_templates_file)
