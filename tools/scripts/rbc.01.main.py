import subprocess
import sys
import os
import fnmatch


def find_fastq_files(directory):
    fastq_files = []
    for root, dirs, files in os.walk(directory):
        for filename in files:
            if fnmatch.fnmatch(filename, '*.fq.gz'):
                fastq_files.append(os.path.join(root, filename))
    return fastq_files


def process_fastq_files(fastq_files_list):
    sample_dict = {}
    
    for filepath in fastq_files_list:
        filename = os.path.basename(filepath)  # 获取文件名，如'CX1680_Raw_1.fq.gz'
        sample_name = filename.split('_')[0]   # 提取样本名称，如'CX1680'
        file_number = filename.split('_')[-1]  # 提取文件编号，如'1.fq.gz'
        file_key = 'fq1' if file_number == '1.fq.gz' else 'fq2'  # 判断文件类型
        
        if sample_name in sample_dict:
            sample_dict[sample_name][file_key] = filepath
        else:
            sample_dict[sample_name] = {file_key: filepath}
    
    return sample_dict

# project_name:  test-autoupload
# upload_dir:  /data/webapp/mysite/pipeline/project_hpa/8833a5c2-0fb1-4c13-b734-6acace996f08

upload_dir = sys.argv[1]
project_name = sys.argv[2]
output_templates_file = sys.argv[3]

project_dir = os.path.join(upload_dir , project_name)
os.makedirs(project_dir, exist_ok=True)

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
software_path = os.path.join(BASE_DIR, "software")
database_path = os.path.join(BASE_DIR, "database")
script_path = os.path.join(BASE_DIR, "scripts")
ref_fa = os.path.join(BASE_DIR, "ref/hg38/Homo_sapiens_assembly38.fasta")

# ## start analysis !!!!!
fastq_list = find_fastq_files(upload_dir)
sample_dict = process_fastq_files(fastq_list)

result_list = []
for sample_name in sample_dict.keys():
    fq1 = sample_dict[sample_name]['fq1']
    fq2 = sample_dict[sample_name]['fq2']
    sample_folder = os.path.join(project_dir, sample_name)
    os.makedirs(sample_folder, exist_ok=True)

    os.system(f"""bwa mem -t 8 -Y -H "@HD\\tVN:1.5\\tGO:none\\tSO:coordinate" -R "@RG\\tID:TEST\\tSM:TEST\\tLB:Target\\tPL:ILLumina\\tPU:HVW2MCCXX:6:none"  {ref_fa}  \
        {fq1} {fq2} | samtools view -buhS -t {ref_fa}.fai - | samtools sort -o {project_dir}/{sample_name}/{sample_name}.Rawsample.hg38.sortedByCoord.bam - """)
    os.system(f"""samtools index {project_dir}/{sample_name}/{sample_name}.Rawsample.hg38.sortedByCoord.bam""")
    os.system(f"""{software_path}/gatk --java-options "-Xmx30G" HaplotypeCaller -R {ref_fa} \
        -I {project_dir}/{sample_name}/{sample_name}.Rawsample.hg38.sortedByCoord.bam \
        -L {database_path}/Blood.gene.bed \
        -O {project_dir}/{sample_name}/{sample_name}.Rawsample.output.raw.vcf \
        --create-output-variant-index false""")
    os.system(f"""{software_path}/gatk VariantFiltration -R {ref_fa}  \
        -V {project_dir}/{sample_name}/{sample_name}.Rawsample.output.raw.vcf \
        -filter "QD < 2.0"  --filter-name "QDFilter" \
        -filter "MQ < 40.0" --filter-name "MQFilter" \
        -filter "FS > 60.0" --filter-name "FSFilter" \
        --create-output-variant-index false \
        --output {project_dir}/{sample_name}/{sample_name}.Rawsample.output.filtered.vcf """)

    os.system(f"""{software_path}/gatk --java-options "-Xmx30G"   Funcotator \
        --data-sources-path {software_path}/funcotator_dataSources.v1.7.20200521s \
        --ref-version hg38  \
        --output-file-format MAF   \
        --reference {ref_fa}      \
        --exclude-field Center \
        --exclude-field Tumor_Sample_Barcode \
        --exclude-field Matched_Norm_Sample_Barcode \
        --exclude-field Match_Norm_Seq_Allele1 \
        --exclude-field Match_Norm_Seq_Allele2 \
        --exclude-field Tumor_Validation_Allele1 \
        --exclude-field Tumor_Validation_Allele2 \
        --exclude-field Match_Norm_Validation_Allele1 \
        --exclude-field Match_Norm_Validation_Allele2 \
        --exclude-field Verification_Status \
        --exclude-field Validation_Status \
        --exclude-field Mutation_Status \
        --exclude-field Sequencing_Phase \
        --exclude-field Sequence_Source \
        --exclude-field Validation_Method \
        --exclude-field Score \
        --exclude-field BAM_File \
        --exclude-field Sequencer \
        --exclude-field Tumor_Sample_UUID \
        --exclude-field Matched_Norm_Sample_UUID \
        --variant {project_dir}/{sample_name}/{sample_name}.Rawsample.output.filtered.vcf \
        --output {project_dir}/{sample_name}/{sample_name}.Rawsample.funcotated.raw.MAF \
        --remove-filtered-variants true""")

    os.system(f"""grep -v "^#" {project_dir}/{sample_name}/{sample_name}.Rawsample.funcotated.raw.MAF \
        > {project_dir}/{sample_name}/{sample_name}.Rawsample.funcotated.MAF.xls""")

    os.system(f"""python3 {script_path}/deal_maf_table.py  -m {project_dir}/{sample_name}/{sample_name}.Rawsample.funcotated.MAF.xls \
                -o {project_dir}/{sample_name}/{sample_name}.Rawsample.funcotated.small.MAF.xls """)

    os.system(f"""python3 {script_path}/maf_add_dbSNP_HVGS.py \
                -m {project_dir}/{sample_name}/{sample_name}.Rawsample.funcotated.small.MAF.xls  \
                -s {database_path}/dbSNP.db.txt \
                -t {output_templates_file} \
                -o {project_dir}/{sample_name}/{sample_name}.Rawsample.funcotated.small.add.dbsnp.MAF.xls""")

    os.system(f"""python3 {script_path}/Blood_group_identification.py \
                -b {database_path}/Blood.Gene.metadata.manually.checked.v2.xlsx \
                -p {database_path}/HPA.Gene.cDNA_Changes.xls \
                -m {project_dir}/{sample_name}/{sample_name}.Rawsample.funcotated.small.add.dbsnp.MAF.xls \
                -o {project_dir}/{sample_name}/{sample_name}.Rawsample.funcotated.brief.table.xlsx """)
    os.system(f"""python3 {script_path}/extractCDS.py \
              {project_dir}/{sample_name}/{sample_name}.Rawsample.funcotated.small.add.dbsnp.MAF.xls \
              {database_path}/RBC.HPA.CDS.60.fasta \
              {project_dir}/{sample_name}/{sample_name}.Rawsample.CDS.fasta""")
    
    result_list.append(f"{project_dir}/{sample_name}/{sample_name}.Rawsample.funcotated.brief.table.xlsx")
    result_list.append(f"{project_dir}/{sample_name}/{sample_name}.Rawsample.CDS.fasta")

##########################################################################
# write file list that need to be packaged into a file.
with open(os.path.join(project_dir, "need_to_be_packaged.txt"), "w") as f:
    for each in result_list:
        f.write(each + "\n")