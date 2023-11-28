#_*_coding:utf-8_*_
"""
@Time : 20221102
@Author: dushiyi
@Email: dushiyi319@163.com
@File: main.py
@IDE : PyCharm
"""

#!/usr/bin/python3.8
import os
import shutil
import re
import argparse
import subprocess
import configparser
import logging
import sys


def _argparse():
    parser = argparse.ArgumentParser(description="Version: 1.0")
    parser.add_argument('-c', '--config', action='store', dest='config', required=True, default='config.txt', help="this is config file")
    parser.add_argument('--quiet', action='store_true', help='Disable logging')
    parser.add_argument('--verbose', action='store_true', help='Enable verbose logging')
    return parser.parse_args()


def make_dir(*dir):
    for each in dir:
        if not os.path.exists(each):
            os.makedirs(each)


def check_cf(cf):
    # result = cf.sections()
    sections = ['data','trim', 'mapping','hla']
    for section in sections:
        if not cf.has_section(section):
            logging.error("\033[1;31mError:\033[0m your config file not has [{}] section".format(section))


def check_fq(fq1,fq2):
    if not os.path.exists(fq1) or not os.path.exists(fq2):
        logging.error("\033[1;31mError:\033[0m can not find fastq : {},{}".format(fq1,fq2))
    if fq1 == fq2:
        logging.error("\033[1;31mError:\033[0m same fastqs, Please Check: {},{}".format(fq1,fq2))
    logging.info("Check Data: {} {} \033[0;36mSuccessful.\033[0m".format(fq1, fq2))


def check_install(soft):
    value = subprocess.getstatusoutput("which " + soft)
    if value[0] != 0 :
        logging.error("\033[1;31mError: \033[0m Can not find software '{}' Please check installed path.".format(soft))
    else:
        logging.info("Check installation path {} : {} \033[0;36mSuccessful.\033[0m".format(os.path.basename(soft), soft))


def section_fastqc(cf, output, project_dir, fq1, fq2):
    fastqc = cf.get("data","fastqc")
    if fastqc in ['True', 'TRUE','true','T']:
        logging.info("Check Module: Fastqc : \033[0;36mTrue\033[0m")
        output.write("\n echo \"# Fastc\"\n")
        fastqc_dir = os.path.join(project_dir, "Fastqc")
        make_dir(fastqc_dir)
        fastqc_soft = cf.get("data", "fastqcsoft")
        check_install(fastqc_soft)
        output.write(f'''if [ ! -f "{fastqc_dir}/{fq2}_1_fastqc.gz" ]; then\n''')
        output.write(f"\t{fastqc_soft} --threads {cf.get('data','fastqc_threads')} --format fastq --outdir {fastqc_dir} --quiet {fq1} {fq2} \n")
        output.write("fi\n")
    else:
        logging.warning("Check Module: Fastqc : \033[0;36mFalse\033[0m")


def section_trim(cf, output, project_dir, fq1, fq2, project_name):
    trimdata = cf.get("trim", "trimdata") 
    if trimdata in ["True","T"]:
        logging.info("Check Module: Trim data : \033[0;36mTrue\033[0m")
        output.write("\necho \"# Trim data\"\n")
        # mkdir for trimed data
        trim_data_dir = os.path.join(project_dir, "downsample")
        make_dir(trim_data_dir)
        seqkit = cf.get("trim", "trimsoft")
        check_install(seqkit)
        trim_string = cf.get("trim", "trim")   # 1M
        adict = {
                'fq1': fq1,
                'fq2': fq2,
                'seqkit': seqkit,
                'trim_data_dir': trim_data_dir,
                'trim_string': trim_string,
                'trim_number': int(float(trim_string.rstrip("M"))*1000000),  # 1000000
                'seqkit_seed': cf.get("trim", "seqkit_seed"),
                'trim_threads': cf.get("trim", "trim_threads"),
                'trim_portion': cf.get("trim", "trim_portion"),
        }
        
        ### sampling
        output.write(f"if [ ! -f \"{trim_data_dir}/Rmdup.Downsample{trim_string}_1.fq\" ]; then\n")
        output.write("\t{seqkit} sample -j {trim_threads} -p {trim_portion} -s {seqkit_seed} {fq1} | {seqkit} head -n {trim_number} >{trim_data_dir}/Rmdup.Downsample{trim_string}_1.fq\n".format(**adict))
        output.write("\t{seqkit} sample -j {trim_threads} -p {trim_portion} -s {seqkit_seed} {fq2} | {seqkit} head -n {trim_number} >{trim_data_dir}/Rmdup.Downsample{trim_string}_2.fq\n".format(**adict))
        output.write("fi\n")
        fq1 = os.path.join(trim_data_dir, "Rmdup.Downsample{}_1.fq".format(trim_string))
        fq2 = os.path.join(trim_data_dir, "Rmdup.Downsample{}_2.fq".format(trim_string))
        sample_name = project_name + ".Downsample"
    else:
        logging.info("Check Module: Trim data : \033[0;36mFalse\033[0m")
        sample_name = project_name +".Rawsample"
    return sample_name, fq1, fq2


def section_mapping(cf, output, project_dir, base_dir, fq1,fq2, sample_name):
    mapping = cf.get("mapping","mapping")
    ref_name = None
    if mapping in ["True",'T','TRUE','true']:
        logging.info("Check Module: Mapping : \033[0;36mTrue\033[0m")
        output.write("\necho \"# BWA\" \n")
        mapping_dir = os.path.join(project_dir, "mapping")
        make_dir(mapping_dir)
        
        bwa = cf.get("mapping", "bwa")
        samtools = cf.get("mapping", "samtools")
        check_install(bwa)
        check_install(samtools)
        # get mapping soft
        adict = {
                "mapping_dir": mapping_dir,
                "sample_name": sample_name,
                "bwa": bwa,
                "bwa_threads": cf.get("mapping", "bwa_threads"),
                "samtools": samtools,
                "samtools_threads": cf.get("mapping", "samtools_threads"),
                "bwa_header": "@HD\\tVN:1.5\\tGO:none\\tSO:coordinate",
                "bwa_RG": "@RG\\tID:TEST\\tSM:TEST\\tLB:Target\\tPL:ILLumina\\tPU:HVW2MCCXX:6:none",
                "ref_fa": os.path.join(base_dir, cf.get("mapping", "bwa_ref")),
        }
        # 20230107
        a = re.match(r'.*/ref/(\w+\-?\w+)/.*',adict["ref_fa"])
        if a is not None:
            ref_name = a.group(1).replace("-","")
        else:
            logging.error("""\033[1;31mError:\033[0m can not find ref_name in [mapping] section, bwa_ref = {}""".format(adict["ref_fa"]))
            sys.exit(1)
        
        # end
        output.write(f"if [ ! -f \"{mapping_dir}/{sample_name}.{ref_name}.sortedByCoord.bam\" ]; then\n")
        output.write("\t{bwa} mem -t {bwa_threads} -Y -H \"{bwa_header}\" -R \"{bwa_RG}\" {ref_fa} {fq1} {fq2} | {samtools} view -buhS -t {ref_fa}.fai - | {samtools} sort -o {mapping_dir}/{sample_name}.{ref_name}.sortedByCoord.bam - \n".format(**{"fq1":fq1,"fq2":fq2,"ref_name":ref_name},**adict))
        output.write("\t{samtools} view -@ {samtools_threads} -bF 4 {mapping_dir}/{sample_name}.{ref_name}.sortedByCoord.bam > {mapping_dir}/{sample_name}.mapped.{ref_name}.bam\n".format(**{"ref_name": ref_name},**adict))
        output.write("\t{samtools} view -@ {samtools_threads} -bf 4 {mapping_dir}/{sample_name}.{ref_name}.sortedByCoord.bam > {mapping_dir}/{sample_name}.Unmapped.{ref_name}.bam\n".format(**adict,**{"ref_name":ref_name}))
        output.write("\t{samtools} fastq -@ {samtools_threads} -1 {mapping_dir}/{sample_name}.mapped.{ref_name}.1.fastq -2 {mapping_dir}/{sample_name}.mapped.{ref_name}.2.fastq -n {mapping_dir}/{sample_name}.mapped.{ref_name}.bam\n".format(**adict,**{"ref_name":ref_name}))
        output.write("\t{samtools} index {mapping_dir}/{sample_name}.{ref_name}.sortedByCoord.bam\n".format(**adict,**{"ref_name":ref_name}))
        output.write("fi\n")
        fq1 = os.path.join(mapping_dir, "{sample_name}.mapped.{ref_name}.1.fastq".format(**adict,**{"ref_name":ref_name}))
        fq2 = os.path.join(mapping_dir, "{sample_name}.mapped.{ref_name}.2.fastq".format(**adict,**{"ref_name":ref_name}))
    elif mapping in ["False","F", "false","FALSE"]:
        logging.info("Check Module: Mapping : \033[0;36mFalse\033[0m")
    else:
        logging.error("\033[1;31mError:\033[0m Unrecognized symbol: '{}' in [mapping] section, mapping = True. only [True or False]".format(mapping))
    return sample_name, fq1, fq2, ref_name


def section_hla(cf,output, project_dir, base_dir, fq1,fq2, sample_name):
    hla = cf.get("hla", "hla")
    if hla in ["True",'T','TRUE','true']:
        logging.info("Check Module: HLA : \033[0;36mTrue\033[0m")
        hla_soft_string = cf.get("hla", "soft")
        if hla_soft_string == "":
            logging.info("\033[1;31mWarning:\033[0m 'soft' section is empty in config file. using default options instead: [hla_scan, hlahd, OptiType]")
            hla_soft_list = ["hla_scan", "hlahd", "OptiType"]
        else:
            hla_soft_list = re.sub(r"\s+", "", hla_soft_string).split(",")
            logging.info("Check HLA options: {}".format(hla_soft_list))

        hla_dir = os.path.join(project_dir, "HLA")
        make_dir(hla_dir)
        config_dict = {
            "fq1": fq1,
            "fq2": fq2,
            "sample_name": sample_name,
            "hla_dir": hla_dir,
            "base_dir": base_dir,
        }

        for each in hla_soft_list:
            if each == "hlahd":
                hlahd = cf.get("hla", "hlahd")
                check_install(hlahd)

                config_dict.update({
                    "hlahd": hlahd,
                    "hlahd_threads": cf.get("hla","hlahd_threads"),
                    "hlahd_ref": cf.get("hla","hlahd_ref"),
                    "split_file": cf.get("hla", "hlahd_gene_split_txt")
                    })
                output.write("\n# hlahd\n")
                output.write("\necho '# Starting hlahd'\n")
                output.write("{hlahd} -t {hlahd_threads} -m 100 -c 0.95 -f {hlahd_ref}/freq_data {fq1} {fq2}  {hlahd_ref}/{split_file} {hlahd_ref}/dictionary/ HLA-HD_Result {hla_dir}\n".format(**config_dict))
                output.write("python3 {base_dir}/scripts/hla.freq.py {base_dir}/database/hla.cwd.xls {hla_dir}/HLA-HD_Result/result/HLA-HD_Result_final.result.txt\n".format(**config_dict))
                output.write("\necho '# Finish hlahd'\n ")
            elif each == "hla_scan":
                hla_scan = cf.get("hla","hla_scan")
                check_install(hla_scan)
                config_dict.update({
                    "hla_scan": hla_scan,
                    "hla_scan_db": cf.get("hla","hla_scan_db"),
                    "hla_scan_threads": cf.get("hla", "hla_scan_threads"),
                })
                output.write("\n# hla_scan\n")
                hla_scan_dir = os.path.join(project_dir, "HLA", "HLAscan_Result")
                make_dir(hla_scan_dir)

                output.write("\necho '# Starting hla_scan'\n")
                hla_genes = ["HLA-A","HLA-B","HLA-C","HLA-DMA","HLA-DMB","HLA-DOA","HLA-DOB","HLA-DPA1","HLA-DPB1","HLA-DQA1","HLA-DQB1","HLA-DRA","HLA-DRB1","HLA-DRB5","HLA-E","HLA-F","HLA-G","MICA","MICB","TAP1","TAP2"]
                for gene in hla_genes:
                    output.write("{hla_scan} -t {hla_scan_threads} -l {fq1} -r {fq2} -d {hla_scan_db} -g {gene} >{hla_scan_dir}/{sample_name}.{gene}.out.txt\n".format(**{"gene": gene, "hla_scan_dir": hla_scan_dir}, **config_dict))
                output.write("python3 {base_dir}/scripts/merge_HLA_result.py {hla_scan_dir}/{sample_name}*.out.txt > {hla_scan_dir}/HLAscan.results.txt\n".format(**{"hla_scan_dir": hla_scan_dir},**config_dict))
                output.write("\necho '# Finish hla_scan'\n")

            elif each == "OptiType":
                opt_dir = os.path.join(hla_dir, "OptiType_Result")
                make_dir(opt_dir)
                OptiType = cf.get("hla","OptiType")
                razers3 = cf.get("hla", "razers3")
                check_install(OptiType)
                check_install(razers3)
                config_dict.update({
                    "OptiType": OptiType,
                    "opt_dir" : opt_dir,
                    "OptiType_ref": cf.get("hla", "OptiType_ref"),
                    "razers3": razers3,
                    "razers3_threads": cf.get("hla", "razers3_threads"),
                })
                output.write("\necho '# Starting OptiType'\n")
                output.write("{razers3} -tc {razers3_threads} -i 95 -m 1 -dr 0 -o {opt_dir}/sample.fished_1.bam {OptiType_ref} {fq1}\n".format(**config_dict))
                output.write("{razers3} -tc {razers3_threads} -i 95 -m 1 -dr 0 -o {opt_dir}/sample.fished_2.bam {OptiType_ref} {fq2}\n".format(**config_dict))
                output.write("{samtools} bam2fq {opt_dir}/sample.fished_1.bam > {opt_dir}/sample.fished_1.fastq\n".format(**{"samtools": cf.get("mapping","samtools")}, **config_dict))
                output.write("{samtools} bam2fq {opt_dir}/sample.fished_2.bam > {opt_dir}/sample.fished_2.fastq\n".format(**{"samtools": cf.get("mapping","samtools")}, **config_dict))
                output.write("python3 {OptiType} -i {opt_dir}/sample.fished_1.fastq {opt_dir}/sample.fished_2.fastq --dna -v -o {opt_dir} -p {sample_name}\n".format(**config_dict))
                output.write("\necho '# Finish OptiType'\n ")
                # output.write("rm {opt_dir}/*.bam {opt_dir}/*.fastq".format(**config_dict))
            else:
                shutil.rmtree(project_dir, ignore_errors=True, onerror=None)
                logging.error("\033[0;36mError: \033[0m Wrong software name:[{}]. must be same as [hla_scan, hla_hd, OptiType]!".format(each)) 
    elif hla in ["F", "FALSE", "false", "False"]:
        logging.info("Check Module: HLA : \033[0;36mFalse\033[0m")
        # output.write('''skip HLA ...''')
    else:
        logging.error("\033[1;31mError: \033[0m Unrecognized symbol: '{}' in [hla] section, hla = True. only [True or False]".format(hla))


def section_gatk(cf, output, project_dir, base_dir, sample_name, ref_version):
    gatk = cf.get("GATK", "GATK")
    if gatk in ["True",'T','TRUE','true']:
        # check mapping section
        if cf.get("mapping","mapping") not in ["True",'T','TRUE','true']:
            logging.error("please config mapping section correctly before using gatk section")
        logging.info("Check Module: GATK : \033[0;36mTrue\033[0m")
        output.write("\necho \"# GATK\" \n")

        gatk_dir = os.path.join(project_dir, "gatk")
        make_dir(gatk_dir)
        
        samtools = cf.get("mapping", "samtools")
        mapping_dir = os.path.join(project_dir, "mapping")
        config_dict = {
            'gatk': cf.get("GATK","gatk_path"),
            'gatk_dir': gatk_dir,
            'mapping_dir': mapping_dir,
            'base_dir':  base_dir,
            'sample_name': sample_name,
            'ref_fa': cf.get("mapping", "bwa_ref"),
            'ref_version': ref_version,
            'gatk_dataSource': cf.get("GATK","functator"),
            'bed': cf.get("GATK", "bed"),
            'blood_table': cf.get("GATK", "blood_table"),
            'blood_metadata': cf.get("GATK", "blood_metadata"),
        }
        # HaplotypeCaller
        output.write(f'''if [ ! -f "{gatk_dir}/{sample_name}.output.raw.vcf" ]; then\n''')
        output.write("""\t{gatk} --java-options \"-Xmx30G\" HaplotypeCaller \\
            -R {ref_fa} \\
            -I {mapping_dir}/{sample_name}.{ref_version}.sortedByCoord.bam \\
            -L {bed} \\
            -O {gatk_dir}/{sample_name}.output.raw.vcf \\
            --create-output-variant-index false\n""".format(**config_dict))
        output.write("fi\n")
        # SelectVariants
        # filter

        output.write(f'''if [ ! -f "{gatk_dir}/{sample_name}.output.filtered.vcf" ]; then\n''')
        output.write("""\t{gatk} VariantFiltration \\
            -R {ref_fa} \\
            -V {gatk_dir}/{sample_name}.output.raw.vcf \\
            -filter \"QD < 2.0\"  --filter-name \"QDFilter\" \\
            -filter \"MQ < 40.0\" --filter-name \"MQFilter\" \\
            -filter \"FS > 60.0\" --filter-name \"FSFilter\" \\
            --create-output-variant-index false \\
            --output {gatk_dir}/{sample_name}.output.filtered.vcf\n\n""".format(**config_dict))
        output.write("fi\n")
        # Funcotator
        output.write(f'''if [ ! -f "{gatk_dir}/{sample_name}.funcotated.raw.MAF" ]; then\n''')
        output.write("""\t{gatk} --java-options \"-Xmx30G\"   Funcotator \\
            --data-sources-path {gatk_dataSource} \\
            --ref-version {ref_version}  \\
            --output-file-format MAF   \\
            --reference {ref_fa}    \\
            --exclude-field Center \\
            --exclude-field Tumor_Sample_Barcode \\
            --exclude-field Matched_Norm_Sample_Barcode \\
            --exclude-field Match_Norm_Seq_Allele1 \\
            --exclude-field Match_Norm_Seq_Allele2 \\
            --exclude-field Tumor_Validation_Allele1 \\
            --exclude-field Tumor_Validation_Allele2 \\
            --exclude-field Match_Norm_Validation_Allele1 \\
            --exclude-field Match_Norm_Validation_Allele2 \\
            --exclude-field Verification_Status \\
            --exclude-field Validation_Status \\
            --exclude-field Mutation_Status \\
            --exclude-field Sequencing_Phase \\
            --exclude-field Sequence_Source \\
            --exclude-field Validation_Method \\
            --exclude-field Score \\
            --exclude-field BAM_File \\
            --exclude-field Sequencer \\
            --exclude-field Tumor_Sample_UUID \\
            --exclude-field Matched_Norm_Sample_UUID \\
            --variant {gatk_dir}/{sample_name}.output.filtered.vcf \\
            --output {gatk_dir}/{sample_name}.funcotated.raw.MAF \\
            --remove-filtered-variants true\n\n""".format(**config_dict))
        output.write("fi\n")

        # remove MAF header
        output.write(f'''if [ ! -f "{gatk_dir}/{sample_name}.funcotated.MAF.xls" ]; then\n''')
        output.write('''\tgrep -v \"^#\" {gatk_dir}/{sample_name}.funcotated.raw.MAF \\
            > {gatk_dir}/{sample_name}.funcotated.MAF.xls\n\n'''.format(**config_dict))
        output.write("fi\n")
        
        # filter MAF table
        output.write('''python3 {base_dir}/scripts/deal_maf_table.py \\
            -m {gatk_dir}/{sample_name}.funcotated.MAF.xls \\
            -o {gatk_dir}/{sample_name}.funcotated.small.MAF.xls \n\n'''.format(**config_dict))
        # add dbsnp HVGS info
        output.write('''python3 {base_dir}/scripts/maf_add_dbSNP_HVGS.py \\
            -m {gatk_dir}/{sample_name}.funcotated.small.MAF.xls \\
            -s {base_dir}/database/dbSNP.db.txt \\
            -t {blood_table} \\
            -o {gatk_dir}/{sample_name}.funcotated.small.add.dbsnp.MAF.xls\n\n'''.format(**config_dict))
        # identify blood group info
        ref_dir = os.path.dirname(config_dict['ref_fa'])
        output.write("""python3 {base_dir}/scripts/Blood_group_identification.py \\
            -b {blood_metadata} \\
            -p {base_dir}/database/HPA.Gene.cDNA_Changes.xls \\
            -m {gatk_dir}/{sample_name}.funcotated.small.add.dbsnp.MAF.xls \\
            -o {gatk_dir}/{sample_name}.funcotated.brief.table.xlsx \n\n""".format(**{'ref_dir':ref_dir},**config_dict))
        output.write("""python3 {base_dir}/scripts/extractCDS.py \\
            {gatk_dir}/{sample_name}.funcotated.small.add.dbsnp.MAF.xls \\
            {base_dir}/database/RBC.HPA.CDS.60.fasta \\
            {gatk_dir}/{sample_name}.CDS.fasta\n\n""".format(**config_dict))
        output.write("echo gatk done!\n")


def main():
    """docstring for main"""
    parser = _argparse()

    if parser.quiet:
        logging.disable(logging.CRITICAL)  # 禁用所有日志输出
    elif parser.verbose:
        logging.basicConfig(level=logging.DEBUG, format="%(asctime)s %(levelname)s %(message)s", datefmt = '%Y-%m-%d %H:%M:%S')  # 设置日志级别为DEBUG，打印详细信息
    else:
        logging.basicConfig(level=logging.DEBUG, format="%(asctime)s %(levelname)s %(message)s", datefmt = '%Y-%m-%d %H:%M:%S')

    cf = configparser.ConfigParser(allow_no_value=True)
    cf.read(parser.config)
    
    ############################### check config file
    check_cf(cf)
    
    # Global variables 
    base_dir    = os.path.dirname(os.path.abspath(__file__))
    project_name = cf.get("data", "project_name")
    project_dir = os.path.join(os.getcwd(), project_name)
    
    make_dir(project_dir)
    
    output = open(os.path.join(project_dir, "work.sh"), "w")
    output.write("""time1=$(date "+%Y-%m-%d %H:%M:%S")\necho $time1 "Start Running."\n""")
    
    ############################### Check fastq
    fq1 = cf.get("data", "fastq1")
    fq2 = cf.get("data", "fastq2")
    check_fq(fq1, fq2)
    
    ### move fastq to rawdata and rename the header.
    # data_path = os.path.join(project_dir, "data")
    # make_dir(data_path)
    # fixed_fq1 = str(os.path.basename(fq1)).replace("fq.gz", "fixed.fq.gz")
    # fixed_fq2 = str(os.path.basename(fq2)).replace("fq.gz", "fixed.fq.gz")
    # output.write('''zcat {fq1} |awk '{{if(NR%4==1){{O=$0;gsub("/1", " 1", O);print O}}else{{print$0}}}}' |gzip > {data_path}/{fixed_fq1}\n'''.format(**{
    #     'fq1':fq1,'data_path':data_path,'fixed_fq1':fixed_fq1}))
    # output.write('''zcat {fq2} |awk '{{if(NR%4==1){{O=$0;gsub("/2", " 2", O);print O}}else{{print$0}}}}' |gzip > {data_path}/{fixed_fq2}\n'''.format(**{
    #     'fq2':fq2,'data_path':data_path,'fixed_fq2':fixed_fq2}))
    # fq1 = os.path.join(data_path, fixed_fq1)
    # fq2 = os.path.join(data_path, fixed_fq2)
    ############################### Fastqc
    section_fastqc(cf, output, project_dir, fq1, fq2)
    
    ############################### Trim
    sample_name, fq1, fq2 = section_trim(cf, output, project_dir, fq1, fq2, project_name)
    # print("sample name will be used in following section: {}".format(sample_name))
    
    ############################### Mapping
    sample_name, fq1, fq2, ref_version = section_mapping(cf, output, project_dir, base_dir, fq1, fq2, sample_name)
    # print("sample name will be used in following section: {}".format(sample_name))
    
    ############################### HLA
    section_hla(cf, output, project_dir, base_dir, fq1,fq2, sample_name)

    ############################### GATK and annotation
    section_gatk(cf, output, project_dir, base_dir, sample_name, ref_version)
    
    ############################### Finished.
    logging.info("""
        -- All Done --

    Run following commands:
        `cd {} `
        `sh work.sh` or `nohup sh work.sh &`
        """.format(cf.get("data", "project_name")))
    output.write("""time1=$(date "+%Y-%m-%d %H:%M:%S")\necho $time1 "Finished Run." \n""")
    output.close()


if __name__ == '__main__':
    main()
