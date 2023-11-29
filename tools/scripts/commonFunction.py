import os

def find_files_by_suffix(suffix, path):
    result = []
    for root, dirs, files in os.walk(path):
        for file in files:
            if file.endswith(suffix):
                result.append(os.path.join(root, file))
    return result

def process_fastq_files(fastq_files_list):
    '''
        将fastq文件按照样本名称进行分组，返回一个字典，
        字典的key为样本名称，value为一个字典，包含两个key，分别为'fq1'和'fq2'，value为fastq文件的路径
        
        sample_dict = {
            'CX1680_Raw': {
                'fq1': '/data/webapp/mysite/pipeline/project_hla/CX1680_Raw_1.fq.gz',
                'fq2': '/data/webapp/mysite/pipeline/project_hla/CX1680_Raw_2.fq.gz'
            },
            'CX1681_Raw': {
                'fq1': '/data/webapp/mysite/pipeline/project_hla/CX1681_Raw_1.fq.gz',
                'fq2': '/data/webapp/mysite/pipeline/project_hla/CX1681_Raw_2.fq.gz'
            }
        }
    '''
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