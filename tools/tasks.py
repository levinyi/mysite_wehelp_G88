from my_celery import app
import os
import rarfile
import zipfile
import shutil
from .models import Result
from django.contrib.auth.models import User
from django.utils import timezone
import subprocess

@app.task()
def update_task_status(status, unique_id, user_id, end_time=None):
    # 更新数据库中的任务状态和数据文件路径
    user = User.objects.get(id=user_id)
    result = Result.objects.get(user=user, unique_id=unique_id)
    result.status = status
    if end_time:
        result.end_time = end_time
    result.save()

@app.task()
def extract_archive(archive_file, project_folder):
    file_extension = os.path.splitext(archive_file)[1].lower()
    if file_extension == '.zip':
        # 解压ZIP文件到指定目录
        with zipfile.ZipFile(archive_file, 'r') as zip_ref:
            for member in zip_ref.namelist():
                filename = os.path.basename(member)
                if not filename:
                    continue
                source = zip_ref.open(member)
                target = open(os.path.join(project_folder, filename), "wb")
                with source, target:
                    shutil.copyfileobj(source, target)
        
    elif file_extension == '.rar':
        # 解压RAR文件到指定目录
        with rarfile.RarFile(archive_file) as rf:
            for member in rf.namelist():
                filename = os.path.basename(member)
                if not filename:
                    continue
                source = rf.open(member)
                target = open(os.path.join(project_folder, filename), "wb")
                with source, target:
                    shutil.copyfileobj(source, target)
    else:
        print("Unsupported archive file type.")


@app.task()
def main_task(user_id, tools_id, unique_id):
    result = Result.objects.get(user=user_id, tools_name=tools_id, unique_id=unique_id)
    project_name = result.project_name
    upload_dir = result.result_path
    zip_file_path = os.path.join(upload_dir, os.path.basename(result.data.name))
    
    # print("project_name: ", project_name) # test_celery
    # print("upload_dir: ", upload_dir) # /data/webapp/mysite/pipeline/project_hla/486f110f-4244-4a6b-9a36-80ae4ef2533e
    # print("zip_file_path: ", zip_file_path) # /data/webapp/mysite/pipeline/project_hla/486f110f-4244-4a6b-9a36-80ae4ef2533e/HLA_TEST_DATA.rar

    # 更新任务状态为"extract_archive" &  
    extract_archive(zip_file_path, upload_dir)
    # 更新任务状态为"in_progress" &  执行Python脚本
    update_task_status('in_progress', unique_id, user_id)
    BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    script_dir = os.path.join(BASE_DIR, f'pipeline/project_{tools_id}')
    # print("BASE_DIR: ", BASE_DIR) # /data/webapp/mysite
    # print("script_dir: ", script_dir) # /data/webapp/mysite/pipeline/project_hla
    if tools_id == 'hpa':
        subprocess.run(f'python {script_dir}/hpa.py {upload_dir} {project_name} {upload_dir}/tableOfBloodGroupSystems.hpa.xls', shell=True)  # 调用python脚本
    elif tools_id == 'hla':
        # print("yes hla script is running")
        subprocess.run(f"python {script_dir}/hla.py {upload_dir} {project_name}", shell=True)
    elif tools_id == 'rbc':
        subprocess.run(f'python {script_dir}/rbc.py {upload_dir} {project_name} {upload_dir}/tableOfBloodGroupSystems.rbc.xls', shell=True)

    # 更新任务状态为"打包中" & 打包文件
    update_task_status('packaging', unique_id, user_id)
    full_zip_path = os.path.join(upload_dir, project_name) + '.zip'
    subprocess.run(f'xargs -a {upload_dir}/{project_name}/need_to_be_packaged.txt zip -q -j {full_zip_path}', shell=True, check=True)
    
    # 更新任务状态为"已完成" & print end_time.
    end_time = timezone.now()
    update_task_status('completed', unique_id, user_id, end_time=end_time)