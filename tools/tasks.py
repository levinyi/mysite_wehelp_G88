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
    print(f"in main task function: {user_id}, {tools_id}, {unique_id}")
    user = User.objects.get(id=user_id)
    result = Result.objects.get(user=user, tools_name=tools_id, unique_id=unique_id)
    project_name = result.proj_name
    data_dir     = result.raw_data 
    result_path  = result.result_path
    print(f"project_name:{project_name},data_dir: {data_dir}, result_path: {result_path}")
    # 更新任务状态为"in_progress" &  执行Python脚本
    update_task_status('running', unique_id, user_id)


    ### 这里要改，以后要将script_dir放到tools/scripts下面
    BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__))) 
    script_dir = os.path.join(BASE_DIR, f'pipeline/project_{tools_id}')
    print("BASE_DIR: ", BASE_DIR) # /data/webapp/mysite
    print("script_dir: ", script_dir) # /data/webapp/mysite/pipeline/project_hla
  
    if tools_id == 'hpa':
        print("YES hpa is running")
        print(f'python {script_dir}/hpa.py {data_dir} {project_name} {result_path}/tableOfBloodGroupSystems.hpa.xls')  # 调用python脚本
        subprocess.run(f'python {script_dir}/hpa.py {data_dir} {project_name} {result_path}/tableOfBloodGroupSystems.hpa.xls', shell=True)  # 调用python脚本
    elif tools_id == 'hla':
        print("yes hla script is running")
        print(f"python {script_dir}/hla.py {data_dir} {project_name} {result_path}")
        subprocess.run(f"python {script_dir}/hla.py {data_dir} {project_name} {result_path}", shell=True)
    elif tools_id == 'rbc':
        subprocess.run(f'python {script_dir}/rbc.py {data_dir} {project_name} {result_path}/tableOfBloodGroupSystems.rbc.xls', shell=True)

    # 更新任务状态为"打包中" & 打包文件
    update_task_status('packaging', unique_id, user_id)
    full_zip_path = os.path.join(result_path, project_name) + '.zip'
    subprocess.run(f'xargs -a {result_path}/{project_name}/need_to_be_packaged.txt zip -q -j {full_zip_path}', shell=True, check=True)
    
    # 更新任务状态为"已完成" & print end_time.
    end_time = timezone.now()
    update_task_status('completed', unique_id, user_id, end_time=end_time)
    print("finished!")

