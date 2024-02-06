import shutil
from django.shortcuts import get_object_or_404, render, redirect
from django.http import HttpResponse, HttpResponseBadRequest, JsonResponse
from .models import RBCPanel, Tools, Result, NGSDataPath
from django.contrib.auth.decorators import login_required
from django.views.decorators.http import require_http_methods
from django.utils import timezone
import os, json, uuid
from .tasks import main_task
from .utils.pagination import Pagination
import subprocess
from datetime import datetime
from decouple import config

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

DATA_DIR = config('DATA_DIR')

# Create your views here.
def index(request):
    return render(request, 'index.html')


def tools_list(request):
    tools = Tools.objects.all()
    return render(request, 'tools/list_tools.html', {'tools': tools})


def save_json_as_xls(json_data, file_path):
    with open(file_path, "w") as f:
        for item in json_data:
            line = "\t".join(str(value) for value in item.values())
            f.write(line + "\n")


@login_required
def tools_use(request, tools_id):
    if request.method == "GET":
        data_objects = NGSDataPath.objects.filter(data_type=tools_id)
        available_files_list = []
        for file in data_objects:
            file_dict = {
                'data_path': file.data_path,
                'data_name': file.data_name,
                'create_time': file.create_time.strftime('%Y-%m-%d'),  # 假设 create_time 是一个 datetime 对象
                'status': file.status,
                'fq_count': file.fq_count,
            }
            available_files_list.append(file_dict)
        
        if tools_id == 'rbc':
            panel_list = RBCPanel.objects.all()
            return render(request, f'tools/tools_{tools_id}_use.html', {'tools_id': tools_id,'available_files':available_files_list, 'panel_list': panel_list})
        else:
            return render(request, f'tools/tools_{tools_id}_use.html', {'tools_id': tools_id,'available_files':available_files_list})
    elif request.method == "POST":
        data = json.loads(request.body)
        project_name = data.get("projectName")
        selectedPath = data.get("selectedPath")
        rawdata_path = NGSDataPath.objects.get(data_name=selectedPath, data_type=tools_id).data_path

        project_base_dir = os.path.join(DATA_DIR, f"{tools_id.upper()}","analysis")
        project_dir = os.path.join(project_base_dir, project_name)
        if not os.path.exists(project_dir):
            os.mkdir(project_dir)

        software_list = None
        output_site = None

        if tools_id == 'hla':
            software_list = ",".join(data.get('selectedSoftware'))
            output_site = data.get('outputOption')  # 'all_sites' or 'partial_sites'
            print(f"output_site: {output_site}")
        elif tools_id == 'hpa':
            json_data = data.get('tableData')
            # file_extension = 'hpa.xls' if tools_id == 'hpa' else 'rbc.xls'
            file_extension = 'hpa.xls'
            file_path = os.path.join(project_dir, f'tableOfBloodGroupSystems.{file_extension}')
            save_json_as_xls(json_data, file_path)
        elif tools_id == 'rbc':
            panel_id = data.get('panel_id')
            print(f"panel_id: {panel_id}")

            panel = RBCPanel.objects.get(id=panel_id)
            panel_file = panel.panel_file.path
            print(f"panel_file: {panel_file}")

            new_file_path = os.path.join(project_dir, "tableOfBloodGroupSystems.rbc.xls")
            shutil.copy(panel_file, new_file_path)

        unique_id = str(uuid.uuid4())
        user_id = request.user.id

        result = Result.objects.create(
            user = request.user, 
            tools_name = tools_id,
            unique_id = unique_id,
            proj_name = project_name,
            result_path = project_dir,
            raw_data = rawdata_path,
            status = 'pending',
        )
        print(f"submit a {tools_id} task: project name : {project_name}")
        
        # 异步处理：
        try:
            main_task.delay(user_id, tools_id, unique_id, software_list=software_list, output_site=output_site)
        except Exception as e:
            response_data = {
                'status': 'error',
                'message': str(e),
            }
            return JsonResponse(response_data, status=500)
        return JsonResponse({"status": "success", "message": "submit successfully."})
    else:
        return JsonResponse({"status":'error',"response":'error Post Method.'})

@login_required
def check_status(request, tools_id):
    # if request.method == "GET":
    queryset = Result.objects.filter(user=request.user, tools_name=tools_id).order_by("-id")
    page_obj = Pagination(request, queryset, page_size=15, page_param="page", plus=5)
    context = {
        "queryset": page_obj.page_queryset,
        "page_string": page_obj.html(),
        "tools_id": tools_id,
    }
    return render(request, 'tools/check_status.html', context)

@login_required
def delete_status(request, unique_id):
    if request.method == 'POST':
        result = Result.objects.get(user=request.user, unique_id=unique_id)
        # delete project folder according to unique_id
        project_folder = result.result_path
        print("delete project_folder: ", project_folder)
        subprocess.run(f"rm -rf {project_folder}", shell=True)
        # delete result object
        tools_id = result.tools_name
        result.delete()
        return JsonResponse({'status': 'success', 'message': 'Result deleted successfully'})
    else:
        return HttpResponseBadRequest("Invalid request method.")


@login_required
def download_result(request, unique_id):
    try:
        result = Result.objects.get(unique_id=str(unique_id))
        if result.status != "completed":
            return HttpResponse("任务尚未完成，无法下载")

        result_path = result.result_path
        print("result_path: ",result_path)
        if not os.path.exists(result_path):
            return HttpResponse("result path have been deleted by someone.")

        zip_file_name = os.path.basename(result.proj_name) + ".zip"
        full_zip_file = os.path.join(result_path, result.proj_name + "_Analysis_Result.zip")
        print("full_zip_file: ", full_zip_file)

        if not os.path.exists(full_zip_file):
            print("file does not exists!")
            return HttpResponse("file does not exists")

        print("starting open zip file to response!")
        with open(full_zip_file, "rb") as f:
            response = HttpResponse(f)
            response['Content-Type'] = 'application/octet-stream'
            response['Content-Disposition'] = f'attachment;filename="{zip_file_name}"'
            return response
    except Result.DoesNotExist:
        return HttpResponse(f"Result with unique ID {unique_id} not found.")

# 存储路径的字典
rawdata_dirs = {
    'hla': str(os.path.join(DATA_DIR, "HLA", "rawdata")),
    'hpa': str(os.path.join(DATA_DIR, "HPA", "rawdata")),
    'rbc': str(os.path.join(DATA_DIR, "RBC", "rawdata")),
}

@require_http_methods(["POST"])
def update_data_path(request):
    # get data_type
    data = json.loads(request.body)
    data_type = data.get('data_type')  # hla, hpa, rbc from front-end
    if not data_type:
        return JsonResponse({'message':'No data type provided'}, status=400)

    if data_type not in rawdata_dirs:
        return JsonResponse({'message': 'Invalid data type'}, status=400)

    base_path = rawdata_dirs[data_type]
    if not os.path.exists(base_path):
        return JsonResponse({'message': 'Base path does not exist'}, status=400)

    try:
        # Remove existing entries for the given data_type
        NGSDataPath.objects.filter(data_type=data_type).delete()

        # Add new entries
        for subdir in os.listdir(base_path):
            subdir_path = os.path.join(base_path, subdir)
            if os.path.isdir(subdir_path):
                fq_gz_files = [file for file in os.listdir(subdir_path) if file.endswith('fq.gz')]
                fq_gz_count = len(fq_gz_files)
                if fq_gz_count > 0:
                    create_time = datetime.fromtimestamp(os.path.getctime(subdir_path)).date()
                    NGSDataPath.objects.update_or_create(
                        data_path=subdir_path,
                        data_name=subdir,
                        data_type=data_type,
                        defaults={'status': 'active', 'create_time': create_time, 'fq_count': fq_gz_count}
                    )
        return JsonResponse({'status': 'success', 'message': 'Data updated successfully'})
    except Exception as e:
        print("try Exception as e:", e)
        return JsonResponse({'message': str(e)}, status=500)


def download_file(request, full_zip_path):
    '''下载任意全路径的文件都可以'''
    # 在点击时检查参数是否有效
    if not full_zip_path:
        return HttpResponseBadRequest("Invalid file path")
    
    file_name = os.path.basename(full_zip_path)
    with open(full_zip_path, 'rb') as f:
        response = HttpResponse(f)
        response['Content-Disposition'] = f'attachment; filename={file_name}'
        return response

def panel_upload(request):
    if request.method == 'POST':
        panel_name = request.POST.get('panel_name')
        panel_file = request.FILES.get('panel_file')
        print(f"panel_name: {panel_name}, panel_file: {panel_file}")
        if not panel_name or not panel_file:
            return JsonResponse({'status': 'error', 'message': 'Invalid parameters'}, status=400)
        try:
            RBCPanel.objects.create(name=panel_name, panel_file=panel_file)
            return JsonResponse({'status': 'success', 'message': 'Panel uploaded successfully'})
        except Exception as e:
            return JsonResponse({'status': 'error', 'message': str(e)}, status=500)
    else:
        return render(request, 'tools/tools.rbc_use.html')

def panel_delete(request, panel_id):
    if request.method == 'POST':
        panel = RBCPanel.objects.get(id=panel_id)
        file_path = panel.panel_file.path
        panel.delete()
        # 同时删除服务器上的文件
        if os.path.exists(file_path):
            os.remove(file_path)
        return JsonResponse({'status': 'success', 'message': 'Panel deleted successfully'})
    else:
        return render(request, 'tools/tools.rbc_use.html')
    
def test(request):
    return render(request, 'test.html')

