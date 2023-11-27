from django.shortcuts import render, redirect
from django.http import HttpResponse, JsonResponse
from .models import Tools, Result, NGSDataPath
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
        # 写入标题行
        headers = "\t".join(json_data[0].keys())
        f.write(headers + "\n")

        # 写入数据行
        for item in json_data:
            line = "\t".join(str(value) for value in item.values())
            f.write(line + "\n")


@login_required
def tools_use(request, tools_id):
    if request.method == "GET":
        available_files = NGSDataPath.objects.filter(data_type=tools_id)
        return render(request, f'tools/tools_{tools_id}_use.html', {'tools_id': tools_id,'available_files':available_files})
    elif request.method == "POST":
        data = json.loads(request.body)
        project_name = data.get("projectName")
        selectedPath = data.get("selectedPath")
        rawdata_path = NGSDataPath.objects.get(data_name=selectedPath).data_path
        project_base_dir = os.path.join(DATA_DIR, f"{tools_id.upper()}","analysis")
        project_dir = os.path.join(project_base_dir, project_name)
        if not os.path.exists(project_dir):
            os.mkdir(project_dir)

        if tools_id in ['hpa','rbc']:
            json_data = data.get('tableData')
            file_extension = 'hpa.xls' if tools_id == 'hpa' else 'rbc.xls'
            file_path = os.path.join(project_dir, f'tableOfBloodGroupSystems.{file_extension}')
            save_json_as_xls(json_data, file_path)
            print("Saved file:", file_path)

        unique_id = str(uuid.uuid4())
        now = timezone.now()
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
        
        # 异步处理：
        main_task.delay(user_id, tools_id, unique_id)

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
    result = Result.objects.get(user=request.user, unique_id=unique_id)
    # delete project folder according to unique_id
    project_folder = result.result_path
    print("delete project_folder: ", project_folder)
    subprocess.run(f"rm -rf {project_folder}", shell=True)
    # delete result object
    tools_id = result.tools_name
    result.delete()
    return redirect(f"/tools/check_status/{tools_id}", tools_id = tools_id)


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

        zip_file_name = os.path.basename(result.project_name) + ".zip"
        full_zip_file = os.path.join(result_path, result.project_name + ".zip")
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
    data_type = data.get('data_type')
    if not data_type:
        return JsonResponse({'message':'No data type provided'}, status=400)
    # print("your data type is : ", data_type)
    if data_type not in rawdata_dirs:
        return JsonResponse({'message': 'Invalid data type'}, status=400)

    base_path = rawdata_dirs[data_type]

    if not os.path.exists(base_path):
        return JsonResponse({'message': 'Base path does not exist'}, status=400)

    try:
        # print("yes, here")
        for subdir in os.listdir(base_path):
            subdir_path = os.path.join(base_path, subdir)
            if os.path.isdir(subdir_path):
                if any(file.endswith('fq.gz') for file in os.listdir(subdir_path)):
                    create_time = datetime.fromtimestamp(os.path.getctime(subdir_path)).date()
                    NGSDataPath.objects.update_or_create(
                        data_path=subdir_path,
                        data_name=subdir,
                        data_type=data_type,
                        defaults={'status': 'active','create_time': create_time}
                    )
        objects_list = NGSDataPath.objects.filter(data_type=data_type)
        # print(len(objects_list))
        files = list(NGSDataPath.objects.filter(data_type=data_type).values_list('data_name', flat=True))
        print('Update database successfully,and found :',files)
        return JsonResponse({'files':files, 'message':'Data paths updated successfully'})
    except Exception as e:
        print("try Exception as e:", e)
        return JsonResponse({'message': str(e)}, status=500)

def test(request):
    return render(request, 'test.html')
