from django.shortcuts import render, redirect
from django.http import HttpResponse, JsonResponse
from .models import Tools, Result
from django.contrib.auth.decorators import login_required
import os, json, uuid
from .tasks import main_task
from .utils.pagination import Pagination
import subprocess

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


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
        return render(request, f'tools/tools_{tools_id}_use.html', {'tools_id': tools_id})
    elif request.method == "POST":
        if tools_id == 'hla':
            unique_id = str(uuid.uuid4())
            upload_dir = os.path.join(BASE_DIR, f'pipeline/project_{tools_id}', unique_id)

            user = request.user
            user_id = user.id
            data_path = request.POST.get("data_path")
            project_name = os.path.basename(os.path.dirname(data_path))
            print(data_path, project_name)

            result = Result.objects.create(
                user = user, 
                tools_name = tools_id,
                unique_id = unique_id,
                project_name = project_name,
                result_path = upload_dir,
                data = data_path,
                status = 'pending',
            )
            print("tools_id: ", tools_id)

        elif tools_id == 'hpa':
            # # 获取JSON数据, 并保存为xls文件
            json_data = request.POST.get('jsonData')
            data = json.loads(json_data)
            file_path = os.path.join(upload_dir, 'tableOfBloodGroupSystems.hpa.xls')
            save_json_as_xls(data, file_path)
        elif tools_id == 'rbc':
            json_data = request.POST.get('jsonData')
            data = json.loads(json_data)
            file_path = os.path.join(upload_dir, 'tableOfBloodGroupSystems.rbc.xls')
            save_json_as_xls(data, file_path)

        # 异步处理：
        main_task.delay(user_id, tools_id, unique_id)

        response_data = {
            'check_available': "True",
            'tools_id': tools_id,
            'project_name': project_name,
            'unique_id': unique_id
        }
        return redirect(f'/tools/check-status/{tools_id}/')
    else:
        return JsonResponse({"status":'error',"response":'error Post Method.'})
'''
# tools_use for chuck upload
@login_required
def tools_use(request, tools_id):
    if request.method == "GET":
        return render(request, f'tools/tools_{tools_id}_use.html', {'tools_id': tools_id})
    elif request.method == "POST":
        user = request.user
        project_name = request.POST.get("project_name")
        chunks = int(request.POST.get("chunks"))
        index = int(request.POST.get("index"))

        if index == 0:
            unique_id = str(uuid.uuid4())
        else:
            unique_id = request.POST.get("unique_id")
        
        upload_dir = os.path.join(BASE_DIR, f'pipeline/project_{tools_id}', unique_id)
        os.makedirs(upload_dir, exist_ok=True)

        chunk = request.FILES['chunk']
        with open(os.path.join(upload_dir, str(index)), 'wb+') as destination:
            for chunk in chunk.chunks():
                destination.write(chunk)

        # 如果是最后一个chunk，那么就合并文件
        if index == chunks - 1:
            # 合并文件
            with open(os.path.join(upload_dir, project_name), 'wb+') as destination:
                for i in range(chunks):
                    with open(os.path.join(upload_dir, str(i)), 'rb') as source:
                        destination.write(source.read())
                    # 删除分片文件
                    os.remove(os.path.join(upload_dir, str(i)))
        
            zip_file = os.path.join(upload_dir, project_name)

            user = request.user
            user_id = user.id
            project_name = request.POST.get("project_name")
            result = Result.objects.create(
                user = user, 
                tools_name = tools_id,
                unique_id = unique_id,
                project_name = project_name, 
                result_path = upload_dir,
                data = zip_file,
                status = 'pending',
            )
            print("tools_id: ", tools_id)

            if tools_id == 'hpa':
                # # 获取JSON数据, 并保存为xls文件
                json_data = request.POST.get('jsonData')
                data = json.loads(json_data)
                file_path = os.path.join(upload_dir, 'tableOfBloodGroupSystems.hpa.xls')
                save_json_as_xls(data, file_path)
            elif tools_id == 'rbc':
                json_data = request.POST.get('jsonData')
                data = json.loads(json_data)
                file_path = os.path.join(upload_dir, 'tableOfBloodGroupSystems.rbc.xls')
                save_json_as_xls(data, file_path)

            # 异步处理：
            main_task.delay(user_id, tools_id, unique_id)

            response_data = {
                'check_available': "True",
                'tools_id': tools_id,
                'project_name': project_name,
                'unique_id': unique_id
            }
            return JsonResponse(response_data)  
        else:
            response_data = {
                'check_available': "False",
                'tools_id': tools_id,
                'project_name': project_name,
                'unique_id': unique_id
            }
            return JsonResponse(response_data)     

'''
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
    return redirect(f"/tools/check-status/{tools_id}", tools_id = tools_id)


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
    


def test(request):
    return render(request, 'test.html')
