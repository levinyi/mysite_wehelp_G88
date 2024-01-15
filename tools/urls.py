from django.urls import path
from . import views

app_name = "tools"

urlpatterns = [
    path('index/', views.index, name="index"),
    path('tools-list/', views.tools_list, name="tools_list"),
    path('tools-use/<str:tools_id>/', views.tools_use, name='tools_use'),
    path('check_status/<str:tools_id>/', views.check_status, name='check_status'),
    path('download_result/<str:unique_id>/', views.download_result, name='download_result'),
    path('delete_status/<str:unique_id>/', views.delete_status, name='delete_status'),
    path('test/', views.test, name="test"),
    path('update_data_path/', views.update_data_path, name='update_data_path'),
    path('download_file/<path:full_zip_path>/', views.download_file, name='download_file'),
    path('panel_upload/', views.panel_upload, name='panel_upload'),
    path('panel_delete/<int:panel_id>/', views.panel_delete, name='panel_delete'),
]
