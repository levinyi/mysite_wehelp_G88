from django.contrib import admin
from .models import Tools, Result, NGSDataPath, RBCPanel

# Register your models here.
class ToolsAdmin(admin.ModelAdmin):
    list_display = ('tools_name','tools_desc','tools_icon')
    list_filter = ('tools_name',)


class ResultAdmin(admin.ModelAdmin):
    list_display = ('unique_id','proj_name', 'tools_name', 'user','status','created_at','end_time', 'result_path')


class NGSDataPathAdmin(admin.ModelAdmin):
    list_display = ('data_name','data_type', 'fq_count', 'data_path', 'create_time', 'status')
    list_filter = ('data_type', 'data_type', 'create_time', 'status')


class RBCPanelAdmin(admin.ModelAdmin):
    list_display = ('name', 'panel_file', 'created_at', 'updated_at')
    list_filter = ('name', 'created_at', 'updated_at')


admin.site.register(Tools, ToolsAdmin)
admin.site.register(Result, ResultAdmin)
admin.site.register(NGSDataPath,NGSDataPathAdmin)
admin.site.register(RBCPanel, RBCPanelAdmin)
