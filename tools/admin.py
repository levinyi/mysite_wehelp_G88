from django.contrib import admin
from .models import Tools, Result

# Register your models here.
class ToolsAdmin(admin.ModelAdmin):
    list_display = ('tools_name','tools_desc','tools_icon')
    list_filter = ('tools_name',)

admin.site.register(Tools, ToolsAdmin)

class ResultAdmin(admin.ModelAdmin):
    list_display = ('unique_id','project_name','user','status','created_at','end_time', 'result_path')

admin.site.register(Result, ResultAdmin)