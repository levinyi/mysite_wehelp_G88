from django.db import models
from django.contrib.auth.models import User


# Create your models here.
class Tools(models.Model):
    """小工具管理"""
    tools_name = models.CharField(verbose_name="小工具", max_length=64)
    tools_id = models.CharField(verbose_name="tools_id", max_length=64)
    tools_desc = models.TextField(verbose_name="功能表述")
    tools_icon = models.ImageField(upload_to="static/img/tools_icon", blank=True)


def generate_upload_path(instance, filename):
    # 使用实例对象的 tools_name 和 unique_id 生成动态路径
    upload_path = f'pipeline/project_{instance.tools_name}/{instance.unique_id}/{filename}'
    return upload_path

class Result(models.Model):
    """store result path"""
    user = models.ForeignKey(User, on_delete=models.CASCADE, null=True)
    TOOLS_CHOICES = [
        ('RBC', 'rbc'),
        ('HLA', 'hla'),
        ('HPA', 'hpa'),
    ]
    tools_name = models.CharField(max_length=10, choices=TOOLS_CHOICES)
    unique_id = models.CharField(max_length=255, unique=True)
    project_name = models.CharField(max_length=255,blank=True, null=True)
    result_path = models.CharField(max_length=255,blank=True, null=True)
    data = models.FileField(verbose_name="data", max_length=255, upload_to=generate_upload_path)
    status = models.CharField(max_length=20, default='pending')
    created_at = models.DateTimeField(auto_now_add=True)
    end_time = models.DateTimeField(null=True, blank=True)

    def __str__(self):
        return self.unique_id
