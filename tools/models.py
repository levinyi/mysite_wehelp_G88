from django.db import models
from django.contrib.auth.models import User


# Create your models here.
class Tools(models.Model):
    """小工具管理"""
    tools_name = models.CharField(verbose_name="小工具", max_length=64)
    tools_id = models.CharField(verbose_name="tools_id", max_length=64)
    tools_desc = models.TextField(verbose_name="功能表述")
    tools_icon = models.ImageField(upload_to="static/img/tools_icon", blank=True)



class Result(models.Model):
    """store result path"""
    user = models.ForeignKey(User, on_delete=models.CASCADE, null=True)
    TOOLS_CHOICES = [
        ('RBC', 'rbc'),
        ('HLA', 'hla'),
        ('HPA', 'hpa'),
    ]
    tools_name  = models.CharField(verbose_name="tool name", max_length=10, choices=TOOLS_CHOICES)
    unique_id   = models.CharField(verbose_name="unique id", max_length=255, unique=True)
    proj_name   = models.CharField(verbose_name="proj name", max_length=255, blank=True, null=True)
    result_path = models.CharField(verbose_name="resu path", max_length=255, blank=True, null=True)
    raw_data    = models.CharField(verbose_name="raw data",  max_length=255, blank=True, null=True)
    status      = models.CharField(verbose_name="status",    max_length=20, default='pending')
    created_at  = models.DateTimeField(verbose_name="start time", auto_now_add=True)
    end_time    = models.DateTimeField(verbose_name="finish time", null=True, blank=True)

    def __str__(self):
        return self.unique_id

class NGSDataPath(models.Model):
    """NGS Data Path for HLA,HPA,RBC"""
    data_path = models.CharField(verbose_name="data path", max_length=332)
    data_name = models.CharField(verbose_name="data name", max_length=332)
    data_type = models.CharField(verbose_name="data type", max_length=332, null=True, blank=True)
    fq_count  = models.SmallIntegerField(verbose_name="fq.gz count", default=2)
    create_time = models.DateField(verbose_name="create time", db_index=True)
    status = models.CharField(verbose_name="status", max_length=332)

