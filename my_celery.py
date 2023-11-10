import os
from celery import Celery, platforms
from django.conf import settings

# 设置默认的Django设置模块
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'mysite.settings')

# 创建Celery应用.创建了一个名为 app 的 Celery 应用程序对象，并将其命名为 rootpath_PicWall。
app = Celery('mysite')

# 加载配置 # 将从 Django 的 settings 中读取以 CELERY_ 开头的配置项，并将其应用到 Celery 应用程序对象。
app.config_from_object(settings, namespace='CELERY')

# 自动发现任务
# 传递给该方法的参数是包含任务模块的列表。在您的配置中，您指定了 'mysite' 和 'mysite.tools'，以便 Celery 可以自动发现和注册这些模块中的任务。
app.autodiscover_tasks([
    'mysite',
    'mysite.tools',
])
platforms.c_FORCE_ROOT = True