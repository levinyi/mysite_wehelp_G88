# Generated by Django 4.0.4 on 2023-11-27 10:05

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('tools', '0001_initial'),
    ]

    operations = [
        migrations.AddField(
            model_name='result',
            name='fq_gz_count',
            field=models.SmallIntegerField(default=2, verbose_name='fq.gz count'),
        ),
    ]