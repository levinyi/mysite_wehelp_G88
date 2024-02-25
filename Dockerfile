FROM python:3.10
LABEL maintainer="dushiyi319@163.com"

# 设置工作目录
WORKDIR /app

# 添加这两行以安装必要的系统依赖
RUN apt-get update && apt-get install -y --no-install-recommends \
    pkg-config \
    default-libmysqlclient-dev \
    default-jdk \
    zip \
    bowtie2 \
    bwa \
    samtools && apt-get clean && rm -rf /var/lib/apt/lists/*

# 验证 软件安装
RUN java -version
RUN bowtie2 --version

COPY . /app/
# 安装项目所需的依赖
RUN pip install -i https://mirrors.aliyun.com/pypi/simple/ --no-cache-dir -r requirements.txt

