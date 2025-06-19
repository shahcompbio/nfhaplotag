FROM python:3.11-slim
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        procps \
        build-essential \
        libbz2-dev \
        liblzma-dev \
        libz-dev \
        libcurl4-openssl-dev \
        libssl-dev \
        && rm -rf /var/lib/apt/lists/*

RUN pip install --no-cache-dir \
    pysam \
    tqdm \
    pandas \
    seaborn \
    matplotlib \
    numpy \
    cyvcf2 \
    click

CMD ["python3"]