FROM informaticsmatters/rdkit-python3-debian:Release_2020_09
USER ${UID}:${GID}

ENV DEBIAN_FRONTEND noninteractive
RUN apt update --allow-releaseinfo-change && \
    apt upgrade -y && \
    apt install -y perl && \
    rm -rf /var/lib/apt/lists/*

COPY requirements.txt .
RUN pip install -r requirements.txt

WORKDIR /app
COPY models/files ./models/files
COPY ddt ./ddt
COPY ligandnet.py best_models.txt ./

ENTRYPOINT ["python3", "ligandnet.py"]
