FROM hassanmohsin/rdkit-openbabel:latest

RUN apt update && apt upgrade -y
RUN apt install -y python3-pip

COPY requirements.txt .
# RUN pip3 install -r requirements.txt

WORKDIR /app
COPY ../ligandnet .

ENTRYPOINT ["python3", "-m", "ligandnet.predict"]

# Web
# COPY app.py .
# ENTRYPOINT [ "python3", "app.py" ]
# EXPOSE 5000