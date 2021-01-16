FROM sirimullalab/ligandnet:v1

COPY requirements.txt ./
RUN pip install -r requirements.txt

COPY app.py ./app.py

ENTRYPOINT [ "python", "app.py" ]
EXPOSE 5000

