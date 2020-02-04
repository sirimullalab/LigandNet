FROM hassanmohsin/ligandnet:latest

COPY requirements.txt ./
RUN pip3 install -r requirements.txt

COPY app.py ./app.py

ENTRYPOINT [ "python3", "app.py" ]
EXPOSE 5000