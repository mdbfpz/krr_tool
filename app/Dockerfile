FROM python:3.10-slim

ENV PYTHONUNBUFFERED=1
ENV PYTHONDONTWRITEBYTECODE=1
ENV PYTHONPATH=/krr_tool

WORKDIR /krr_tool

# Install rdfox_db dependencies first (cached unless rdfox_db/requirements.txt changes)
COPY rdfox_db/requirements.txt rdfox_db/requirements.txt
RUN pip install --no-cache-dir -r rdfox_db/requirements.txt

# Install app dependencies (cached unless app/requirements.txt changes)
COPY app/requirements.txt app/requirements.txt
RUN pip install --no-cache-dir -r app/requirements.txt

# Copy the rest of the source code
COPY . .

WORKDIR /krr_tool/app

EXPOSE 8000

CMD ["python", "processor/processor.py"]