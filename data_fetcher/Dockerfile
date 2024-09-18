# Use the official Python image from the Docker Hub
FROM python:3.10-slim

WORKDIR /krr_tool/data_fetcher
# Copy and install requirements
COPY data_fetcher/requirements.txt .
RUN pip install --no-cache-dir -r ./requirements.txt

# Copy the rest of the application code
COPY data_fetcher/ .