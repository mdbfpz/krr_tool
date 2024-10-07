# Use the official Python image from the Docker Hub
FROM python:3.10-slim

# Set environment variables to prevent Python from writing .pyc files and to buffer output
ENV PYTHONUNBUFFERED=1
ENV PYTHONDONTWRITEBYTECODE=1

# Set the working directory inside the container
WORKDIR /krr_tool
# Copy the entire krr_tool directory (from the build context) into the container
COPY . .
# Install dependencies for other services
RUN pip install --no-cache-dir -r /krr_tool/data_fetcher/requirements.txt -r /krr_tool/rdf_converter/requirements.txt
RUN pip install --no-cache-dir -r /krr_tool/rdfox_db/requirements.txt -r /krr_tool/rdfox_db/requirements.txt

# Add krr_tool to PYTHONPATH so that data_fetcher can be imported by processor
ENV PYTHONPATH=/krr_tool

# Set the working directory to processor
WORKDIR /krr_tool/processor

# Command to run the processor
CMD ["python", "processor.py"]