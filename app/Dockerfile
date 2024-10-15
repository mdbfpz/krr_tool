# Use the official Python image from the Docker Hub
FROM python:3.10-slim

# Set environment variables to prevent Python from writing .pyc files and to buffer output
ENV PYTHONUNBUFFERED=1
ENV PYTHONDONTWRITEBYTECODE=1

# Set the working directory inside the container (to the root project folder)
WORKDIR /krr_tool

# Copy the entire krr_tool directory into the container
COPY . .

# Set PYTHONPATH to include the krr_tool directory so all subdirectories are available as modules
ENV PYTHONPATH=/krr_tool

# Install dependencies from rdfox_db
RUN pip install --no-cache-dir -r rdfox_db/requirements.txt

# Set the working directory to app/processor
WORKDIR /krr_tool/app

# Install dependencies for the app service
RUN pip install --no-cache-dir -r requirements.txt

EXPOSE 8000

# Command to run the processor
CMD ["python", "processor/processor.py"]