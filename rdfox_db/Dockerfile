# Use the official Python image from the Docker Hub
FROM python:3.10-slim

# Set the working directory inside the container
WORKDIR /krr_tool/rdfox_db
# Install any needed packages specified in requirements.txt
COPY rdfox_db/requirements.txt .
RUN pip install --no-cache-dir -r ./requirements.txt

# Copy the rest of the application code
COPY rdfox_db/ .
