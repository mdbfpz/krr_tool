# Use the official Python image from the Docker Hub
FROM python:3.10-slim

# Set the working directory inside the container
WORKDIR /krr_tool/rdf_converter
# Install any needed packages specified in requirements.txt
COPY rdf_converter/requirements.txt .
RUN pip install --no-cache-dir -r ./requirements.txt

# Copy the rest of the application code
COPY rdf_converter/ .
