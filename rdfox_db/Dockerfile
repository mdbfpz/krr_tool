# First stage: Use the Python image to install dependencies
FROM python:3.12 AS dependencies

# Set the working directory in the Python stage
WORKDIR /krr_tool/rdfox_db

# Install curl and other necessary packages
RUN apt-get update && apt-get install -y curl && rm -rf /var/lib/apt/lists/*

# Copy the requirements file
COPY rdfox_db/requirements.txt .

# Install Python dependencies in the first stage
RUN pip3 install --no-cache-dir -r requirements.txt

# Second stage: Use the official RDFox image
FROM oxfordsemantic/rdfox

# Set the working directory inside the container
WORKDIR /krr_tool/rdfox_db

# Copy curl from the dependencies stage to the RDFox stage
COPY --from=dependencies /usr/bin/curl /usr/bin/curl

# Switch back to the rdfox user for security
USER rdfox

# Create the RDFox home directory
RUN mkdir -p /home/rdfox/.RDFox && chown -R rdfox:rdfox /home/rdfox/.RDFox

# Copy the RDFox license into the container
COPY rdfox_db/src/RDFox.lic /opt/RDFox/RDFox.lic

# Copy the initialization script into the image
COPY rdfox_db/start.sh /home/start.sh

# Copy the rest of the application code
COPY rdfox_db/ .

# Expose the RDFox REST API port
EXPOSE 12110

# Switch back to the rdfox user for security
USER root

# Create the required directory for RDFox server with correct ownership
RUN mkdir -p /var/lib/RDFox && chown -R rdfox:rdfox /var/lib/RDFox

# Make the initialization script executable
RUN chmod +x /home/start.sh

ENTRYPOINT ["/bin/bash", "/home/start.sh"]

USER rdfox