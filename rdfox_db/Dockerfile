# First stage: Use the Python image to install dependencies
FROM python:3.12 AS dependencies

# Set the working directory in the Python stage
WORKDIR /krr_tool/rdfox_db

# Install curl and other necessary packages (you can keep this; not used in final stage)
RUN apt-get update && apt-get install -y curl && rm -rf /var/lib/apt/lists/*

# Copy the requirements file
COPY rdfox_db/requirements.txt .

# Install Python dependencies in the first stage
RUN pip3 install --no-cache-dir -r requirements.txt


# Second stage: Use the official RDFox image
FROM oxfordsemantic/rdfox

# Set the working directory inside the container
WORKDIR /krr_tool/rdfox_db

# --- IMPORTANT FIX: install curl in this stage (do NOT copy from dependencies) ---
USER root
RUN apt-get update && apt-get install -y curl && rm -rf /var/lib/apt/lists/*

# Create RDFox directories and set correct ownership
RUN mkdir -p /home/rdfox/.RDFox /var/lib/RDFox \
    && chown -R rdfox:rdfox /home/rdfox/.RDFox /var/lib/RDFox

# Copy the RDFox license into the container
COPY rdfox_db/src/RDFox.lic /opt/RDFox/RDFox.lic

# Copy the initialization script into the image
COPY rdfox_db/start.sh /home/start.sh
RUN chmod +x /home/start.sh

# Copy the rest of the application code
COPY rdfox_db/ .

# Expose the RDFox REST API port
EXPOSE 12110

# Switch back to the rdfox user for security
USER rdfox

ENTRYPOINT ["/bin/bash", "/home/start.sh"]
