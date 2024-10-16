#!/usr/bin/env bash

# Exit if any command fails
set -e

# Function to check if RDFox is healthy
check_health() {
    local retries=10
    local delay=1
    for attempt in $(seq 1 $retries); do
        # Check the health endpoint
        response=$(curl -s -o /dev/null -w "%{http_code}" "${RDFOX_ENDPOINT}/health")
        
        # Print the status code
        echo "Health check attempt $attempt: HTTP status code $response"

        if [[ "$response" == "200" || "$response" == "204" ]]; then
            echo "RDFox is ready."
            return 0  # RDFox is ready
        fi
        
        echo "RDFox not ready, retrying in $delay seconds..."
        sleep $delay
    done

    echo "RDFox did not become ready in time."
    return 1  # RDFox is not ready
}

# Declare that RDFox setup (init + daemon start) is not complete yet
echo "0" > /var/lib/RDFox/rdf_setup_complete.txt

# Check if RDFox has already been initialized by looking for the 'server.params' file
if [ -f /var/lib/RDFox/server.params ]; then
    echo "RDFox already initialized, skipping initialization."
else
    echo "Initializing RDFox..."
    /opt/RDFox/RDFox -server-directory /var/lib/RDFox -persistence file -role "${RDFOX_ROLE}" -password "${RDFOX_PASSWORD}" init
fi

echo "Starting RDFox daemon..."
/opt/RDFox/RDFox -server-directory /var/lib/RDFox daemon &

# Wait for RDFox to be ready
if check_health; then   
    
    # If everything completed successfully
    echo "1" > /var/lib/RDFox/rdf_setup_complete.txt

else
    echo "Failed to start RDFox, exiting."
    exit 1
fi

# Keep the script running
while true; do
    sleep 60  # Sleep for 60 seconds before checking again (or do nothing)
done