#!/usr/bin/env bash

# Exit if any command fails
set -e

# Declare that RDFox setup is not complete yet
echo "0" > /var/lib/RDFox/rdf_setup_complete.txt

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

echo "Initializing RDFox..."
/opt/RDFox/RDFox -server-directory /var/lib/RDFox -persistence file -role "${RDFOX_ROLE}" -password "${RDFOX_PASSWORD}" init &&

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