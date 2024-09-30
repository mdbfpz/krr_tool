#!/usr/bin/env bash
/opt/RDFox/RDFox -server-directory /var/lib/RDFox -persistence file -role "$RDFOX_ROLE" -password "${RDFOX_PASSWORD}" init &&
/opt/RDFox/RDFox -server-directory /var/lib/RDFox daemon &&
/opt/RDFox/RDFox dstore create test_store &&
/opt/RDFox/RDFox active test_store &&
/opt/RDFox/RDFox set output out