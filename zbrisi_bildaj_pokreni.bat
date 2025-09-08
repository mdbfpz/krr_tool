@echo off
echo Stopping any running containers...
docker compose down

REM Remove the RDFox v00001 directory if it exists
SET "folder=rdfox_db\src\data\v00001"
IF EXIST "%folder%" (
    echo Removing %folder% ...
    rmdir /s /q "%folder%"
) ELSE (
    echo No v00001 directory found to remove.
)

REM Run docker compose up --build from the project root
echo Running docker compose up --build ...
docker compose up --build
