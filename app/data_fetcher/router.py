import httpx
from fastapi import FastAPI, HTTPException, Request


class DataFetcher:
    def __init__(self, converter_queue):
        """Initialize the DataFetcher object with the shared queue."""

        # Define expected content types
        self.content_types = {
            "json": "application/json",
            "xml": ["application/xml", "text/xml"]
        }
        self.converter_queue = converter_queue
        self.app = FastAPI()

        # Define the 'home' GET endpoint
        @self.app.get("/")
        async def index():
            return {
                "info": "Use '/receive_data/json' or '/receive_data/xml' endpoint to send data of the corresponding type."
            }

        # Define a single POST endpoint for receiving data in JSON/XML format
        @self.app.post("/receive_data/{data_type}")
        async def receive_data(data_type: str, request: Request):
            """Receive data from a client and add it to the converter queue."""

            try:
                client_ip = request.client.host

                # TODO: implement this later
                """# Check if the client IP is allowed
                if client_ip != os.getenv("ALLOWED_PROXY_IP", "127.0.0.1"):
                    raise HTTPException(status_code=httpx.codes.FORBIDDEN, detail="Unauthorized proxy.")"""
                

                # Check if the data_type is valid
                if data_type not in self.content_types.keys():
                    raise HTTPException(
                        status_code=httpx.codes.UNSUPPORTED_MEDIA_TYPE,
                        detail="Invalid data type specified."
                    )
                
                # Validate the request's Content-Type
                expected_content_type = self.content_types[data_type]
                content_type = request.headers.get("Content-Type")
                if isinstance(expected_content_type, list):
                    if content_type not in expected_content_type:
                        raise HTTPException(
                            status_code=httpx.codes.UNSUPPORTED_MEDIA_TYPE,
                            detail=f"Content type must be one of: {', '.join(expected_content_type)}."
                        )
                else:
                    if content_type != expected_content_type:
                        raise HTTPException(
                            status_code=httpx.codes.UNSUPPORTED_MEDIA_TYPE,
                            detail=f"Content type must be {expected_content_type}."
                        )
                
                # Extract data from the request
                received_data = await (request.json() if data_type == "json" else request.body())
                
                # Add received data to the shared queue
                await self.converter_queue.add_to_queue(
                    {
                        "data_type": data_type,
                        "data": received_data
                    }
                )

                return {"status": "Data received successfully."}
            
            except HTTPException as http_exc:
                # Raise the HTTP exception to return the correct status code
                raise http_exc
            except Exception as e:
                # Handle generic exceptions and return an error response with 400 status
                raise HTTPException(
                    status_code=httpx.codes.BAD_REQUEST,
                    detail=f"Invalid data format: {str(e)}"
                )