"""import os
import httpx


class DataFetcher:
    def __init__(self):
        #Initialize the DataFetcher with external API endpoints.

        # self.lfv_api_endpoint = os.getenv("API_ENDPOINT_LFV", "default_lfv_endpoint")
        self.polaris_api_endpoint = os.getenv("API_ENDPOINT_POLARIS", "default_polaris_endpoint")
        # print(self.polaris_api_endpoint)

       
    async def fetch_from_polaris(self):
        #Fetch data from Polaris API endpoint.

        async with httpx.AsyncClient() as client:
            try:
                response = await client.get(self.polaris_api_endpoint)
                response.raise_for_status()
                data = response.json()
                print("Data fetched from Polaris:", data)
                return data
            except httpx.HTTPStatusError as e:
                print(f"Error fetching data from Polaris: {e}")
                return None

    async def run(self):
        #Run scheduled tasks.

        # lfv_data = self.fetch_from_lfv()
        polaris_data = await self.fetch_from_polaris()

        return {
            # "lfv_data": lfv_data,
            "polaris_data": polaris_data
        }
"""

import os
import httpx
import uvicorn
from fastapi import FastAPI, HTTPException, Request
from processor.processor import DataQueue


class DataFetcher:
    def __init__(self, converter_queue: DataQueue):
        """Initialize the DataFetcher object."""

        self.app = FastAPI()
        self.allowed_proxy_ip = os.getenv("ALLOWED_PROXY_IP", "127.0.0.1")  # Set the allowed proxy IP
        self.converter_queue = converter_queue

        # POST endpoint - receive data from Polaris
        @self.app.post("/receive_data")
        async def receive_data(request: Request):
            """Receive data from a client and add it to the converter queue."""

            client_ip = request.client.host
            if client_ip != self.allowed_proxy_ip:
                raise HTTPException(status_code=httpx.codes.FORBIDDEN, detail="Unauthorized proxy.")
            
            data = await request.json()
            await self.converter_queue.add_to_queue(data)  # Add received data to the shared converter queue
            
            return {"status": "Data received successfully", "data_size": len(data)}

    async def run(self):
        """Run the FastAPI server to listen for incoming POST requests."""

        config = uvicorn.Config(self.app, host="0.0.0.0", port=8000)
        server = uvicorn.Server(config)
        await server.serve()