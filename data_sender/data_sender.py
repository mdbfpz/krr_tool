import os
import httpx


class DataSender:
    def __init__(self):
        """Initialize the DataSender with external API endpoints."""

        self.polaris_api_endpoint = os.getenv("API_ENDPOINT_POLARIS", "default_polaris_endpoint")
        # print(f"API Endpoint: {self.api_endpoint}")

    async def send_data(self, data):
        """Send data back to the Polaris API endpoint via POST request."""

        print("Sending data not implemented yet...")
        return

        async with httpx.AsyncClient() as client:
            try:
                response = await client.post(self.api_endpoint, json=data)
                response.raise_for_status()
                print("Data sent back to Polaris API successfully.")
            except httpx.HTTPStatusError as e:
                print(f"Error sending data to API: {e}")

    async def run(self, data):
        """Run scheduled tasks."""

        # Send data back to the API
        if data:
            await self.send_data(data)
