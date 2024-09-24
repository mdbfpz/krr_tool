import os
import httpx


class DataFetcher:
    def __init__(self):
        """Initialize the DataFetcher with external API endpoints."""

        # self.lfv_api_endpoint = os.getenv("API_ENDPOINT_LFV", "default_lfv_endpoint")
        self.polaris_api_endpoint = os.getenv("API_ENDPOINT_POLARIS", "default_polaris_endpoint")
        # print(self.polaris_api_endpoint)

       
    async def fetch_from_polaris(self):
        """Fetch data from Polaris API endpoint."""

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
        """Run scheduled tasks."""

        # lfv_data = self.fetch_from_lfv()
        polaris_data = await self.fetch_from_polaris()

        return {
            # "lfv_data": lfv_data,
            "polaris_data": polaris_data
        }
