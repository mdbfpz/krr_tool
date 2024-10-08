import httpx
from rdfox_db.core import RDFoxDB


class RDFoxQuery(RDFoxDB):
    def __init__(self, connection_id):
        """ Inherit selected properties from RDFoxDB and add specific ones for Query. """

        super().__init__()
        self.connection_id = connection_id
        self.headers_output_turtle = self.headers.copy().update({"Accept": "text/turtle; charset=UTF-8"})
        self.headers_output_csv = self.headers.copy().update({"Accept": "text/csv; charset=UTF-8"})

    async def import_turtle_data(self, ttl_data: str):
        """Import Turtle data into the RDFox database asynchronously."""

        async with httpx.AsyncClient() as client:
            try:
                print("Importing here: ", f"{self.endpoint}/datastores/{self.data_store}/content?connection={self.connection_id}")
                import_resp = await client.post(
                    f"{self.endpoint}/datastores/{self.data_store}/content?connection={self.connection_id}",
                    data=ttl_data,
                    headers=self.headers
                )
                if import_resp.status_code in self.valid_responses:
                    print("Turtle data imported successfully.")
                else:
                    print("Turtle data import failed.")
                return

            except httpx.HTTPStatusError as e:
                print(f"HTTP error occurred: {e.response.status_code} - {e.response.text}")
                return
            except httpx.RequestError as e:
                print(f"Request error occurred: {e}")
                return
            except Exception as e:
                print(f"An unexpected error occurred: {e}")
                return
    
    async def select_all(self):
        """Select all facts from the RDFox data store."""

        async with httpx.AsyncClient() as client:
            try:
                # Change query parameters based on select type
                import_resp = await client.get(
                    #f"{self.endpoint}/datastores/{self.data_store}/content?connection={self.connection_id}",
                    #f"{self.endpoint}/datastores/{self.data_store}/content?default",
                    f"{self.endpoint}/datastores/{self.data_store}/content",
                    headers=self.headers_output_turtle
                )
                print("Fetched all from the datastore:\n", import_resp.content)
                if import_resp.status_code in self.valid_responses:
                    print("Query 'select_all' succeeded.")
                else:
                    print("Query 'select_all' failed.")
                return

            except httpx.HTTPStatusError as e:
                print(f"HTTP error occurred: {e.response.status_code} - {e.response.text}")
                return
            except httpx.RequestError as e:
                print(f"Request error occurred: {e}")
                return
            except Exception as e:
                print(f"An unexpected error occurred: {e}")
                return
    