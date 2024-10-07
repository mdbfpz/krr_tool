import os
import httpx
import base64
import asyncio
import time


class RDFoxDB:
    def __init__(self):
        """Initialize the RDFDB with RDFox connection details."""

        self.endpoint = os.getenv("RDFOX_ENDPOINT", "http://localhost:12110")
        self.role = os.getenv("RDFOX_ROLE", "admin")
        self.password = os.getenv("RDFOX_PASSWORD", "password")
        self.data_store = os.getenv("RDFOX_DSTORE", "test_store")
        self.headers = self._create_import_headers()

        self.setup_completed_file = os.getenv("RDFOX_SETUP_COMPLETED_FILE")
        self.data_store_id = None
        self.connection_id = None
        self.rdfox_ready_flag = False

    def _create_import_headers(self):
        """Create headers for RDFox API requests."""

        headers = {
            "Content-Type": "text/turtle",
            "Authorization": f"Basic {self._encode_credentials()}"
        }

        return headers
    
    def _encode_credentials(self):
        credentials = f"{self.role}:{self.password}"
        return base64.b64encode(credentials.encode()).decode()
    
    
    def _rdfox_setup_completed(self, retries: int = 10, delay: int = 1):
        """Check if the RDFox database is up and running, and data store created."""
        
        # Check if RDFox setup is completed
        for _ in range(retries):
            with open(self.setup_completed_file, "r") as file:
                rdfox_setup_status = file.read().strip()

            if rdfox_setup_status != "1":
                time.sleep(delay)
            else:
                return True
        return False
    
    async def _create_data_store(self, retries: int = 10, delay: int = 1):
        """Create a data store in the RDFox DB."""

        for _ in range(retries):
            if not self._rdfox_setup_completed():
                await asyncio.sleep(delay)
            else:
                async with httpx.AsyncClient() as client:
                    # Create data store
                    try:
                        print("Creating data store...")
                        response = await client.post(
                            f"{self.endpoint}/datastores/{self.data_store}",
                            headers={
                                "Authorization": f"Basic {self._encode_credentials()}"
                            }
                        )

                        print("Create DS response: ", response.headers)
                        dstore_id = response.headers.get("etag").replace('"', "").split("-")[0]
                        dstore_name = response.headers.get("location").split("/")[-1]

                        if dstore_name == self.data_store:
                            self.data_store_id = dstore_id
                            test_resp = await client.get(
                                f"{self.endpoint}/datastores/{self.data_store}/info?component-info=extended",
                                headers={
                                    "Authorization": f"Basic {self._encode_credentials()}"
                                }
                            )

                            print("Dstore id: ", self.data_store_id)
                            print(f"Data store {dstore_name} was successfully created.")

                            await self.bring_dstore_online()

                            return
                        
                        return 

                    except httpx.HTTPStatusError as e:
                        print(f"HTTP error occurred: {e.response.status_code} - {e.response.text}")
                    except httpx.RequestError as e:
                        print(f"Request error occurred: {e}")
                    except Exception as e:
                        print(f"An unexpected error occurred: {e}")

        print("Max retries exceeded, failed to create data store.")
    
    async def bring_dstore_online(self):
        async with httpx.AsyncClient() as client:
            # Bring data store online
            try:
                print("Bringing data store online...")
                response = await client.patch(
                    f"{self.endpoint}/datastores/{self.data_store}?operation=bring-online",
                    headers={
                        "Authorization": f"Basic {self._encode_credentials()}"
                    }
                )
                if response.status_code in [200, 204]:
                    print("Data store brought online successfully.")
                else:
                    print("Not able to bring data store online.")
                # print("Bringing DS online response: ", response.content)
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
    
    def _data_store_created(self):
        return bool(self.data_store_id)

    async def _establish_connection(self, retries: int = 10, delay: int = 1):
        """Establish a connection to the RDFox data store."""

        for _ in range(retries):
            if not self._data_store_created():
                await asyncio.sleep(delay)
            else:
                async with httpx.AsyncClient() as client:
                    # Establish a connection
                    try:
                        print("Establishing data store connection...")
                        print("path: ", f"{self.endpoint}/datastores/{self.data_store}/connections")
                        conn_response = await client.post(
                            f"{self.endpoint}/datastores/{self.data_store}/connections",
                            headers={
                                "Authorization": f"Basic {self._encode_credentials()}"
                            }
                        )
                        print("Conn response headers: ", conn_response.headers)
                        print("Conn response content: ", conn_response.content)
                        conn_id = conn_response.headers.get("location").split("/")[-1]
                        if conn_id:
                            self.connection_id = conn_id
                            print("Connection established.")
                        print("Connection ID: ", self.connection_id)
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

        print("Max retries exceeded, failed to establish connection.")

    async def _conn_established(self, retries: int = 10, delay: int = 1):
        """Check if a connection to the RDFox is established."""

        for _ in range(retries):
            if not self.connection_id:
                await asyncio.sleep(delay)
            else:
                return True

                """async with httpx.AsyncClient() as client:
                    try:                       
                        # Check if the connection was established
                        conn_response = await client.get(
                            f"{self.endpoint}/datastores/{self.data_store}/connections",
                            headers={
                                "Authorization": f"Basic {self._encode_credentials()}"
                            }
                        )
                        if not self.connection_id:
                            self.connection_id = self.decode_conn_response(conn_response)
                        return True

                    except httpx.HTTPStatusError as e:
                        print(f"HTTP error occurred: {e.response.status_code} - {e.response.text}")
                        return False
                    except httpx.RequestError as e:
                        print(f"Request error occurred: {e}")
                        return False
                    except Exception as e:
                        print(f"An unexpected error occurred: {e}")
                        return False"""
        
        print("Max retries exceeded, failed to check connection.")
        return False
    
    def decode_conn_response(self, response):
        conn_response_decoded = response.content.decode("utf-8")
        conn_id = None
        
        # Split the response into lines and extract names
        lines = conn_response_decoded.split("\n")
        if len(lines) > 1:  # Check if there's data after the header
            conn_ids = lines[1:]  # Exclude the first line (header)
            conn_id = [id.strip('"') for id in conn_ids if id][0] # Clean response, use only non-blank IDs without quotes
        
        return conn_id
    
    async def rdfox_ready(self):
        if not await self._conn_established():
            self.rdfox_ready_flag = False
        else:
            self.rdfox_ready_flag = True

        return self.rdfox_ready_flag
    
    async def run(self):
        await self._create_data_store()
        await self._establish_connection()