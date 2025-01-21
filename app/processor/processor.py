import asyncio
import logging
from rdfox_db.core import RDFoxDB
from rdfox_db.queries.queries import RDFoxQuery
from app.data_fetcher.router import DataFetcher
from app.data_sender.data_sender import DataSender
from app.rdf_converter.rdf_converter import RDFConverter
from app.reverse_rdf_converter.reverse_rdf_converter import ReverseRDFConverter
import uvicorn
from .utils import DataQueue, GeodesicService


class Processor:
    """Processor microservice that orchestrates the pipeline."""
    
    def __init__(self):
        self.converter_queue = DataQueue()
        self.rdfdb_queue = DataQueue()
        self.reverse_converter_queue = DataQueue()
        self.sender_queue = DataQueue()

        self.rdfox_db = RDFoxDB()

    async def fetch_and_enqueue(self):
        """
        Run the FastAPI server to listen for incoming POST requests.
        Fetch data and add it to the converter queue when the data is received.
        """

        config = uvicorn.Config(self.data_fetcher.app, host="0.0.0.0", port=8000)
        server = uvicorn.Server(config)
        await server.serve()

    async def convert_and_enqueue(self):
        """Convert data to RDF and add it to the RDFox database queue."""

        while True:
            received_data = await self.converter_queue.get_from_queue()
            print("Got from queue: ", received_data["data_type"], "\n", received_data["data"] )
            turtle_data = self.rdf_converter.run(
                received_data["data"].decode("utf-8") # Convert bytes to string using UTF-8 encoding
            )
            await self.rdfdb_queue.add_to_queue(turtle_data)
    
    async def insert_and_enqueue(self):
        """Insert data into the RDFox database asynchronously and add it to the reverse RDF converter queue."""

        # Wait for RDFox to be ready
        while await self.rdfox_db.rdfox_ready():
            turtle_data = await self.rdfdb_queue.get_from_queue()
            await self.rdfox_queries.import_turtle_data(turtle_data)
            await self.reverse_converter_queue.add_to_queue(turtle_data)
            # Test if /content retrieval works:
            await self.rdfox_queries.select_all()

    async def reverse_convert_and_enqueue(self):
        """Convert RDF data back to JSON and add it to the data sender queue."""

        while True:
            rdf_graph_data = await self.reverse_converter_queue.get_from_queue()
            data_to_send = self.reverse_rdf_converter.run(rdf_graph_data)
            await self.sender_queue.add_to_queue(data_to_send)

    async def send_data(self):
        """Send data back to the external API."""

        while True:
            data_to_send = await self.sender_queue.get_from_queue()
            await self.data_sender.run(data_to_send)

    async def run(self):
        """Run the data processing pipeline."""

        await self.rdfox_db.run()
        
        self.data_fetcher = DataFetcher(self.converter_queue)
        self.rdf_converter = RDFConverter()
        self.reverse_rdf_converter = ReverseRDFConverter()
        self.data_sender = DataSender()
        self.rdfox_queries = RDFoxQuery(self.rdfox_db.connection_id)

        tasks = [
            self.fetch_and_enqueue(),
            self.convert_and_enqueue(),
            self.insert_and_enqueue(),
            self.reverse_convert_and_enqueue()
            #self.send_data()
        ]
        await asyncio.gather(*tasks)

if __name__ == "__main__":
    """
    logging.basicConfig(level=logging.DEBUG)
    httpx_log = logging.getLogger("httpx")
    httpx_log.setLevel(logging.DEBUG)
    """
    processor = Processor()
    asyncio.run(processor.run())
