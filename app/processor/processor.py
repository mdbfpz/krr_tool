import asyncio
import logging
from rdfox_db.core import RDFoxDB
from rdfox_db.queries.queries import RDFoxQuery
from app.data_fetcher.router import DataFetcher
from app.data_sender.data_sender import DataSender
from app.rdf_converter.rdf_converter import RDFConverter
from app.reverse_rdf_converter.reverse_rdf_converter import ReverseRDFConverter
import uvicorn


class DataQueue:
    """Manages the queue for storing data asynchronously."""
    
    def __init__(self):
        self.queue = asyncio.Queue()
        self.condition = asyncio.Condition()

    async def add_to_queue(self, data: dict): # TODO: rename to 'put'
        """Add data to the queue asynchronously."""

        async with self.condition:
            await self.queue.put(data)
            self.condition.notify()  # Notify any waiting tasks

    async def get_from_queue(self) -> dict:  # TODO: rename to 'get'
        """Retrieve data from the queue asynchronously."""

        async with self.condition:
            while self.queue.empty():
                await self.condition.wait()  # Wait for new data
            return await self.queue.get()

class Processor:
    """Processor class that orchestrates the pipeline."""
    
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
            fetched_data = await self.converter_queue.get_from_queue()
            print("Got from queue: ", fetched_data["data_type"], "\n", fetched_data["data"] )
            data_record = fetched_data["data"]
            data_record_type = fetched_data["data_type"]

            print("Previous timestamp: ", self.rdf_converter.last_timestamp)
            if data_record_type == "json":
                # 1. case - FIXM data
                if "flight_plan" in data_record.keys():
                    self.rdf_converter.convert_fixm_data(data_record)
                    # TODO: remove this maping from here, it should be done in the RDFConverter
                    # TODO: check what happens if there are two same timestamps for xml and json, will they overwrite each other?
                # 2. case - Event data
                else:
                    self.rdf_converter.convert_event_data(data_record)
                    # TODO: check what happens if there are two same timestamps for xml and json, will they overwrite each other?
                
            turtle_data = self.rdf_converter.serialize()
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
