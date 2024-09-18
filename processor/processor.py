from data_fetcher.data_fetcher import DataFetcher
from data_sender.data_sender import DataSender
from rdf_converter.rdf_converter import RDFConverter
from reverse_rdf_converter.reverse_rdf_converter import ReverseRDFConverter
from rdfox_db.rdfox_db import RDFoxDB
import asyncio


class DataQueue:
    """Manages the queue for storing data asynchronously."""
    
    def __init__(self):
        self.queue = asyncio.Queue()
        self.condition = asyncio.Condition()

    async def add_to_queue(self, data: dict):
        """Add data to the queue asynchronously."""

        async with self.condition:
            await self.queue.put(data)
            self.condition.notify()  # Notify any waiting tasks

    async def get_from_queue(self) -> dict:
        """Retrieve data from the queue asynchronously."""

        async with self.condition:
            while self.queue.empty():
                await self.condition.wait()  # Wait for new data
            return await self.queue.get()

class Processor:
    """Processor microservice that orchestrates the pipeline."""
    
    def __init__(self):
        self.converter_queue = DataQueue()
        self.rdfdb_queue = DataQueue()
        self.reverse_converter_queue = DataQueue()
        self.sender_queue = DataQueue()

        self.data_fetcher = DataFetcher()
        self.rdf_converter = RDFConverter()
        self.rdfox_db = RDFoxDB()
        self.reverse_rdf_converter = ReverseRDFConverter()
        self.data_sender = DataSender()

    async def fetch_and_enqueue(self):
        """Fetch data and add it to the converter queue periodically."""

        while True:
            polaris_data = await self.data_fetcher.run()
            if polaris_data:
                await self.converter_queue.add_to_queue(polaris_data)
            await asyncio.sleep(1)

    async def convert_and_enqueue(self):
        """Convert data to RDF and add it to the RDFox database queue."""

        while True:
            polaris_data = await self.converter_queue.get_from_queue()
            rdf_graph_data = self.rdf_converter.run(polaris_data)
            await self.rdfdb_queue.add_to_queue(rdf_graph_data)

    async def insert_and_enqueue(self):
        """Insert data into the RDFox database asynchronously and add it to the reverse RDF converter queue."""

        while True:
            rdf_graph_data = await self.rdfdb_queue.get_from_queue()
            await self.rdfox_db.insert_graph_data(rdf_graph_data)
            await self.reverse_converter_queue.add_to_queue(rdf_graph_data)

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

        tasks = [
            self.fetch_and_enqueue(),
            self.convert_and_enqueue(),
            self.insert_and_enqueue(),
            self.reverse_convert_and_enqueue(),
            self.send_data()
        ]
        await asyncio.gather(*tasks)

if __name__ == "__main__":
    processor = Processor()
    asyncio.run(processor.run())
