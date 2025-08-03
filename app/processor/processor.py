import asyncio
import logging
import json
import os
from rdfox_db.core import RDFoxDB
from rdfox_db.queries.queries import RDFoxQuery
from app.data_fetcher.router import DataFetcher
from app.data_sender.data_sender import DataSender
from app.rdf_converter.rdf_converter import RDFConverter
from app.reverse_rdf_converter.reverse_rdf_converter import ReverseRDFConverter
from app.modules.conflicts import ConflictDetection
import uvicorn
import json


def load_static_data(): # TODO: make this a class method
    """Load static json data from container"""
    
    path = "/krr_tool/app/data/aixmData.json"
    
    try:
        with open(path, 'r', encoding='utf-8') as f:
            data = json.load(f)
        print(f"Uspješno učitano {len(data)} zapisa iz: {path}")
        return data
    except FileNotFoundError:
        print("JSON datoteka nije pronađena ni u jednoj od očekivanih lokacija:")        

    return None


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

        cd_data = []  # List to store conflict detection data for Javier testing only
        while True:
            fetched_data = await self.converter_queue.get_from_queue()
            # print("Got from queue: ", fetched_data["data_type"], "\n", fetched_data["data"] )
            data_record = fetched_data["data"]
            data_record_type = fetched_data["data_type"]

            # Convert the fetched data to RDF and update the data repo in the meantime
            self.rdf_converter.record_data_to_rdf(data_record, data_record_type)

            # Use repo last state and pass info to the CD module
            timestamp = self.rdf_converter.last_timestamp
            out = self.conflict_detection.detect(self.rdf_converter.data_repository, timestamp)
            # This block is for Javier to work with the CD data
            """if len(next(iter(out.values()))) != 0:
                cd_data.append(out)
            if len(cd_data) == 100:
                with open("/krr_tool/outputs/output_for_javier.json", "w") as f:
                    print("Writing CD data to output_for_javier.json")
                    json.dump({"data": cd_data}, f, indent=2)"""

            # Update the last repo state with the CD findings
            self.rdf_converter.update_cd_findings(self.conflict_detection.detections, timestamp)
            
            # print("Repo: ", self.rdf_converter.data_repository)
            turtle_data = self.rdf_converter.serialize()
            # print("Turtle data: ", turtle_data)
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
        self.rdfox_queries = RDFoxQuery(self.rdfox_db.connection_id)
        
        self.data_fetcher = DataFetcher(self.converter_queue)
        self.rdf_converter = RDFConverter()
        self.conflict_detection = ConflictDetection()
        self.reverse_rdf_converter = ReverseRDFConverter()
        self.data_sender = DataSender()

        self.static_data = load_static_data()
                
        if self.static_data is not None:
            self.rdf_converter._process_aixm_data(self.static_data)
            turtle_data = self.rdf_converter.serialize()
            await self.rdfox_queries.import_turtle_data(turtle_data)
        else:
            print("Statički podaci nisu učitani")
            
        tasks = [
            self.fetch_and_enqueue(),
            self.convert_and_enqueue(),
            self.insert_and_enqueue(),
            # self.reverse_convert_and_enqueue()
            # self.send_data()
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