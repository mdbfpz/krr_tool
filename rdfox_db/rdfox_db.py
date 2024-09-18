import os
from rdflib import Graph, Namespace
from rdflib.namespace import RDF, RDFS, XSD
import httpx


EX = Namespace("http://example.org/")

class RDFoxDB:
    def __init__(self):
        """Initialize the RDFDB with RDFox connection details."""

        self.rdfox_endpoint = os.getenv("RDFOX_ENDPOINT", "http://localhost:12110")
        self.rdfox_username = os.getenv("RDFOX_USERNAME", "admin")
        self.rdfox_password = os.getenv("RDFOX_PASSWORD", "password")

    def _initialize_graph(self):
        """Initialize the RDF graph."""
        self.graph = Graph()
        self.graph.bind("rdf", RDF)
        self.graph.bind("rdfs", RDFS)
        self.graph.bind("xsd", XSD)
        self.graph.bind("ex", EX)

    async def insert_graph_data(self, rdf_graph_data):
        """Store RDF graph data in RDFox."""

        if isinstance(rdf_graph_data, Graph):
            self.graph = rdf_graph_data
        else:
            # If not already a Graph object, attempt to parse the input data
            self._initialize_graph()
            self.graph.parse(data=rdf_graph_data, format='turtle')

        turtle_data = self.graph.serialize(format="turtle")

        print("Skipping data insertion...not implemented yet!")

        """# Send RDF data to RDFox asynchronously
        async with httpx.AsyncClient() as client:
            try:
                response = await client.post(
                    f"{self.rdfox_endpoint}/datastore",
                    auth=(self.rdfox_username, self.rdfox_password),
                    headers={"Content-Type": "application/x-turtle"},
                    data=turtle_data
                )
                response.raise_for_status()
                print("Data successfully sent to RDFox.")
            except httpx.HTTPStatusError as e:
                print(f"Failed to store data in RDFox: {e}")"""
